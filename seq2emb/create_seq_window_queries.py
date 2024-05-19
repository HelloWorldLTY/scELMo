
"""
Create embedding query from DNA sequence window and regions of interest.
=========================================
Copyright 2023 GlaxoSmithKline Research & Development Limited. All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=========================================
To compute the embeddings of regions of interest (roi) for sequence windows
with larger sequence inputs (100kb +) it is worth organising them into
queries where each query specifies the sequence window and the roi positions.
This script will take regions of interest, stitch them into patches if desired
and construct sequence windows adhering to chromosome boundaries and create
queries for the sequence model. Genomic position and the id of the prediction
bin with which the rois are intersecting are listed: 0-based! A subsequent
script calculating DNA sequence embeddings can then use the bin_ids to
extract the embeddings of interest.

Stitching: if enabled will group the rois based on a supplied grouping
variable (e.g. a gene name or id) Rois with the same grouping id will be
grouped into patches, Patches are split if they stretch over more than the
supplied threshold (50kb default) and sequence windows are constructed over
the center of patches. The position and bin id of the rois are listed in
a comma separated string.

..Notes::
* rois only accept a single position, if larger regions of interests should
be supplied then please center the coordinate first.

* By default, this script will center the sequence windows on rois or at the
center of a stitched patch. Thus allowing prediction with a maximal sequence
context reach for every roi.

* If the prediction bins of the sequence window are even, then the center of
the sequence window is covered by the border of two bins. In that case the
sequence window is shifted by half a bin size to center the rois within a bin
and avoid overlapping borders as much as possible.

* A single strand per query sequence will be assigned. The strands of the
overlapping regions of interest per query sequence are assessed and the
strand occurring most frequently is used. If different strands occur equally
the first occurring is used. If no strand is provided the + strand is used
by default.

..Input::
    1) a plain text file of region of interest file, where every line
    specifies a region of interest. Supplied via --in. Common formats are bed
    files or vcf file without header.

    chr1	65418	65419	ENST00000641515.2	.	+	ENSG00000186092.7	OR4F5
    chr1	451677	451678	ENST00000426406.4	.	-	ENSG00000284733.2	OR4F29
    chr1	686653	686654	ENST00000332831.5	.	-	ENSG00000284662.2	OR4F16
    chr1	923922	923923	ENST00000616016.5	.	+	ENSG00000187634.13	SAMD11

..Arguments::
  -h, --help            Show this help message and exit
  --in IN_FILE          Input file in bed/vcf or related format. Requires the
                        following columns: chromosome position - 0 or 1 based
                        index specify via --position_base accordingly default=1
                        group_id - an id based on which related regions of
                        interest can be grouped and stitched together, e.g. a
                        gene name or better gene ID if no group_id is supplied
                        the roi will be tagged by chromosome_position and no
                        stitching is performed additional_id - an additional
                        id column which is not used for grouping but added as
                        additional identifier for the patch/roi, e.g. a
                        canonical gene name. If none supplied (set to 0) the
                        group_id will be repeated in the output default=0
  --ref_genome REF_GENOME
                        Path to reference genome fasta file, needs to be indexed
                        and match the coordinate system in the query file.
  --out OUT_QUERY_FILE  Path and name for output file
                        default=./tss_query_file_seq_model.tsv
  --chromosome_col CHROMOSOME_COL
                        Chromosome column default=1.
  --position_col POSITION_COL
                        Position column default=2. Where to find the base pair
                        position of the ROI.Note: the position is treated as
                        1-based by default (change via --position_base)
  --group_id_col group_id_COL
                        Column number of the name/id to use for grouping the
                        regions. default=0. If grouping in patches is turned on
                        this id will be used as basis for grouping. If no
                        grouping id is supplied one   will be created from
                        chrom and position also meaning no grouping/stitching
                        is being performed.
  --additional_id_col additional_id_COL
                        Column number for an additional id column. Will be
                        treated as additional id/name and attached to th query
                        file without processing. E.g. a Gene name or ID.if set
                        to 0 will copy the group id Default=0
  --strand_col STRAND_COL
                        Column number where to find the strand information if
                        present. Set to 0 to ignore strand information. Strand
                        information is used to orient the ROI patch to
                        downstream guide from which strand sequence is
                        extracted. Default will be the + strand unless
                        otherwise specified. Default=0
  --position_base POS_BASE   Position index 0 or 1 based default 1.
  --seq_window_size SEQ_WINDOW_SIZE
                        Size in bp of the DNA sequence window for the desired
                        sequence model (Enformer input as default)
                        default=196608.
  --predict_window_size PREDICT_WINDOW_SIZE
                        Size in bp of the sequence window (central to the
                        sequence window) over which the sequence model
                        actually predicts targets. (Enformer prediction window
                        as default) default=114688.
  --predict_bin_size PREDICT_BIN_SIZE
                        Bin size in bp that the sequence model predicts
                        targets. Default = 128 Used to highlight which
                        prediction bins the rois are intersecting with.
  --threshold STRETCH_THRESHOLD
                        Threshold in bp above which roi patches will be split
                        up default=50000.
  --stitch              Flag to stitch the supplied regions into patches when
                        within the specified threshold
  --no-stitch

..Usage::
python ./create_seq_window_queries.py \
    --in example_files/gexample.bed \
    --chrom_sizes example_files/hg38.chrom.sizes \
    --out ./query_example.tsv \
    --chromosome_col 1\
    --position_col 3\
    --position_base 1 \
    --group_id_col 7\
    --additional_id_col 8\
    --stitch

..Output:: Output is a tab separated query file that lists the chrom start
end, strand of the sequence window the ids of the stitched patch and the
grouping and additional_id, the center of the sequence window the number of
regions of interest within the distance between multiple rois in the sequence
and the strand position and bin id of the rois, comma separated if multiple
ones are available.

chr	seq_start	seq_end	seq_strand	patch_id	group_id	add_id	center
num_roi	stretch	strands_roi	positions_roi	bins_roi
chr1	1	196481	+	OR4F5_0	OR4F5	ENSG00000186092.7	98241	1
0	['+']	65419	191
chr1	353438	549918	-	OR4F29_0	OR4F29	ENSG00000284733.2	451678	1
0	['-']	451678	447
chr1	588414	784894	-	OR4F16_0	OR4F16	ENSG00000284662.2	686654	1
0	['-']	686654	447
chr1	826587	1023067	+	SAMD11_0	SAMD11	ENSG00000187634.13	924827	2
1808	['+', '+']	923923,925731	440,454
chr1	861016	1057496	-	NOC2L_0	NOC2L	ENSG00000188976.11	959256	1
0	['-']	959256	447
chr1	862344	1058824	+	KLHL17_0	KLHL17	ENSG00000187961.15	960584	1
0	['+']	960584	447
chr1	868252	1064732	+	PLEKHN1_0	PLEKHN1	ENSG00000187583.11	966492	2
20 	['+', '+']  966482,966502	447
chr1	883865	1080345	-	PERM1_0	PERM1	ENSG00000187642.10	982105	2
24	 ['-', '-']	982093,982117	447
chr1	901836	1098316	-	HES4_0	HES4	ENSG00000188290.11	1000076	3
191	['-', '-', '-']	999981,1000097,1000172	446,447,448
"""
import argparse
import logging
from sys import exit

import pandas as pd
from pysam import FastaFile

import seq2cells.utils.script_utils as script_utils

parser = argparse.ArgumentParser(
    description="Create sequence model input windows with "
    "regions of interest intersected, "
    "specifying their bin positions , "
    "clustered if desired and with augmentations as desired."
)
parser.add_argument(
    "--in",
    dest="in_file",
    type=str,
    required=True,
    help="Input file in bed/vcf or related format.\n"
    "Requires the following columns:\n"
    "\tchromosome\n"
    "\tposition - 0 or 1 based index specify via "
    "--position_base accordingly (default=1)\n"
    "\tgroup_id - an id based on which related "
    "regions of interest can be "
    "grouped and stitched together, e.g. a gene name or better gene ID "
    "if no group_id is supplied the roi will be tagged by chromosome_position"
    "and no stitching is performed\n"
    "\tadditional_id - an additional id column which is not used for grouping "
    "but added as additional identifier for the patch/roi, "
    "e.g. a canonical gene name. If none supplied (set to 0) "
    "the group_id will be repeated in the output"
    "default=0\n",
)
parser.add_argument(
    "--ref_genome",
    dest="ref_genome",
    type=str,
    required=True,
    help="Path to reference genome fasta file, needs to be indexed and "
    "match the coordinate system in the query file.",
)
parser.add_argument(
    "--out",
    dest="out_query_file",
    default="./query_file_seq_model.tsv",
    type=str,
    required=False,
    help="Path and name for output file default=./tss_query_file_seq_model.tsv",
)
parser.add_argument(
    "--chromosome_col",
    dest="chromosome_col",
    default=1,
    type=int,
    required=False,
    help="Chromosome column default=1.",
)
parser.add_argument(
    "--position_col",
    dest="position_col",
    default=2,
    type=int,
    required=False,
    help="Position column default=2. Where to find the "
    "base pair position of the ROI\n"
    ".Note: the position is treated as "
    "1-based by default (change via --position_base)",
)
parser.add_argument(
    "--group_id_col",
    dest="group_id_col",
    default=0,
    type=int,
    required=False,
    help="Column number of the name/id to use for grouping the regions. "
    "default=0. If grouping in patches is turned on this id will "
    "be used as basis for grouping. If no grouping id is supplied one "
    "will be created from chrom and position also meaning "
    "no grouping/stitching is being performed,",
)
parser.add_argument(
    "--additional_id_col",
    dest="additional_id_col",
    default=0,
    type=int,
    required=False,
    help="Column number for an additional id column. "
    "Will be treated as additional id/name and attached to "
    "th query file without processing. E.g. a Gene name or ID."
    "if set to 0 will copy the group id Default=0 ",
)
parser.add_argument(
    "--strand_col",
    dest="strand_col",
    default=0,
    type=int,
    required=False,
    help="Column number where to find the strand information if present. Set "
    "to 0 to ignore strand information. Strand information is used to "
    "orient the ROI patch to downstream guide from which strand sequence is "
    "extracted. Default will be the + strand unless otherwise specified. "
    "Default=0 ",
)
parser.add_argument(
    "--position_base",
    dest="position_base",
    default=1,
    type=int,
    required=False,
    help="Position index 0 or 1 based default 1.",
)
parser.add_argument(
    "--seq_window_size",
    dest="seq_window_size",
    default=196608,
    type=int,
    required=False,
    help="Size in bp of the DNA sequence window for the desired sequence model "
    "(Enformer input as default) default=196608.",
)
parser.add_argument(
    "--predict_window_size",
    dest="predict_window_size",
    default=114688,
    type=int,
    required=False,
    help="Size in bp of the sequence window (central to the sequence window)"
    " over which the sequence model actually predicts targets. "
    "(Enformer prediction window as default) default=114688.",
)
parser.add_argument(
    "--predict_bin_size",
    dest="predict_bin_size",
    default=128,
    type=int,
    required=False,
    help="Bin size in bp that the sequence model predicts targets. "
    "Default = 128"
    " Used to highlight which prediction bins the ROI's are intersecting with.",
)
parser.add_argument(
    "--threshold",
    dest="stretch_threshold",
    default=50000,
    type=int,
    required=False,
    help="Threshold in bp above which roi patches will be split up. " "default=50000.",
)

# arguments for patching regions or not
parser.add_argument(
    "--stitch",
    dest="to_stitch",
    action="store_true",
    help="Flag to stitch the supplied regions into patches "
    "when within the specified threshold",
)
parser.add_argument("--no-stitch", dest="to_stitch", action="store_false")
parser.set_defaults(to_stitch=False)
parser.add_argument(
    "--debug", dest="debug", action="store_true", help="Flag switch on debugging mode."
)
parser.set_defaults(debug=False)


def store_patches(
    df_roi: pd.DataFrame,
    chromosome_col: int,
    position_col: int,
    position_base: int,
    group_id_col: int,
    additional_id_col: int,
    strand_col: int,
    to_stitch: bool,
):
    """Create a dict and store all roi per grouping id

    Parameters
    ----------
    df_roi: pd.DataFrame
        Dataframe containing the information accessed by the individual
        column identifiers: chromosome, position, group_id and additional_id
    chromosome_col: int
        Column index indicating where the chromosome information is located.
    position_col: int
        Column index indicating where the chromosome information is located.
    position_base: int [0,1]
        Indicating which base the position information is encoded in. 0 or 1
        based.
    group_id_col: int
        Column index indicating which column is to be used as grouping
        variable for patches.
    additional_id_col: int
        Column index indicating which column is to be used as additional
        naming variable not used for grouping the patches.
    strand_col: int
        Column index indicating which column is to be used as strand
        variable for patches
    to_stitch: bool
        Indicate if to perform stitching.

    Returns
    -------
    dict
        Return the a dictionary with regions of interest stored as sequence
        query patches.
    """
    dic_roi = {}
    for index_roi in range(df_roi.shape[0]):
        pos = int(df_roi.iloc[index_roi][position_col])
        if position_base == 0:  # handle 0-based indexed positions
            pos += 1
        chrom = df_roi.iloc[index_roi][chromosome_col]

        # if no group_id supplied use chr_pos and perform no stitching
        if group_id_col < 0:
            group_id = f"{chrom}_{str(pos)}"
        elif group_id_col > 0:
            group_id = df_roi.iloc[index_roi][group_id_col]
        else:
            logger.warning(
                "Please specify a valid column number " "for the " "group_id!"
            )
            exit()

        if additional_id_col < 0:
            # if no additional column supplied reuse group_id
            add_id = group_id
        else:
            add_id = df_roi.iloc[index_roi][additional_id_col]

        patch_counter = 0
        patch = f"{group_id}_{patch_counter}"

        if strand_col < 0:
            strand = "."
        else:
            strand = df_roi.iloc[index_roi][strand_col]

        if patch in dic_roi:
            # skip position if already listed
            if pos not in dic_roi[patch]["pos"]:
                dic_roi[patch]["pos"].append(pos)
                dic_roi[patch]["strands_roi"].append(strand)
            # if stitching suppressed add with new patch name
            if not to_stitch:
                patch_counter += 1
                dic_roi[f"{group_id}_{patch_counter}"] = {
                    "add_id": add_id,
                    "group_id": group_id,
                    "chr": chrom,
                    "pos": [pos],
                    "strands_roi": [strand],
                }
        else:
            # add new entry
            dic_roi[patch] = {
                "add_id": add_id,
                "group_id": group_id,
                "chr": chrom,
                "pos": [pos],
                "strands_roi": [strand],
            }
    return dic_roi


def calc_patch_stats(dic_roi: dict):
    """Calculate start, end, distance, number of roi and center

    Parameters
    ----------
    dic_roi: dict
        Needs to be a dictionary as produced from
        store_patches
        storing chr, start, stop, naming and group ids as well as the
        position of the intersecting regions of interest

    Returns
    -------
    dict
        Return a dictionary with the same contents as input but
        start, end, distance, number of roi and center statistics
        for each patch are added to each entry.
    """
    for k, inner_dict in dic_roi.items():
        inner_dict["num_roi"] = len(inner_dict["pos"])
        inner_dict["start"] = min(inner_dict["pos"])
        inner_dict["end"] = max(inner_dict["pos"])
        inner_dict["stretch"] = inner_dict["end"] - inner_dict["start"]
        inner_dict["center"] = (
            int((inner_dict["end"] - inner_dict["start"]) / 2) + inner_dict["start"]
        )

    return dic_roi


def assign_patch_strand(dic_roi: dict) -> dict:
    """Assign single strand per query patch

    Parameters
    ----------
    dic_roi: dict
        Needs to be a dictionary as produced from
        store_patches -> calc_patch_stats
        storing chr, start, stop, naming and group ids as well as the
        position of the intersecting regions of interest

    Returns
    -------
    dict
        Return a dictionary with the same contents as input but
        a single strand assigned to the query,

    Note
    ----
        Should mostly be the same strand for each ROI patch but decide and warn
        if not uniform.
        If '.' is the most occurring strand it will be set to '+'.
    """
    for key, inner_dict in dic_roi.items():
        temp_strands = inner_dict["strands_roi"]
        temp_strands_str = ",".join(temp_strands)
        strand_decision = script_utils.get_most_common(temp_strands)
        if strand_decision == ".":
            strand_decision = "+"
        if strand_decision not in ["+", "-"]:
            logger.error(
                "Unsupported strand specified please use '+', '-' or "
                f"'.'! Unsupported entry: {key} has strands "
                f"{temp_strands_str}"
            )
            exit()

        if not any([entry == strand_decision for entry in temp_strands]):
            logger.warning(
                "Warning: Found more than one strand in the stitched "
                "regions of interest."
                "Assigning a single strand per patch majority > "
                "first_occurring. '.' will be treated as '+' "
                "Make sure this is your intended behavior."
                f"Here: Supplied strands -> {temp_strands} strand decision -> "
                f"{strand_decision} ..."
            )
        inner_dict["seq_strand"] = strand_decision

    return dic_roi


def split_long_patches(dic_roi: dict, threshold: int) -> dict:
    """Split roi patches stretching more than a threshold

    Parameters
    ----------
    dic_roi: dict
        Needs to be a dictionary as produced from
        store_patches -> calc_patch_stats -> assign_patch_strand
        storing chr, start, stop, naming and group ids as well as the
        position of the intersecting regions of interest
    threshold: int
        Size threshold above which the stretches spanned by intersecting
        regions of interest are split into multiple patches.

    Returns
    -------
    dict
        Return the a dictionary with the same contents as input but the
        patches stretching over threshold split up into separate entries.
        Sequence window queries modified accordingly to center each patch on
        a sequence query.

    Notes
    -----
    Take the entries that span a larger than threshold window ,
    bin it into threshold sized windows
    re-intersect the regions of interest within and save all bins containing
    ROIs as new ROI patches.
    Remove the original roi_patch name thus everything with tag _0
    is not split everything else is split (_1, _2, etc.).
    """
    for g in list(dic_roi):
        if dic_roi[g]["stretch"] > threshold:
            # split
            rang = range(dic_roi[g]["start"], dic_roi[g]["end"], threshold)
            patch_counter = 0  # keep track of produced additional patches
            temp_pos = dic_roi[g]["pos"]
            temp_strand = dic_roi[g]["seq_strand"]
            temp_strands_roi = dic_roi[g]["strands_roi"]
            temp_group_id = dic_roi[g]["group_id"]
            temp_add_id = dic_roi[g]["add_id"]
            temp_chrom = dic_roi[g]["chr"]

            for patch_index in rang:
                temp_start = patch_index
                temp_end = patch_index + threshold
                sub_pos = [
                    x for x in temp_pos if x >= temp_start if x <= temp_end
                ]  # check for roi_in within range
                sub_strand = temp_strands_roi[
                    0 : len(sub_pos)
                ]  # check for roi_in within
                # range
                temp_pos = [
                    x for x in temp_pos if x not in sub_pos
                ]  # set to remaining pos
                temp_strands_roi = temp_strands_roi[
                    len(sub_strand) :
                ]  # set to remaining
                # strands
                # if roi_in within
                if len(sub_pos) > 0:
                    patch_counter += 1
                    # create new patch with intersecting regions
                    patch = dic_roi[g]["group_id"] + "_" + str(patch_counter)
                    temp_min = min(sub_pos)
                    temp_max = max(sub_pos)
                    dic_roi[patch] = {
                        "add_id": temp_add_id,
                        "group_id": temp_group_id,
                        "chr": temp_chrom,
                        "pos": sub_pos,
                        "seq_strand": temp_strand,
                        "strands_roi": sub_strand,
                        "num_roi": len(sub_pos),
                        "start": temp_min,
                        "end": temp_max,
                        "stretch": temp_max - temp_min,
                        "center": int((temp_max - temp_min) / 2) + temp_min,
                    }
            del dic_roi[g]

    return dic_roi


def get_dataframe_from_dic(dic_roi: dict) -> pd.DataFrame:
    """Convert a dictionary of regions of interest to a pandas dataframe

    Parameters
    ----------
    dic_roi: dict
        Needs to be a dictionary as produced from
        store_patches -> calc_patch_stats -> split_long_patches
        storing chr, start, stop, naming and group ids as well as the
        position of the intersecting regions of interest

    Returns
    -------
    pd.DataFrame
        Data frame of the patches of interest containing columns:
        add_id, group_id, patch_id, chr, start, end,
        stretch, num_roi, center, positions_roi
    """
    tss_patches = []
    ids = []
    names = []
    chroms = []
    starts = []
    ends = []
    strands = []
    stretches = []
    nums = []
    centers = []
    positions = []
    strands_roi = []

    for roi_query in dic_roi:
        tss_patches.append(roi_query)
        ids.append(dic_roi[roi_query]["add_id"])
        names.append(dic_roi[roi_query]["group_id"])
        chroms.append(dic_roi[roi_query]["chr"])
        starts.append(dic_roi[roi_query]["start"])
        ends.append(dic_roi[roi_query]["end"])
        strands.append(dic_roi[roi_query]["seq_strand"])
        stretches.append(dic_roi[roi_query]["stretch"])
        nums.append(dic_roi[roi_query]["num_roi"])
        centers.append(dic_roi[roi_query]["center"])
        positions.append(dic_roi[roi_query]["pos"])
        strands_roi.append(dic_roi[roi_query]["strands_roi"])

    query_data_frame = pd.DataFrame(
        {
            "add_id": ids,
            "group_id": names,
            "patch_id": tss_patches,
            "chr": chroms,
            "start": starts,
            "end": ends,
            "seq_strand": strands,
            "stretch": stretches,
            "num_roi": nums,
            "center": centers,
            "positions_roi": positions,
            "strands_roi": strands_roi,
        }
    )
    return query_data_frame


def calc_seq_windows(
    df_roi: pd.DataFrame,
    ref_genome: str,
    sequence_window_size: int,
    predict_window_size: int,
    bin_size: int,
) -> pd.DataFrame:
    """Add sequence window columns to the pandas dataframe.

    Parameters
    ----------
    df_roi: pd.DataFrame
        Needs to be a pandas dataframe as produced from 'get_dataframe_from_dic'
        above with chr and center column at the least.
    ref_genome: str
        Path to reference genome fasta file. used to extract chromosome sizes.
        Needs and .fai index file present.
    sequence_window_size: int
        Integer specifying size of the sequence window to be extracted for the
        sequence model.
    predict_window_size: int
        Size of the window over which the sequence model actually predicts
        targets. (bin_size * number_of_bins)
     bin_size: int
        Size in bp of the prediction bins of the sequence model.

    Returns
    -------
    pd.DataFrame
        Data frame of the patches of interest with additional columns
        'seq_start' 'seq_end' that indicate the DNA sequence model that
        should be queried to cover the patch with regions of interest for
        maximal sequence context. Accounting for chromosome size limits.

    Notes
    -----
        The bins (and whole seq window) is shifted by half a bin size
        if the prediction window size and the bin size of the model lead to
        two prediction bins meeting in the center of a sequence window.
    """
    # init fasta file access
    with FastaFile(ref_genome) as fafile:
        # load all chromosomes
        chromosomes = tuple("chr" + str(i) for i in range(1, 23)) + ("chrX", "chrY")

        # get chrom lengths and assemble in dataframe
        sizes_temp = [fafile.get_reference_length(c) for c in chromosomes]
        chr_sizes = pd.DataFrame({"ref": chromosomes, "length": sizes_temp})

    # calculate extension window
    seq_extend = int(sequence_window_size / 2)
    # two bin being the center
    seq_extend_left = seq_extend
    seq_extend_right = seq_extend

    # if number of bins is even --> shift to avoid the border of
    if (predict_window_size / bin_size) % 2 == 0:
        logger.info("Shifting by half a bin to to the left ...")
        # shift by half a bin to center bin on central TSS'
        seq_extend_left += int(bin_size / 2)
        seq_extend_right -= int(bin_size / 2)

    # 1) extend from center to sequence window
    df_roi["seq_start"] = df_roi["center"] - seq_extend_left
    df_roi["seq_end"] = df_roi["center"] + seq_extend_right - 1

    # 2) shift the ones that are below.chrom start
    temp_offsets = df_roi.loc[df_roi["seq_start"] <= 0, "seq_start"].values * -1 + 1
    df_roi.loc[df_roi["seq_start"] <= 0, "seq_end"] += temp_offsets
    df_roi.loc[df_roi["seq_start"] <= 0, "center"] += temp_offsets
    df_roi.loc[df_roi["seq_start"] <= 0, "seq_start"] += temp_offsets

    # 3) shift left what is over the chromosome end
    for c in chromosomes:
        temp_chr_size = chr_sizes.loc[chr_sizes["ref"] == c]["length"].values[0]
        temp_offsets = (
            df_roi.loc[
                (df_roi["chr"] == c) & (df_roi["seq_end"] > temp_chr_size), "seq_end"
            ].values
            - temp_chr_size
        )
        df_roi.loc[
            (df_roi["chr"] == c) & (df_roi["seq_end"] > temp_chr_size), "seq_start"
        ] -= temp_offsets
        df_roi.loc[
            (df_roi["chr"] == c) & (df_roi["seq_end"] > temp_chr_size), "center"
        ] -= temp_offsets
        df_roi.loc[
            (df_roi["chr"] == c) & (df_roi["seq_end"] > temp_chr_size), "seq_end"
        ] -= temp_offsets

    return df_roi


def calc_seq_roi_bins(
    df_seq_windows: pd.DataFrame, seq_seq_window: int, seq_predict_window: int
):
    """Calculate the indices of the bins in which rois are located.

    Parameters
    ----------
    df_seq_windows: pd.DataFrame
        Comma separated string of bins of interest. Or a single integer.
    seq_seq_window: int
        Integer specifying the size of the DNA window the sequence
        model uses for its predictions.
    seq_predict_window: int
        Integer specifying the size of the window over which the sequence
        model actually predicts output.

    Returns
    -------
    npd.DataFrame
        Data frame of the patches of interest with additional column 'bins_roi'
        indicating the bin indices (0-based) in which the regions of interest
        are located.

    """
    pred_offset = int((seq_seq_window - seq_predict_window) / 2)
    seq_tss = []

    for index_seq_window in range(df_seq_windows.shape[0]):
        seq_s = df_seq_windows.iloc[index_seq_window]["seq_start"]
        tss_ps = df_seq_windows.iloc[index_seq_window]["positions_roi"]
        # correct for 1 based positions --> -1
        tss_bins = [((p - seq_s - pred_offset) - 1) // 128 for p in tss_ps]
        tss_bins_filter = [p for p in tss_bins if p >= 0 if p <= 895]
        seq_tss.append(tss_bins_filter)

    df_seq_windows["bins_roi"] = seq_tss

    df_seq_windows = df_seq_windows[df_seq_windows["bins_roi"].map(len) != 0]

    return df_seq_windows


if __name__ == "__main__":
    # fetch arguments
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # if no group_id column supplied then toggle off stitching
    if args.group_id_col <= 0:
        args.to_stitch = False

    # switch columns to 0-based index for python indexing
    args.chromosome_col -= 1
    args.position_col -= 1
    args.additional_id_col -= 1
    args.group_id_col -= 1
    args.strand_col -= 1

    # calculate number of bins per prediction window
    # and if even toggle on seq window shift
    if args.predict_window_size % args.predict_bin_size != 0:
        logger.warning(
            "Warning: The specified prediction window is not evenly "
            "divisible by the specified bin size."
        )

    # check that sensible position_base was supplied
    assert args.position_base in [0, 1], (
        f"position_base not in [0, 1] got " f"{args.position_base}"
    )

    # LOAD FILES  =======================
    # load roi_in file
    roi_in = pd.read_csv(args.in_file, sep="\t", header=None, comment="#")

    # Process ===========================
    # create a group id from chrom and pos if none specified
    if args.group_id_col <= 0:
        if args.additional_id_col > 0:
            roi_in = roi_in.iloc[
                :, [args.chromosome_col, args.position_col, args.additional_id_col]
            ]
        else:
            roi_in = roi_in.iloc[:, [args.chromosome_col, args.position_col]]
        # add group_id column
        args.group_id_col = roi_in.shape[1]
        roi_in[args.group_id_col] = (
            roi_in.iloc[:, args.chromosome_col].astype(str)
            + "_"
            + roi_in.iloc[:, args.position_col].astype(str)
        )

    # make sure it's sorted by gene_id
    roi_in = roi_in.sort_values(
        [args.group_id_col, args.chromosome_col, args.position_col]
    )
    # store roi_in in dictionary
    roi_dict = store_patches(
        roi_in,
        chromosome_col=args.chromosome_col,
        position_col=args.position_col,
        group_id_col=args.group_id_col,
        additional_id_col=args.additional_id_col,
        strand_col=args.strand_col,
        position_base=args.position_base,
        to_stitch=args.to_stitch,
    )  # store the roi_in in a roi_in patch dictionary
    roi_dict = calc_patch_stats(roi_dict)  # calc stats : num of roi_in ,
    roi_dict = assign_patch_strand(roi_dict)  # assign single strand per patch
    # stretch etc.
    roi_dict = split_long_patches(
        roi_dict, args.stretch_threshold
    )  # split up roi_patches above 50 kb in stretch
    # convert dictionary to pandas dataframe
    roi_df = get_dataframe_from_dic(roi_dict)

    # add sequence model sequence windows accounting for chromosome sizes
    roi_df = calc_seq_windows(
        roi_df,
        ref_genome=args.ref_genome,
        sequence_window_size=args.seq_window_size,
        predict_window_size=args.predict_window_size,
        bin_size=args.predict_bin_size,
    )
    # convert the genomic roi positions to relative positions
    # of the sequence model output bins
    # this will NOT mirror bins for minus strand regions
    roi_df = calc_seq_roi_bins(
        roi_df,
        seq_seq_window=args.seq_window_size,
        seq_predict_window=args.predict_window_size,
    )
    # sort the data frame before saving
    roi_df = roi_df.sort_values(["chr", "seq_start"]).reset_index(drop=True)

    # select and arrange columns
    roi_df = roi_df[
        [
            "chr",
            "seq_start",
            "seq_end",
            "seq_strand",
            "patch_id",
            "group_id",
            "add_id",
            "center",
            "num_roi",
            "stretch",
            "positions_roi",
            "bins_roi",
            "strands_roi",
        ]
    ]

    # convert list of positions to comma separated string
    # collapse duplicated ones
    list_of_rois = []
    for i in range(roi_df.shape[0]):
        roi_positions = roi_df.iloc[i]["positions_roi"]
        roi_positions = script_utils.unique(roi_positions)
        roi_string = ",".join([str(e) for e in roi_positions])
        list_of_rois.append(roi_string)
    roi_df = roi_df.drop(["positions_roi"], axis=1)
    roi_df["positions_roi"] = list_of_rois

    # convert list of bins_roi to comma separated string
    # collapse duplicated ones
    list_of_roi_bins = []
    for i in range(roi_df.shape[0]):
        roi_bin_positions = roi_df.iloc[i]["bins_roi"]
        roi_bin_positions = script_utils.unique(roi_bin_positions)
        roi_bin_string = ",".join([str(e) for e in roi_bin_positions])
        list_of_roi_bins.append(roi_bin_string)
    roi_df = roi_df.drop(["bins_roi"], axis=1)
    roi_df["bins_roi"] = list_of_roi_bins

    # save file (explicit stats file)
    roi_df.to_csv(args.out_query_file, sep="\t", index=False)

    logger.info("Wrote results to " + args.out_query_file)
