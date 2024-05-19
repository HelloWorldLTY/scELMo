"""
Intersect queries with enformer region split
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
..Notes::
* By default the Enformer input sequences are of length 196,608 bp.
 These regions were taken from the Basenji2 work with regions of length
 131,072 bp and extended by 32,768 bp to each side.
 The 131,072 bp sequences were shared by the authors via
 (https://console.cloud.google.com/storage/browser/basenji_barnyard/data).
 By default this script trims the shared sequences to the central
 114,688 bp, because Enformer is only trained to predict over
 those  896 * 128 bp bins of each sequence window.
 The pruning can be disabled via the --no_prune flag. This will intersect
 the TSS with the 131,072 bp sequences.
 Alternatively, using --extend flag the sequence windows can be extended to
 the full 196,608 bp.
* Temporary bed files are temporarily stored in the current working directory
    or the temp_directory specified via --temp-dir

..Input::
1) Query file as produced by create_seq_window_queries.py a raw text
file including a header column that specified the sequence windows to be
processed by the seq model and the positions of the regions of interest
within that sequence to be extracted (roi). Positions and bins of the rois
are comma separated in one string. (.txt or .tsv file). Example format:

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

2) Enformer sequence windows used in the publications. Retrieved from:
https://console.cloud.google.com/storage/browser/basenji_barnyard/data
example format:
chr18   936578  1051266 train
chr4    113639139       113753827       train
chr11   18435912        18550600        train
chr16   85813873        85928561        train
chr3    158394380       158509068       train
chr7    136791743       136906431       train
chr8    132166506       132281194       valid
chr21   35647195        35761883        valid
chr16   24529786        24644474        test
chr8    18655640        18770328        test


..Arguments::
  -h, --help            show this help message and exit
  --query QUERY         A region of interest seq model query file as produced
                        from pre-processing step (1) create_seq_windows.py
  --enf_seqs ENF_SEQS   Path to Enformer sequence coordinate file indicating,
                        chr, start, end and set for each Enformer sequence.
  --query_pos_base QUERY_POS_BASE
                        If the center position in the query file is listed in 1
                        or 0-based format. Default=1
  --out_name OUT_NAME   Name / Path of the prefix for the output files.
  --temp_dir TEMP_DIR   Working directory to store temporary bed files. Default
                        will store in current working directory.
  --no_prune            Flag disable pruning the Enformer windows to the
                        114,688 bp regions over which Enformer actually
                        predicts targets.
  --extend              Flag to enable extending of the shared sequence windows
                        to full Enformer size: 196,608 bp
  --report_full_query   Flag indicate that the full query entries with Enformer
                        sequences intersect should be reported.By default only
                        the gene IDs are saved.
  --strip               Flag indicate to strip version and patch indices from
                        gene ID.

..Usage::

python intersect_queries_with_enformer_regions.py \
--query query_gencode_v41_protein_coding_canonical_tss_hg38_nostitch.tsv \
--enf_seqs sequences.bed \
--strip

..Output::
For each set [train, test, valid], will output a .tsv file specifying which
queries were intersected with which Enformer sequence window as indicated by
the gene ID. Retrun only gene IDs or full queries as controlled by
--retrun_full_query flag.
"""

# Imports ==================================
import argparse
import logging
import subprocess
from datetime import datetime
from os import getcwd, remove
from os.path import join
from typing import Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# define arguments ==============================
parser = argparse.ArgumentParser(
    description="Extract targets over regions of interest from Enformer data."
)
# input
parser.add_argument(
    "--query",
    dest="query",
    type=str,
    required=True,
    help="A region of interest seq model query file as produced from "
    "pre-processing step (1) create_seq_windows.py",
)
parser.add_argument(
    "--enf_seqs",
    type=str,
    required=True,
    help="Path to Enformer sequence coordinate file indicating, chr, start, "
    "end and set for each Enformer sequence.",
)
parser.add_argument(
    "--query_pos_base",
    type=int,
    default=1,
    help="If the center position in the query file is listed in 1 or 0-based "
    "format. Default=1",
)
parser.add_argument(
    "--out_name",
    default="query_enf_intersect",
    type=str,
    required=False,
    help="Name / Path of the prefix for the output files.",
)
parser.add_argument(
    "--temp_dir",
    dest="temp_dir",
    default=".",
    type=str,
    required=False,
    help="Working directory to store temporary bed files. Default will store "
    "in current working directory.",
)
parser.add_argument(
    "--no_prune",
    action="store_true",
    default=False,
    help="Flag disable pruning the Enformer windows to the 114,688 bp regions "
    "over which Enformer actually predicts targets.",
)
parser.add_argument(
    "--extend",
    action="store_true",
    default=False,
    help="Flag to enable extending of the shared sequence windows to "
    "full Enformer size: 196,608 bp",
)
parser.add_argument(
    "--report_full_query",
    action="store_true",
    default=False,
    help="Flag indicate that the full query entries with Enformer sequences "
    "intersect should be reported."
    "By default only the gene IDs are saved.",
)
parser.add_argument(
    "--strip",
    action="store_true",
    default=False,
    help="Flag indicate to strip version and patch indices from gene ID.",
)


# Functions =======================================
def load_enformer_seqs(
    seqs_file: str, prune: bool = True, extend: bool = False
) -> pd.DataFrame:
    """Load, prune and sort Enformer sequence windows

    Arguments
    ---------
    seqs_file: str
        Path to a sequences file. Expects a bed-like format
        with columns: chr, start, stop, set where set indicates which
        set the entry responds to (train, test, valid?). Coordinates are
        expected in bed format as shared by the Basenji2/Enformer authors.
    prune: bool
        If to prune Enformer windows to the central 114,688 bp. Default = True.
    extend: bool
        If extend to full sized Enformer windows of length 196,608 bp.
        Default = False.

    Returns
    -------
    seq_df: pd.DataFrame
        Sequence windows sorted, optionally pruned or extended in padnas
        dataframe.
    """
    seq_df = pd.read_csv(
        seqs_file,
        header=None,
        sep="\t",
        names=["chr", "start", "end", "set"],
    )

    assert seq_df["end"][0] - seq_df["start"][0] == 131_072, (
        "Expects sequences of length 131,072 bp as shared by the "
        "Basenji2/Enformer authors. Stopping ..."
    )

    if prune:
        to_remove = 8_192
        seq_df["start"] += to_remove
        seq_df["end"] -= to_remove

    if extend:
        to_add = 32_768
        seq_df["start"] -= to_add
        seq_df["end"] += to_add

    seq_df = seq_df.sort_values(by=["chr", "start"])

    return seq_df


def save_temp_bed_files(query: pd.DataFrame, seqs: pd.DataFrame) -> None:
    """Save temporary bed files for intersection

    Arguments
    ---------
    query: pd.DataFrame
        Sequence queries loaded as pd.DataFrame.
        Expects: "chr" "start" "end" "patch_id" "query_index"
        "positions_roi" columns
    seqs: str
        Enformer sequences optionally pruned or extended
        stored in pd.Dataframe. Expects "chr" "start" "end" "set" columns.

    Returns
    -------
    Tuple(query_bed_file, seq_bed_file, intersect_file)
        Path to temporary bed files.

    Notes
    -----
    Writes two bed files to disk:
     * the Enformer sequence windows (pruned/extended)
     * a bed file of the sequence query regions of interest
     i.e. TSS as 1 bp long bed regions. Also stores the positions of the
     regions of interest.
    """
    timestamp = str(datetime.now()).replace(" ", "_")

    query_bed_file = join(args.temp_dir, f"temp_query_{timestamp}.bed")
    seq_bed_file = join(args.temp_dir, f"temp_seq_{timestamp}.bed")
    intersect_file = join(args.temp_dir, f"temp_intersect_{timestamp}.bed")

    # make a bed dataframe from the center of the query file
    query_bed = query.loc[:, ["chr", "center", "patch_id", "positions_roi"]]

    if args.query_pos_base == 1:
        query_bed["start"] = query_bed.loc[:, "center"] - 1
        query_bed["end"] = query_bed.loc[:, "center"]
    else:  # pos base 0
        query_bed["start"] = query_bed.loc[:, "center"]
        query_bed["end"] = query_bed.loc[:, "center"] + 1

    query_bed["query_index"] = query_bed.index

    # patch ID is used to distinguish between stitched TSS patches of the
    # same gene
    query_bed = query_bed[
        ["chr", "start", "end", "patch_id", "query_index", "positions_roi"]
    ]

    # save temp bed files for intersect run
    query_bed.to_csv(query_bed_file, header=False, index=False, sep="\t")
    seqs.to_csv(seq_bed_file, header=False, index=False, sep="\t")

    return (query_bed_file, seq_bed_file, intersect_file)


def intersect_query_with_seqs(query: pd.DataFrame, seqs: pd.DataFrame) -> pd.DataFrame:
    """Intersect query centers with Enformer seqs

    Arguments
    ---------
    query: pd.DataFrame
        Sequence queries loaded as pd.DataFrame.
        Expects: "chr" "start" "end" "patch_id" "query_index"
        "positions_roi" columns
    seqs: str
        Enformer sequences optionally pruned or extended
        stored in pd.Dataframe. Expects "chr" "start" "end" "set" columns.

    Returns
    -------
    intersect: pd.DataFrame
        Left out join style intersection between queries and Enformer
        sequences in padnas DataFrame.

    Notes
    -----
     Invokes system bedtools in working directory.
    """
    query_file, seq_file, intersect_file = save_temp_bed_files(query_df, sequences_df)

    intersect_command = (
        f"bedtools intersect -wa -wb -loj -a {query_file} -b "
        f"{seq_file} >{intersect_file}"
    )

    subprocess.call(intersect_command, shell=True, cwd=getcwd())

    intersect = pd.read_csv(
        intersect_file,
        header=None,
        sep="\t",
        names=[
            "chr",
            "start",
            "end",
            "patch_id",
            "query_index",
            "positions_roi",
            "enf_chr",
            "enf_start",
            "enf_end",
            "enf_set",
        ],
    )

    for t in [query_file, seq_file, intersect_file]:
        try:
            remove(t)
        except FileNotFoundError:
            logger.error(
                "Temporary file {t} was missing."
                "Check for error in the intersection and potentially re-run!"
                "Printing first five lines:"
            )
            logger.error(intersect.iloc[:5, :])

    return intersect


def report_missed_queries(num_queries: int, intersect: pd.DataFrame) -> None:
    """Report number of queries not intersecting with Enformer sequence

    Arguments
    ---------
    num_queries: int
        Total number of queries provided.
    intersect: pd.DataFrame
        Left out join style intersection between queries and Enformer
        sequences in padnas DataFrame.

    Notes
    -----
    Only reports number of missing genes/tss/queries to logger info.
    """
    num_non_intersect = len(
        pd.unique(
            intersect.loc[intersect["enf_start"] == -1, ["patch_id", "enf_start"]][
                "patch_id"
            ]
        )
    )

    logger.info(
        f"No intersect found for {num_non_intersect} out of" f" {num_queries} queries."
    )


def select_most_central_intersect_per_query(intersect: pd.DataFrame) -> pd.DataFrame:
    """Select the most central intersect per query

    Arguments
    ---------
    intersect: pd.DataFrame
        Left out join style intersection between queries and Enformer
        sequences in padnas DataFrame.

    Returns
    -------
    intersect: pd.DataFrame
        Left out join style intersection between queries and Enformer
        sequences in padnas DataFrame, but with single entry selected per
        query.

    Notes
    -----
    If a query (ROI/TSS) intersect with multiple Enformer sequences will
    select intersect entry where the query is the most central.
    """
    intersect["enf_center"] = (
        intersect["enf_end"] - intersect["enf_start"] - 1 + intersect["enf_start"]
    )

    intersect["query_enf_distance"] = abs(intersect["enf_center"] - intersect["end"])

    total_intersects = intersect.shape[0]

    intersect = intersect.loc[intersect.groupby("patch_id").query_enf_distance.idxmin()]

    logger.info(
        f"Selected {intersect.shape[0]} of {total_intersects} total " f"intersections."
    )

    return intersect


def split_set(
    df: pd.DataFrame, split_col: str = "enf_set"
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Split a data frame into train/test/valid.

    Arguments
    ---------
    df: pd.DataFrame
        Data frame with at least one column specified by 'split_col' that
        contains 'train', 'valid' and 'test' entries
    split_col: str,
        Name of the column based on which to split the data frame entries.

    Returns
    -------
    Tuple of three dataframes with the train, valid and test
    entries respectively.

    Notes
    -----
    Expects 'train', 'valid' and 'test' entries in the specified split_col.
    """
    test_df = df.loc[df[split_col] == "test"].copy()
    valid_df = df.loc[df[split_col] == "valid"].copy()
    train_df = df.loc[df[split_col] == "train"].copy()

    assert len(train_df) > 0, f"No 'train' entries found in column {split_col}!"
    assert len(test_df) > 0, f"No 'test' entries found in column {split_col}!"
    assert len(valid_df) > 0, f"No 'valid' entries found in column {split_col}!"

    return (train_df, test_df, valid_df)


if __name__ == "__main__":

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # parse no_prune flag to keep pruning as default behavior
    to_prune = True
    if args.no_prune or args.extend:
        to_prune = False

    assert args.query_pos_base in [0, 1], (
        "Only query_pos_base 0 or 1 are " "supported. Stopping ..."
    )

    query_df = pd.read_csv(args.query, sep="\t")

    sequences_df = load_enformer_seqs(args.enf_seqs, to_prune, args.extend)

    intersect_df = intersect_query_with_seqs(query_df, sequences_df)

    report_missed_queries(query_df.shape[0], intersect_df)

    # remove queries without an intersect
    intersect_df = intersect_df.loc[intersect_df["enf_start"] != -1, :]

    # select most central intersect per query
    intersect_df = select_most_central_intersect_per_query(intersect_df)

    train_df, test_df, valid_df = split_set(intersect_df)

    if args.strip:
        # remove version and patch indices
        train_df["patch_id"] = train_df["patch_id"].str.replace(
            "\\.\\w+", "", regex=True
        )
        valid_df["patch_id"] = valid_df["patch_id"].str.replace(
            "\\.\\w+", "", regex=True
        )
        test_df["patch_id"] = test_df["patch_id"].str.replace("\\.\\w+", "", regex=True)

    if args.report_full_query:
        train_df.to_csv(args.out_name + "_train.tsv", sep="\t", header=True)
        valid_df.to_csv(args.out_name + "_valid.tsv", sep="\t", header=True)
        test_df.to_csv(args.out_name + "_test.tsv", sep="\t", header=True)

    else:
        train_df = train_df["patch_id"]
        valid_df = valid_df["patch_id"]
        test_df = test_df["patch_id"]

        train_df.to_csv(
            args.out_name + "_train.txt", sep="\t", header=False, index=False
        )
        valid_df.to_csv(
            args.out_name + "_valid.txt", sep="\t", header=False, index=False
        )
        test_df.to_csv(args.out_name + "_test.txt", sep="\t", header=False, index=False)
