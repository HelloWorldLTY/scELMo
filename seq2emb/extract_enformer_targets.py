"""
Extract Enformer targets over regions of interest.
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
1) Read queries and Enformer sequence windows
2) Create temporary bed files for the intersection
3) Perform the intersection with pybedtools
4) Clean up intersect
5) Extract the enformer targets over qery regions per set
6) Save queries and targets

* Enformer TFRecords have been restored as hdf5 files using
    (./research/notebooks/store_enformer_data_pytorch_friendly.ipynb)
* sequence length 114688 because Enformer only predicts over the central
    896 * 128bp bins of each sequence window
* temporary bed files are temporarily stored in the current working directory
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

2) HDF5 (.h5) file containing the Enformer targets in the following structure:
    * ['targets_test']
    * ['targets_valid']
    * ['targets_test']

3) Enformer sequence windows used in the publications. Retrieved from:
https://console.cloud.google.com/storage/browser/basenji_barnyard/data
BUT trimmed to the central 114688 bp per region over which Enformer actually
predicts output and stored in pandas df friendly format.
Supply: train, test and valid sequences
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
  --enf_train_coords ENF_TRAIN_COORDS
                        Path to Enformer sequence coordinate file indicating,
                        chr, start, end and set for the training set for each
                        Enformer sequence. Their length should match the
                        number of prediction bins * bin size.
  --enf_valid_coords ENF_VALID_COORDS
                        Path to Enformer sequence coordinate file indicating,
                        chr, start, end and set for the validation set for
                        each Enformer sequence. Their length should match the
                        number of prediction bins * bin size.
  --enf_test_coords ENF_TEST_COORDS
                        Path to Enformer sequence coordinate file indicating,
                        chr, start, end and set for the test set for each
                        Enformer sequence. Their length should match the
                        number of prediction bins * bin size.
  --enf_data ENF_DATA   Path to Enformer data stored as .h5 file. Requires the
                        keys ['targets_test' , 'targets_valid',
                        'targets_test'] to be present.
  --targets TARGETS_IN  String listing the indices and ranges of indices (all
                        0-based) indicating which targets to select from the
                        predicted output. Only relevant if target saving if
                        toggled on. Ranges are split and the last index per
                        range is extracted as well. E.g "4675:5312" referring
                        to all Enformer CAGE-seq target ids.
  --out_name OUT_NAME   Prefix for output files.
                        Default="out_enf_targets_extracted"
  --store_text          Flag to indicating to store extracted targets as tab
                        separated (tsv) file(s).
  --gz                  Flag indicating that text/tsv file should be gzipped.
  --store_npy           Flag to indicating to store extracted targets as .npy
                        file(s).
  --store_h5            Flag to indicating to store extracted targets as h5py
                        file(s)
  --predict_bin_size PREDICT_BIN_SIZE
                        Bin size in bp that the sequence model predicts
                        targets. Default = 128
  --tar_aggr TAR_AGGR   Select an aggregation method to summarize multiple
                        target outputs of interest per query if present. E.g.
                        multiple TSS. Valid values: 'mean', 'median', 'sum',
                        'min', max'. default='sum
  --add_bins ADD_BINS   Number of neighbouring bins to extract surrounding the
                        bins of interest. Will extract +- <add_bins> bin
                        windows and aggregate them using the selected
                        function. default=0
  --temp_dir TEMP_DIR   Working directory to store temporary bed files.
                        Default will store in current working directory.
  --debug               Flag switch on debugging mode.

..Usage::
python ${SCRIPT_DIR}/extract_enformer_targets.py \
  --query ${QUERY} \
  --enf_train_coords ${TRAIN_COORDS} \
  --enf_valid_coords ${VALID_COORDS} \
  --enf_test_coords ${TEST_COORDS} \
  --enf_data ${ENF_DATA} \
  --out_name "${OUT_DIR}/enf_targets_extracted_all_unique_tss" \
  --targets '0:5,12' \
  --tar_aggr 'sum' \
  --add_bins 0 \
  --store_text \
  --store_npy \
  --store_h5

..Output::
For each set [train, test, valid], will output a .tsv file specifying which
queries were intersected with which Enformer sequence window. Importantly the
"query_index" match the index of the query dataframe (tsv) going into this
script. So it may be used to retrieve the embeddings and predicted targets
produced by `calc_embeddings_and_targets.py`. Example format:

chr     start   end     patch_id        query_index     positions_roi   enf_chr
enf_start       enf_end enf_set enf_set_index   enf_center      query_enf_distance
chr1    924827  924828  ENSG00000187634.13_0    3       923924,925732   chr1
856575  971264  train   11833   971263  46435
chr1    959256  959257  ENSG00000188976.11_0    4       959257  chr1    856575
971264  train   11833   971263  12006
chr1    960584  960585  ENSG00000187961.15_0    5       960585  chr1    856575
971264  train   11833   971263  10678

Depending on the storage options specified the extracted targets are stored
in a tab seperated .tsv file, stored as numpy array .npy file or as hdf5
file. The .tsv and .npy files are stored one for each set. The hdf5 file is
stored in one file with the keys: ['targets_test', 'targets_valid',
'targets_test']
"""


# Imports ==================================
import argparse
import logging
from datetime import datetime
from os import remove
from os.path import exists, join
from typing import Tuple, Union

import h5py
import numpy as np
import pandas as pd
import pybedtools

from seq2cells.utils.data_utils import aggregate_ndarray, get_roi_bins
from seq2cells.utils.script_utils import split_index_range_mix

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
    "--enf_train_coords",
    dest="enf_train_coords",
    type=str,
    required=True,
    help="Path to Enformer sequence coordinate file indicating, chr, start, "
    "end and set for the training set for each Enformer sequence. Their "
    "length should match the number of prediction bins * bin size.",
)
parser.add_argument(
    "--enf_valid_coords",
    dest="enf_valid_coords",
    type=str,
    required=True,
    help="Path to Enformer sequence coordinate file indicating, chr, start, "
    "end and set for the validation set for each Enformer sequence. "
    "Their length should match the number of prediction bins * bin size.",
)
parser.add_argument(
    "--enf_test_coords",
    dest="enf_test_coords",
    type=str,
    required=True,
    help="Path to Enformer sequence coordinate file indicating, chr, start, "
    "end and set for the test set for each Enformer sequence. Their length "
    "should match the number of prediction bins * bin size.",
)
parser.add_argument(
    "--enf_data",
    dest="enf_data",
    type=str,
    required=True,
    help="Path to Enformer data stored as .h5 file. Requires the keys ["
    "'targets_test' , 'targets_valid', 'targets_test'] to be present.",
)
parser.add_argument(
    "--targets",
    dest="targets_in",
    type=str,
    required=True,
    help="String listing the indices and ranges of indices (all 0-based) "
    "indicating which targets to select from the predicted output. Only "
    "relevant if target saving if toggled on. Ranges are "
    "split and the last index per range is extracted as well. "
    'E.g "4675:5312" referring to all Enformer CAGE-seq target ids.',
)
# output
parser.add_argument(
    "--out_name",
    dest="out_name",
    type=str,
    default="out_enf_targets_extracted",
    help='Prefix for output files. Default="out_enf_targets_extracted"',
)
# define what datatypes to store
parser.add_argument(
    "--store_text",
    dest="store_text",
    action="store_true",
    default=False,
    help="Flag to indicating to store extracted targets as tab separated (tsv) "
    "file(s).",
)
parser.add_argument(
    "--gz",
    dest="gz",
    action="store_true",
    default=False,
    help="Flag indicating that text/tsv file should be gzipped.",
)
parser.add_argument(
    "--store_npy",
    dest="store_npy",
    action="store_true",
    default=False,
    help="Flag to indicating to store extracted targets as .npy file(s).",
)
parser.add_argument(
    "--store_h5",
    dest="store_h5",
    action="store_true",
    default=False,
    help="Flag to indicating to store extracted targets as h5py file(s)",
)
parser.add_argument(
    "--predict_bin_size",
    dest="predict_bin_size",
    default=128,
    type=int,
    required=False,
    help="Bin size in bp that the sequence model predicts targets. Default = 128",
)
parser.add_argument(
    "--tar_aggr",
    dest="tar_aggr",
    default="sum",
    type=str,
    required=False,
    help="Select an aggregation method to summarize multiple target "
    "outputs of interest per query if present. E.g. multiple TSS. "
    "Valid values: 'mean', 'median', 'sum', 'min', max'. default='sum",
)
parser.add_argument(
    "--add_bins",
    dest="add_bins",
    default=0,
    type=int,
    required=False,
    help="Number of neighbouring bins to extract surrounding the bins of "
    "interest. Will extract +- <add_bins> bin windows and aggregate them "
    "using the selected function. default=0",
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
    "--debug",
    dest="debug",
    action="store_true",
    default=False,
    help="Flag switch on debugging mode.",
)


def extract_enf_targets_by_intersect(
    intersect_query: pd.DataFrame,
    enf_data_file: str,
    targets: Union[str, int],
    filter: str,
    bin_size: int,
    aggregate: str,
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Split roi patches stretching more than a threshold

    Parameters
    ----------
    intersect_query: pd.DataFrame
        Dataframe of the intersection of Enformer sequence windows with the
        query region data frame. The following columns are required.
        'enf_set' 'query_index' 'enf_start'
        'enf_end' 'positions_roi' 'enf_set_index'
    enf_data_file: str
        Path to .h5 file from where the Enformer data can be read.
        Requires the following keys present:
        ['targets_test' , 'targets_valid', 'targets_test']
    targets: Union[str, int]
        Comma separated string of bins of interest or range of bins.
        Unlike python ranges the last entry of the range will be used as well.
        Or a single integer.
    filter: str ['train', 'test' , 'valid']
        String used to filter the enf seq dataset into train test and
        validation set.
    bin_size: int
        Bin size in bp of the prediction bins.
    aggregate: str ['sum', 'mean', 'median', 'max', 'min']
        Str indicating which aggregation function to use to
        aggregate the targets of multiple ROI per query.

    Returns
    -------
    intersect_query: pd.DataFrame
        A pandas data frame indicating which queries have been intersected
        with enformer windows and which enformer windows have been matched
        and targets extracted from.
    temp_array: np.ndarray
        Numpy array of shape [num_queries, num_targets] containing the
        extracted Enformer targets.

    Notes
    -----

    """
    assert bin_size > 0, "The bin_size needs to be a positive integer!"
    assert filter in ["train", "test", "valid"], (
        "Please choose filter from " "['train', 'test', 'valid']"
    )
    assert aggregate in ["sum", "mean", "median", "max", "min"], (
        "Aggregation method needs to be valid: ['sum', 'mean', 'median', "
        "'max', 'min'] "
    )

    # split up the desired targets to list
    targets_list = split_index_range_mix(targets)

    # set a handle to filter train test valid or not
    if filter:
        intersect_query = intersect_query.loc[intersect_query["enf_set"] == filter]

    # sort by enf_index
    intersect_query = intersect_query.sort_values(by=["query_index"])

    # get length of query and targets to init empty numpy array
    # for storing the extracted data
    temp_array = np.zeros([intersect_query.shape[0], len(targets_list)])

    # get max bins covered (which is the nuber of target bins)
    max_bin = int(
        (intersect_query.iloc[0]["enf_end"] - intersect_query.iloc[0]["enf_start"] - 1)
        / 128
    )

    with h5py.File(enf_data_file, "r") as f:

        extract_counter = 0
        for i in intersect_query.index:
            # get roi positions as integers
            pos_to_get = get_roi_bins(
                bins_in=intersect_query["positions_roi"][i], add_bins=args.add_bins
            )
            # get the Enformer window start position correacted
            # from bed 0-based back to 1-based
            enf_start = intersect_query["enf_start"][i] + 1

            # convert to bins
            bins_to_get = [(p - enf_start) // bin_size for p in pos_to_get]
            # if bins are over 895 --> they extend past the enf seq window
            # and no better match could be found there must be some
            # intersecting bins otherwise it would have been filtered earlier
            # so get rid of all bins over 865
            bins_to_get = [n for n in bins_to_get if n < max_bin]

            # get enf seq window index
            # read it from h5 file
            enf_extract = f[f"targets_{filter}"][intersect_query["enf_set_index"][i]]

            # select the query bins
            enf_extract = enf_extract[bins_to_get, :][:, targets_list]

            # aggregrate using the designated method along axis = 0
            enf_extract = aggregate_ndarray(enf_extract, axis=0, method=aggregate)
            # write to numpy array
            temp_array[extract_counter, :] = enf_extract

            extract_counter += 1  # count up index

            if extract_counter % 500 == 0:
                logger.info(f"{extract_counter} queries extracted ...")

    return intersect_query.reset_index(drop=True), temp_array


if __name__ == "__main__":

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.INFO)

    # asserts ========================
    # check if hdf5 file check_is_not_none
    # stop if so, do not overwrite!
    if args.store_h5:
        out_file_h5 = f"{args.out_name}_targets.h5"
        assert not exists(out_file_h5), (
            "The specified HDF (.h5) file already check_is_not_none. "
            "This script will not overwrite it. Stopping ..."
        )

    # 1) Read queries and Enformer sequence windows ==================
    query_df = pd.read_csv(args.query, sep="\t")
    test_df = pd.read_csv(
        args.enf_test_coords,
        header=None,
        sep="\t",
        names=["chr", "start", "end", "set"],
    )
    valid_df = pd.read_csv(
        args.enf_valid_coords,
        header=None,
        sep="\t",
        names=["chr", "start", "end", "set"],
    )
    train_df = pd.read_csv(
        args.enf_train_coords,
        header=None,
        sep="\t",
        names=["chr", "start", "end", "set"],
    )
    # concat and sort
    seq_df = pd.concat([train_df, test_df, valid_df])
    seq_df = seq_df.sort_values(by=["chr", "start"])

    # 2) Create temporary bed files for the intersection =================
    # get a unique time stamp to avoid mixtures
    timestamp = str(datetime.now())
    query_bed_file = join(args.temp_dir, f"temp_query_{timestamp}.bed")
    seq_bed_file = join(args.temp_dir, f"temp_seq_{timestamp}.bed")
    intersect_file = join(args.temp_dir, f"temp_intersect_{timestamp}.bed")
    # make a bed dataframe from the center
    query_bed = query_df.loc[:, ["chr", "center", "patch_id", "positions_roi"]]
    query_bed["start"] = query_bed.loc[:, "center"] - 1
    query_bed["query_index"] = query_bed.index
    query_bed = query_bed[
        ["chr", "start", "center", "patch_id", "query_index", "positions_roi"]
    ]

    seq_bed = seq_df
    seq_bed["bed_start"] = seq_bed["start"] - 1
    seq_bed["set_index"] = seq_bed.index
    seq_bed = seq_bed[["chr", "bed_start", "end", "set", "set_index"]]

    # 3) Perform the intersection with pybedtools =======================
    # save temp bed files
    query_bed.to_csv(query_bed_file, header=False, index=False, sep="\t")
    seq_bed.to_csv(seq_bed_file, header=False, index=False, sep="\t")
    # intersect the two
    queries = pybedtools.BedTool(query_bed_file)
    seqs = pybedtools.BedTool(seq_bed_file)
    queries.intersect(seqs, wa=True, wb=True, loj=True).saveas(intersect_file)
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
            "enf_set_index",
        ],
    )

    if exists(query_bed_file):
        remove(query_bed_file)
    else:
        logger.error(
            "The query bed file seems to have not been written! "
            "Check for error and re-run!"
        )
    if exists(seq_bed_file):
        remove(seq_bed_file)
    else:
        logger.error(
            "The Enformer sequence bed file seems to have not been "
            "written! Check for error and re-run!"
        )
    if exists(intersect_file):
        remove(intersect_file)
    else:
        logger.error(
            "The intersection file seems to have not been written! "
            "Check for error and re-run!"
        )

    # 3) Count and inform about TSS / genes missed in the intersect =======
    num_non_intersect = len(
        pd.unique(
            intersect.loc[intersect["enf_start"] == -1, ["patch_id", "enf_start"]][
                "patch_id"
            ]
        )
    )
    queries_total = query_df.shape[0]
    logger.warning(
        f"No intersect found for {num_non_intersect} out of"
        f" {queries_total} queries."
    )

    # 4) Clean up intersect =============================
    # remove queries without and
    intersect = intersect.loc[intersect["enf_start"] != -1, :]
    # select most central intersect per query
    intersect["enf_center"] = (
        intersect["enf_end"] - intersect["enf_start"] - 1 + intersect["enf_start"]
    )
    intersect["query_enf_distance"] = abs(intersect["enf_center"] - intersect["end"])
    intersect = intersect.loc[intersect.groupby("patch_id").query_enf_distance.idxmin()]

    # 5) Extract the Enformer targets over query regions per set =========
    # 6) Save queries and targets =======================================
    # looping over test, valid and train set
    for set in ["test", "valid", "train"]:

        logger.info(f"Extracting {set} targets ...")
        out_queries, out_targets = extract_enf_targets_by_intersect(
            intersect,
            enf_data_file=args.enf_data,
            targets=args.targets_in,
            filter=set,
            bin_size=args.predict_bin_size,
            aggregate=args.tar_aggr,
        )

        # Save queries for which targets could be extracted
        queries_file = f"{args.out_name}_queries_intersect_{set}.tsv"

        if args.gz:
            logger.info("gzipping queries intersect output ...")
            queries_file = queries_file + ".gz"

        out_queries.to_csv(queries_file, sep="\t", index=False)

        if args.store_text:
            targets_file = f"{args.out_name}_targets_{set}.tsv"
            if args.gz:
                logger.info("gzipping targets text output ...")
                targets_file = targets_file + ".gz"
            np.savetxt(targets_file, out_targets, delimiter="\t", index=False)

        if args.store_npy:
            targets_file = f"{args.out_name}_targets_{set}.npy"
            if args.gz:
                logger.info("gzipping targets npy output ...")
                targets_file = targets_file + ".gz"
            np.save(targets_file, out_targets)

        # save as h5py
        if args.store_h5:
            out_file_h5 = f"{args.out_name}_targets.h5"
            dataset_string = f"targets_{set}"
            with h5py.File(out_file_h5, "a") as h:
                h.create_dataset(dataset_string, data=out_targets)

        # TODO implement a toggle to write everything to a single set
        #   collapsing test, valid, train for re-sampling later
