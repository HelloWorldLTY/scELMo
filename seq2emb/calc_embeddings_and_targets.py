"""
Calculate embeddings and targets given sequence queries.
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
This script will take a query file as produced by
create_seq_window_queries.py and compute embeddings and predicted targets
over the regions of interest.

Main idea here is that regions of interest (ROI) are always centered on the
sequence model query window as much as possible to allow a balanced, maximal
sequence context for each prediction.

Ideally only a single region of interest or regions very close together are
supplied per query. Larger sets should be split in the prior pre-processing
step. E.g. split multiple clusters of TSS more than ~ 50 kb apart into
separate entities for summary later.

Embeddings and targets from multiple ROIs are aggregated according to
specified methods. Default: Embeddings - average, Targets sum up.
Better aggregation strategies are to be explored.

..Notes::
* Shift augmentations are chosen randomly from the selected range of bp shifts
selected a single bp shift if wanting to precisely control for that

* Reverse complement augmentation only applied if --rc is turned on. Randomly
reverse complements half the input sequences.

* Strandedness is taken into account '+' strand --> use forward strand
* '-' --> toggle --rc flag

..Input::

1) Query file as produced by create_seq_window_queries.py a raw text
file including a header column that specified the sequence windows to be
processed by the seq model and the positions of the regions of interest
within that sequence to be extracted (roi). Positions and bins of the rois
are comma separated in one string. (.txt or .tsv file)

example format:
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

 2) reference genome in fasta format. Needs to be indexed (same name file
 with .fa.fai ending present)

..Arguments::
  --in_query IN_QUERY   A region of interest seq model query file as produced
                        from pre-processing step (1) create_seq_windows.py
  --ref_genome REF_GENOME
                        Path to reference genome fasta file, needs to be
                        indexed and match the coordinate system in the query
                        file.
  --position_base POSITION_BASE
                        Position index 0 or 1 based default 1.
  --out_name OUT_NAME   Prefix for output files. Default="out_seq_model"
  --targets TARGETS_IN  String listing the indices and ranges of indices (all
                        0-based) indicating which targets to select from the
                        predicted output. Only relevant if target saving is
                        toggled on. Ranges are split and the last index per
                        range is extracted as well. E.g. "4675:5312" covers all
                        CAGE-seq data from Enformer.
                        both embeddings and targets are computed and saved.
  --store_text          Flag to indicate that results should be stored in a
                        tab separated .txt/.tsv file.
                        (tsv) file(s) will store an embedding and target file
                        separately.
  --store_h5            Flag to indicating to store results as h5py file(s)
                        (pandas). The embeddings are stored in the 'emb' slot
                        and the targets in 'tar'.
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
                        targets. Default = 128
  --emb_dim EMBEDDING_DIM
                        Dimension of the embedding produced. Default = 3072
  --emb_aggr EMB_AGGR   Select an aggregation method to summarize multiple
                        embeddings of interest per query if present. E.g.
                        multiple TSS. Also applies if extracting neighbouring
                        bins using --add_bins. The aggregation will be applied
                        along the instance dimension, the embedding dimension
                        remains unchanged. Valid values: 'mean', 'median',
                        'sum'. default='mean'
  --tar_aggr TAR_AGGR   Select an aggregation method to summarize multiple
                        target outputs of interest per query if present. E.g.
                        multiple TSS. Also applies if extracting neighbouring
                        bins using --add_bins. The aggregation will be applied
                        along the instance dimension the embedding dimension
                        remains unchanged. Valid values: 'mean', 'median',
                        'sum', 'min', max'. default='sum'
                        'sum', 'min', max'. default='sum'
  --shift_augs SHIFT_AUGS [SHIFT_AUGS ...]
                        List specifying the range of bp shifts to create
                        augmentations. Use single integer for a fixed
                        augmentation otherwise a random value in the range
                        will be selected. In theory every range of shifts is
                        accepted but will become less meaningful if shifting
                        by equal or more than 1x bin_size. E.g. for standard
                        Enformer work more than [-127, 127]. Leave blank to
                        not apply shifts.
                        Control aggregation with --emb_aggr. default=0
  --add_bins ADD_BINS   Number of neighbouring bins to extract surrounding the
                        bins of interest. Will extract +- <add_bins> bin
                        windows and aggregate them using the selected function.
                        default=0
  --rc_force            Flag to force convert sequence(s) to reverse
                        complement.
  --rc_aug              Flag to use reverse complement, randomly in 50 perc of
                        instances.
  --cached_model CACHED_MODEL
                        Optionally, specify a path to cached (Enformer) model.
                        If set to None will load Enformer model from the
                        public hugging face repo. Note this may clock up your
                        home dir .cache. default='None'
--debug                 Flag switch on debugging mode. Runs no seq model forward
                        pass but produces random numbers in the right format.

..Usage::
python $SCRIPT_DIR/calc_embeddings_and_targets.py \
  --in_query query_example_enformer_head2.tsv \
  --ref_genome $GENOME \
  --chromosome_sizes $CHROM_SIZES \
  --out_name test_emb_out \
  --position_base 1 \
  --add_bins 0 \
  --store_text \
  --store_h5 \
  --targets '4,9:11,0'

..Output:: Output depending on which was specified are one or two tab separated
text files storing the embeddings and targets and/or an hdf5 file storing the
embedding and target pandas data frames under the 'emb' and 'tar' handle
respectively. The header columns in the target text file / data frame
correspond to the selected target ids (0-based). E.g. Enformer targets.

"""


import argparse
import logging

import numpy as np
import pandas as pd
import torch
from enformer_pytorch import Enformer

import seq2cells.utils.data_seq as data_seq
from seq2cells.utils.data_utils import (
    aggregate_ndarray,
    get_roi_bins,
    mirror_roi_bins,
)
from seq2cells.utils.script_utils import split_index_range_mix

# for debugging
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="Calculate seq model embeddings and save targets if desired."
)
parser.add_argument(
    "--in_query",
    dest="in_query",
    type=str,
    required=True,
    help="A region of interest seq model query file as produced from "
    "pre-processing step (1) create_seq_windows.py",
)
parser.add_argument(
    "--ref_genome",
    dest="ref_genome",
    type=str,
    required=True,
    help="Path to reference genome fasta file, needs to be indexed and match "
    "the coordinate system in the query file.",
)
parser.add_argument(
    "--position_base",
    dest="position_base",
    default=1,
    type=int,
    required=False,
    help="Position index 0 or 1 based default 1.",
)
# output
parser.add_argument(
    "--out_name",
    dest="out_name",
    type=str,
    default="out_seq_model",
    help='Prefix for output files. Default="out_seq_model"',
)
parser.add_argument(
    "--targets",
    dest="targets_in",
    type=str,
    required=False,
    help="String listing the indeces and ranges of indices (all 0-based) "
    "indicating which targets to select from the predicted output. Only "
    "relevant if target saving if toggled on. Ranges are "
    "split and the last index per range is extracted as well. "
    'E.g "4675:5312" referring to all Enformer CAGE-seq target ids.',
)
# define what datatypes to store
parser.add_argument(
    "--store_text",
    dest="store_text",
    default=False,
    action="store_true",
    help="Flag to indicate that results should be stored in a "
    "tab separated .txt/.tsv file.",
)
parser.add_argument(
    "--gz",
    dest="gz",
    action="store_true",
    default=False,
    help="Flag indicating that text/tsv file should be gzipped.",
)
parser.add_argument(
    "--store_h5",
    dest="store_h5",
    default=False,
    action="store_true",
    help="Flag indicating to store results as h5py file(s) (pandas). The "
    "embeddings are stored in the 'emb' slot and the targets in 'tar'.",
)
# seq model specifics
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
    help="Bin size in bp that the sequence model predicts targets. Default = 128",
)
parser.add_argument(
    "--emb_dim",
    dest="embedding_dim",
    default=3072,
    type=int,
    required=False,
    help="Dimension of the embedding produced. Default = 3072",
)
# control aggregation
parser.add_argument(
    "--emb_aggr",
    dest="emb_aggr",
    default="mean",
    type=str,
    required=False,
    help="Select an aggregation method to summarize multiple embeddings of "
    "interest per query if present. E.g. multiple TSS. Also applies if "
    "extracting neighbouring bins using --add_bins. The aggregation will be "
    "applied along the instance dimension, the embedding dimension remains "
    "unchanged. Valid values: 'mean', 'median', 'sum'. default='mean'",
)
parser.add_argument(
    "--tar_aggr",
    dest="tar_aggr",
    default="sum",
    type=str,
    required=False,
    help="Select an aggregation method to summarize multiple target "
    "outputs of interest per query if present. E.g. multiple TSS. Also "
    "applies if extracting neighbouring bins using --add_bins. "
    "The aggregation will be applied along the instance dimension the "
    "embedding dimension remains unchanged. Valid values: 'mean', "
    "'median', 'sum', 'min', max'. default='sum'",
)
# augments
parser.add_argument(
    "--shift_augs",
    dest="shift_augs",
    nargs="+",
    default=[0],
    type=int,
    required=False,
    help="List specifying the range of bp shifts to create augmentations. "
    "Use single integer for a fixed augmentation otherwise a random "
    "value in the range will be selected. In theory every range of shifts is "
    "accepted but will become less meaningful if shifting by equal or more "
    "than 1x bin_size. E.g. for standard Enformer work more than [-127, "
    "127]. Leave blank to not apply shifts. Control aggregation with "
    "--emb_aggr. default=0",
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
    "--rc_force",
    dest="reverse_complement_force",
    default=False,
    action="store_true",
    help="Flag to force convert sequence(s) to reverse complement.",
)
parser.add_argument(
    "--rc_aug",
    dest="reverse_complement_aug",
    default=False,
    action="store_true",
    help="Flag to use reverse complement, randomly in 50 perc of instances.",
)
parser.add_argument(
    "--cached_model",
    dest="cached_model",
    default="None",
    type=str,
    required=False,
    help="Optionally, specify a path to cached (Enformer) model. If set to "
    "None will load Enformer model from the public hugging "
    "face repo. Note this may clock up your home dir .cache. "
    "default='None'",
)
# test flag using random numbers and not Enformer slow inference
parser.add_argument(
    "--debug",
    dest="debug",
    default=False,
    action="store_true",
    help="Flag switch on debugging mode. Runs no seq model forward pass but "
    "produces random numbers in the right format.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # make sure a valid embedding and target
    # aggregation method was provided
    # TODO implement embedding aggregation - use central one when extend
    #  around a single roi
    assert args.emb_aggr in ["sum", "mean", "median"], (
        "emb_aggr method not valid, " "check help!"
    )
    assert args.tar_aggr in [
        "sum",
        "mean",
        "median",
        "min",
        "max",
    ], "emb_aggr method not valid, check help!"

    # if reverse_complement_force and reverse_complement_aug turned on
    # overwrite augment and only use force, throw warning
    if args.reverse_complement_force and args.reverse_complement_aug:
        logger.warning(
            "--rc_force and --rc_aug were selected. rc_force will "
            "overwrite rc_aug and apply the reverse complement to "
            "every sequence!"
        )
        args.reverse_complement_aug = False

    # load files  ===========================
    logger.info("Reading queries ...")
    query = pd.read_csv(args.in_query, sep="\t", comment="#")

    # prep some fixtures
    # variable if supplied as list or ar range
    selected_targets = split_index_range_mix(args.targets_in)
    num_targets = len(selected_targets)
    num_queries = query.shape[0]

    # load model ============================
    # select device / gpu?
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # load a pretrained Enformer with DeepMind trained weights
    logger.info("Loading Enformer public model ...")
    if args.cached_model == "None":
        # use the public link to download the model to cache folder
        model = Enformer.from_pretrained("EleutherAI/enformer-official-rough")
    else:
        model = Enformer.from_pretrained(args.cached_model)

    # move to respective device
    logger.info("Moving model to device ...")
    model.to(device)

    # pre-process query =====================
    logger.info("Creating genome interval dataset ...")
    # create GenomeIntervalDataset
    genome_interval_dataset = data_seq.GenomeIntervalDataset(
        query=query,
        fasta_file=args.ref_genome,
        context_length=args.seq_window_size,
        filter_df_fn=data_seq.identity,
        query_to_fasta_map=dict(),
        return_seq_indices=False,
        shift_augs=args.shift_augs,
        rc_force=args.reverse_complement_force,
        rc_aug=args.reverse_complement_aug,
        return_augs=True,
        pos_base=args.position_base,
    )

    # initialize emtpy embedding and targets data frames
    emb_collect = np.zeros([num_queries, args.embedding_dim])
    tar_collect = np.zeros([num_queries, num_targets])

    # for each query
    logger.info("Running inference on the queries ...")
    with torch.no_grad():
        if not args.debug:
            model.eval()
        for idx in range(num_queries):
            # run the predictions ===================
            # call interval dataset
            seqs, strand, shift_augs, rc_aug = genome_interval_dataset[idx]
            if args.debug:
                # dummy numpy arrays for writing the post process and aggregate
                tar = np.random.rand(896, 5313)
                emb = np.random.rand(869, 3072)
            else:
                # Enformer inference run
                seqs = seqs.to(device)
                # get embeddings or targets or both
                tar, emb = model(seqs, return_embeddings=True)
                logger.info("Processed another ..")
                # move to cpu as numpy array
                emb = emb.cpu().numpy()
                tar = tar["human"].cpu().numpy()

            # post process ==========================
            # subset targets to the selected ones
            tar = tar[:, selected_targets]
            # get all bins of interest
            bins_of_interest = get_roi_bins(
                query["bins_roi"][idx], add_bins=args.add_bins
            )
            # mirror the bins of interest on the center of the sequence regions
            # if the strand is minus and rc_aug == False or if strand is plus
            # and rc_aug == True
            if (not rc_aug and strand == "-") or (rc_aug and strand == "+"):
                bins_of_interest = mirror_roi_bins(
                    bins_of_interest, args.predict_bin_size, args.predict_window_size
                )
            temp_emb = emb[bins_of_interest, :]
            temp_tar = tar[bins_of_interest, :]

            # aggregate emb and targets per query using the selected methods
            if temp_emb.shape[0] > 1:
                temp_emb = aggregate_ndarray(temp_emb, axis=0, method=args.emb_aggr)
            else:
                temp_emb = temp_emb[0, :]

            # targets
            if temp_tar.shape[0] > 1:
                temp_tar = aggregate_ndarray(temp_tar, axis=0, method=args.tar_aggr)
            else:
                temp_tar = temp_tar[0, :]

            # add to collection array
            emb_collect[idx, :] = temp_emb
            tar_collect[idx, :] = temp_tar

    # save ==================================
    logger.info("Converting and saving output ...")
    # convert to pandas dataframe
    emb_collect = pd.DataFrame(emb_collect)
    tar_collect = pd.DataFrame(tar_collect)
    # replace header with selected targets
    tar_collect.set_axis(selected_targets, axis=1, inplace=True)

    # save as text
    if args.store_text:
        if args.gz:
            out_file_emb_text = args.out_name + "_emb.tsv.gz"
            out_file_tar_text = args.out_name + "_tar.tsv.gz"
            logger.info("gzipping text output ...")
        else:
            out_file_emb_text = args.out_name + "_emb.tsv"
            out_file_tar_text = args.out_name + "_tar.tsv"
        emb_collect.to_csv(out_file_emb_text, sep="\t", index=False)
        tar_collect.to_csv(out_file_tar_text, sep="\t", index=False)

    # save as h5py
    if args.store_h5:
        out_file_h5 = args.out_name + ".h5"
        emb_collect.to_hdf(out_file_h5, key="emb", mode="w")
        tar_collect.to_hdf(out_file_h5, key="tar", mode="a")
