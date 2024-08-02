"""
Add precomputed gene/seq embeddings to AnnData object
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
..Input::
1) Query file as produced by create_seq_window_queries.py a tsv file
including a header column that specified the sequence windows to be
processed by the seq model and the positions of the regions of interest
within that sequence to be extracted (roi). Positions and bins of the rois
are comma separated in one string. (.txt or .tsv file). Example format:

chr	seq_start	seq_end	seq_strand	patch_id	group_id	add_id\
center num_roi	stretch	strands_roi	positions_roi	bins_roi
chr1	1	196608	+	ENSG00000186092.7_0	ENSG00000186092.7	OR4F5\
98369	1	0	['+']	65420	191
chr1	353311	549918	-	ENSG00000284733.2_0	ENSG00000284733.2	OR4F29\
451679	1	0	['-']	451679	448
chr1	588287	784894	-	ENSG00000284662.2_0	ENSG00000284662.2	OR4F16\
686655	1	0	['-']	686655	448

2) HDF5 the gene / TSS / ROI embeddings:
 Requires the embeddings stored in the 'emb' slot.
 Produced by 'calc_embeddings_and_targets.py'

3) Gene sets (test, validation)
Two raw txt files with IDs of the regions of interest / TSS / genes to split
into the test and validation set respectively. One ID per line. All other IDs
will be used as training set.

..Arguments::
  -h, --help            show this help message and exit
  --query QUERY_FILE    A region of interest seq model query file as produced
                        from pre-processing step (1) create_seq_windows.py
  --emb EMBEDDINGS_FILE
                        Path to .h5 file storing the precomputed sequence
                        embeddings in the 'emb' slot. The order must match
                        the order of the query df!
  --anndata ANNDATA_FILE
                        Path to .h5 stored AnnData object. Usually single cell
                        data. Must be in the format (cell x genes).
  --valid_ids VALID_IDS
                        Path to gene / ROI ids to use for the validation set.
                        Must match the gene/ROI IDs in the AnnData object.
  --test_ids TEST_IDS   Path to gene / ROI ids to use for the test set.
                        Must match the gene/ROI IDs in the AnnData object.
  --out_name OUT_NAME   Prefix for output hdf5 stored AnnData object.
                        Default="anndata_with_seq_embeddings"
  --reset_cell_index RESET_CELL_INDEX
                        New column name under which to store the current
                        cell index in the AnnData object when resetting
                        the cell index by Cell_1, Cell_2, ...
                        Default = None -> no resetting of the cell index.
  --strip               Flag to strip the gene/ROI/TSS IDs from versions.
                        Remove '\_\d+'
  --make_unique         Indicate to make the gene names in the AnnData object
                        unique. Usually necessary if gene symbols are stored
                        rather than gene IDs.
  --fetch_ids           Indicate to fetch gene IDs from the query file.
                        Will match the AnnData gene (.var) index to the 'add_id'
                        column assuming is hosts the gene symbols and fetch the
                        'group_id' column assuming it hosts unique gene ids
                        e.g. Ensembl IDs.
  --gene_id_col GENE_ID_COL
                        Only required if fetching gene IDs for the AnnData
                        object. Name of the gene column under which to store
                        the gene IDs in addition to the IDs becoming the
                        new index. Default='ens_id'
..Usage::

python add_embeddings_to_anndata.py \
    --query query_gencode_v41_protein_coding_canonical_tss_hg38_nostitch.tsv \
    --anndata single_cell_data.h5ad \
    --emb precomputed_embeddings.h5 \
    --valid_ids valid_ids.txt \
    --test_ids test_ids.txt \
    --strip

..Output::
HDF5 stored AnnData object (same as supplied) but transposed to (genes x cells)
and with the sequence embeddings added in the .obsm layers.
"""
import argparse
import logging

# Imports ==================================
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

# PARAMETERS ===============================
# define arguments ==============================
parser = argparse.ArgumentParser(
    description="Combine sequence embeddings and AnnData object."
)
# input
parser.add_argument(
    "--query",
    dest="query_file",
    type=str,
    required=True,
    help="A region of interest seq model query file as produced from "
    "pre-processing step (1) create_seq_windows.py",
)
parser.add_argument(
    "--emb",
    dest="embeddings_file",
    type=str,
    required=True,
    help="Path to .h5 file storing the precomputed sequence embeddings"
    "in the 'emb' slot. The order must match the order of the query df!",
)
parser.add_argument(
    "--anndata",
    dest="anndata_file",
    type=str,
    required=True,
    help="Path to .h5 stored AnnData object. Usually single cell data. "
    "Must be in the format (cell x genes).",
)
parser.add_argument(
    "--valid_ids",
    type=str,
    required=True,
    help="Path to gene / ROI ids to use for the validation set. Must match "
    "the gene/ROI IDs in the AnnData object.",
)
parser.add_argument(
    "--test_ids",
    type=str,
    required=True,
    help="Path to gene / ROI ids to use for the test set. Must match "
    "the gene/ROI IDs in the AnnData object.",
)
# output
parser.add_argument(
    "--out_name",
    dest="out_name",
    type=str,
    default="anndata_with_seq_embeddings",
    help="Prefix for output hdf5 stored AnnData object."
    'Default="anndata_with_seq_embeddings"',
)
# options
parser.add_argument(
    "--reset_cell_index",
    type=str,
    default=None,
    help="New column name under which to store the current cell index in the "
    "AnnData object when resetting the cell index by Cell_1, Cell_2, "
    "... Default = None -> no resetting of the cell index.",
)
# strip
parser.add_argument(
    "--strip",
    dest="strip",
    action="store_true",
    default=False,
    help='Flag to strip the gene/ROI/TSS IDs from versions. Remove "\\_\\d+"',
)
parser.add_argument(
    "--make_unique",
    action="store_true",
    default=False,
    help="Indicate to make the gene names in the AnnData object unique. "
    "Usually necessary if gene symbols are stored rather than gene IDs. ",
)
parser.add_argument(
    "--fetch_ids",
    action="store_true",
    default=False,
    help="Indicate to fetch gene IDs from the query file. Will match the "
    "AnnData gene (.var) index to the 'add_id' column assuming is "
    "hosts the gene symbols and fetch the 'group_id' column assuming it "
    "hosts unique gene ids e.g. Ensembl IDs.",
)
parser.add_argument(
    "--gene_id_col",
    type=str,
    default="ens_id",
    help="Only required if fetching gene IDs for the AnnData object. Name of "
    "the gene column under which to store the gene IDs in addition to "
    "the IDs becoming the new index. Default='ens_id'",
)


# FUNCTIONS ================================
def deduplicate_obs(
    scad: ad.AnnData, gene_symbol_col: str = "gene_symbol"
) -> ad.AnnData:
    """Remove all observations that have a duplicated index / gene name.

    Arguments
    ---------
    scad: ad.AnnData
        AnnData object in (gene x cell) orientation.
    gene_symbol_col: str
        Name of the .obs column under which to store the current index
        of the AnnData object. Default='gene_symbols'

    Returns
    -------
    scad: ad.AnnData
        AnnData object with duplicate entries removed. The provided .obs
        index is stored as an .obs column with the provided name.

    Notes
    -----
    Single cell data tends to be stored in AnnData objects either indexed by
    geme symbols (gene names) or by gene IDs. Gene names tend to be not
    unique requiring us to make them unique and optionally fetch the matching gene IDs.
    """
    # detach anndata object
    scad = scad.copy()

    idx = scad.obs.index.values

    # store in new .obs column
    scad.obs[gene_symbol_col] = idx

    scad.obs_names_make_unique()

    duplicates = scad.obs.duplicated(subset=[gene_symbol_col], keep="first")
    scad = scad[~duplicates, :]

    return scad


def fetch_gene_id(
    scad: ad.AnnData, query: pd.DataFrame, gene_id_col: str = "ens_id"
) -> ad.AnnData:
    """Swap gene symbols from an anndata object to ensembl ids

    Arguments
    ---------
    scad: ad.AnnData
        AnnData object in (gene x cell) orientation. Needs to  have
        deduplicated gene index (.obs.index).
    query: pd.DataFrame
        Sequence query dataframe containing for every gene (TSS, ROI) the
        matching sequence query for DNA models. For this function only
        the following columns are required:
            'group_id' - ID e.g. gene ID used to group , normally this is
            a unique gene ID e.g. Ensembl ID. This ID is used to connect
            regions of interest with the single cell data so the ID system
            must match.
            'add_id' - generated from preprocessing script stored gene names
            / symbols not used for as unique identifiers. This column must
            match the .obs.index in <scad>.
    gene_id_col: str
        Name of the .obs column under which to store the gene ids. The ids
        are stored as .obs column and than set as index.

    Returns
    -------
    scad: ad.AnnData
        AnnData object as provided but with the gene index (.obs.index)
        replaced by gene IDs as fetched from the query data frame.

    Notes
    -----
    Requires deduplication of gene symbols first!
    """
    # detach anndata object
    scad = scad.copy()

    # make a look-up dictionary from query df
    ens_dict = {
        name: group_id for name, group_id in query[["add_id", "group_id"]].values
    }

    # match gene symbols to gene ids
    scad.obs[gene_id_col] = [
        ens_dict[x] if x in ens_dict.keys() else "not_avail" for x in scad.obs.index
    ]

    # remove genes here where no gene id available from query dataframe
    scad = scad[scad.obs[gene_id_col] != "not_avail", :]

    # set gene id as new index
    scad.obs.index = scad.obs[gene_id_col]

    return scad


def prepare_anndata(
    anndata_file: str,
    gene_id_col: str,
    fetch_ids: bool,
    make_unique: bool,
) -> ad.AnnData:
    """Read and prepare the anndata object

    Arguments
    ---------
    anndata_file: str
        Path to AnnData object stored as .h5 file.
    gene_id_col: str
        Name of the column hast hosts the unique gene (region) ids to use.
    fetch_ids: bool
        If to grep unique gene IDs for the gene names stored in the AnnData
        object from the query file. The query file must host a unique gene ID.
        If the provided AnnData object only has non-unique gene names e.g.
        gene symbols, unique gene IDs can be fetched form the query file.
        However duplicate gene symbol entries are lost and only the
        respective first gene symbol entry is used.
    make_unique bool
        Toggle if to make the gene (region) entries of the AnnData object
        unique before progressing.


    Returns
    -------
     scad: ad.AnnData
        AnnData object in format (genes x cells), with deduplicates and gene
        ID fetched indices. Transposed from the input (cell x genes) format.
        And gene IDs stripped if desired.
    """
    scad = sc.read_h5ad(anndata_file)

    # check if adata gene entries are unique or if specified to make them unique
    if not make_unique:
        assert len(scad.var.index) == len(scad.var.index.unique()), (
            "Gene entries of the AnnData object are not unique and the "
            "--make_unique flag was not provided. "
            "Either use the --make_unique flag or make the gene index unique "
            "manually. Stopping ... "
        )

    # transpose anndata cell x gene -> gene x cell
    scad = scad.T.copy()

    if make_unique:
        logger.info(
            "Will make observations (gene names) unique. Adding -1, "
            "-2 etc. to subsequent appearances. This may lead to "
            "ambiguity when fetching gene ids from the query file as "
            "only the first appearance of a gene will be unaltered."
        )
        # if indicated make gene_symbols unique
        scad.obs_names_make_unique()

        # remove left over duplicates
        scad = deduplicate_obs(scad)

    if fetch_ids:
        # grep gene IDs from query data frame matched by gene symbol
        scad = fetch_gene_id(scad, query_df, gene_id_col)

    # strip gene ID
    if args.strip:
        scad.obs.index = scad.obs.index.str.replace("\\.\\w+", "", regex=True)

    return scad


def report_zero_frac(scad: ad.AnnData) -> None:
    """Report the fraction of zero entries in the obs matrix."""
    num_zeros = np.sum(adata.to_df() == 0).sum()
    total = adata.shape[0] * adata.shape[1]
    perc_zeros = np.round((num_zeros / total), decimals=5) * 100
    logger.info(
        f"In the supplied AnnData object {perc_zeros} % of " f"observations are zero."
    )


def add_seq_embeddings_to_anndata(
    scad: ad.AnnData,
    query: pd.DataFrame,
    emb_file: str,
    reset_var_index_name: str = None,
) -> ad.AnnData:
    """Combine a single cell anndata frame and sequence query dataframe.

    Arguments
    ---------
    scad: ad.AnnData
        Single cell AnnData object (genes x cells). Genes as (.obs) are
        indexed by gene ID matching the gene ID used in the sequence query
        data frame.
    query: pd.DataFrame
        Sequence query dataframe containing for every gene (TSS, ROI) the
        matching sequence query for DNA models. Requires the following
        columns to be present:
            'chr' - the chromosome
            'seq_start' - start of the sequence window
            'seq_end' - end of the sequence window
            'group_id' - ID e.g. gene ID used to group , normally this is
            a unique gene ID e.g. Ensembl ID. This ID is used to connect
            regions of interest with the single cell data so the ID system
            must match.
            'patch_id' - generated from preprocessing script to separate
            multiple stitched ROI patches e.g. from the same gene
            'center' - center of the sequence window
            'num_roi' - Number of regions of interest (ROI) (e.g. TSS) within the
            sequence window
            'stretch' - max distance between all ROI in the sequence window
            'positions_roi' - position of the region(s) of interest within a
            sequence window
            'bins_roi' - index/indices of the sequence window(s)/bin(s)
            (e.g. Enformer 128 bp windows)
            'strands_roi' - strands on which each ROI is located
    emb_file: str,
        Path to embedding file in hdf5 format with precomputed embeddings
        stored in the 'emb' slot. The order must match the order of the query df!
    reset_var_index_name: str = None
        The cell index here stored in (AnnData.var) is often the UMI barcode.
        Provide a new .var column name to store the current .var index and
        reset the .var index to Cell_1, Cell_2, ....
        Default = None -> no resetting of the cell index.

    Returns
    -------
    scad: ad.AnnData
        AnnData object subset to genes matching the provided sequence query
        dataframe with the sequence embeddings stored in the .obsm slot
        'seq_embedding' and with meta data added.

    Notes
    -----
    Match the sequence query IDs to the gene IDs in the AnnData object.
    Subset Anndata to the intersection.
    Optionally reset the Anndata index.
    Drops duplicates from the sequence query.
    Add sequence query meta data to AnnData object.
    """
    # detach anndata object
    scad = scad.copy()

    # read embeddings from file
    embeddings = pd.read_hdf(emb_file, "emb", "r")
    embeddings.index = query["group_id"]

    # Intersect query group_ids (e.g. gene IDs)
    # with geneIDs in single cell object
    intersect = set(scad.obs.index).intersection(query["group_id"].values)
    logger.info(
        f"Overlap is {len(intersect)} of {len(scad.obs.index)} enf TSS and "
        f"{len(query['group_id'].values)} single cell genes"
    )

    # subset AnnDataFrame
    scad = scad[list(intersect), :].copy()

    # optionally reset index to Cell_X
    if reset_var_index_name:
        scad.var[reset_var_index_name] = scad.var.index
        scad.var.index = [f"Cell_{i}" for i in range(scad.shape[1])]

    # drop duplicates from query df and embeddings
    subquery = query.loc[query["group_id"].isin(intersect), :]
    subemb = embeddings.loc[embeddings.index.isin(intersect), :]

    dups = sum(subquery.duplicated(subset=["group_id"]).values)
    logger.info(
        f"Dropping {dups} gene duplicate entries from query and embeddings. "
        f"Using first occurrence with ID matching single cell data!"
    )

    subquery = subquery.drop_duplicates(subset=["group_id"])
    subemb = subemb[~subemb.index.duplicated(keep="first")]

    # sort the subset query_df by the anndata index
    subquery.index = subquery["group_id"]
    subquery = subquery.loc[scad.obs.index]

    # add query meta data
    scad.obs["chr"] = subquery["chr"].values
    scad.obs["seq_start"] = subquery["seq_start"].values
    scad.obs["seq_end"] = subquery["seq_end"].values
    scad.obs["patch_id"] = subquery["patch_id"].values
    scad.obs["center"] = subquery["center"].values
    scad.obs["num_roi"] = subquery["num_roi"].values
    scad.obs["stretch"] = subquery["stretch"].values
    scad.obs["positions_roi"] = subquery["positions_roi"].values
    scad.obs["bins_roi"] = subquery["bins_roi"].values
    scad.obs["strands_roi"] = subquery["strands_roi"].values

    # remove index name from embeddings
    subemb.index.name = None

    # convert numeric column numbers to string
    subemb.columns = [str(c) for c in subemb.columns]

    # add embeddings as observations matrix
    scad.obsm["seq_embedding"] = subemb.loc[scad.obs.index]

    return scad


def label_train_test_val(
    test_ids: str, valid_ids: str, scad: ad.AnnData, split_label: str = "enf_set"
) -> ad.AnnData:
    """Label observations as train/test/valid set.

    Arguments
    ---------
    test_ids: str
        Path to test ids. A raw rext file listing a unique identifier per
        gene / region of interest. Must match the identifiers in the .obs
        index of the provided AnnData object.
    valid_ids: str
        Path to test ids. A raw rext file listing a unique identifier per
        gene / region of interest. Must match the identifiers in the .obs
        index of the provided AnnData object.
    scad: ad.AnnData
        AnnData object with unique gene IDs as .obs index.
        Those need to match to the supplied test and validation ids.
    split_label: str
        Name of the .obs column under which to store the train, test,
        valid assignment.

    Returns
    -------
    scad: ad.AnnData
        AnnData object with additional column names after 'split_label' with
        'train' 'test' and 'valid' assigned to each gene (.obs) entry.
    """
    # read test and validation ids to use
    test_ids = pd.read_csv(test_ids, header=None)
    valid_ids = pd.read_csv(valid_ids, header=None)

    if args.strip:
        test_ids[0] = test_ids[0].str.replace("\\.\\w+", "", regex=True)
        valid_ids[0] = valid_ids[0].str.replace("\\.\\w+", "", regex=True)

    # assign train / test /valid
    to_set = [
        "test"
        if i in test_ids.values
        else "valid"
        if i in valid_ids.values
        else "train"
        for i in adata.obs.index
    ]
    scad.obs[split_label] = to_set

    return scad


if __name__ == "__main__":

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Prepare sequence queries
    query_df = pd.read_csv(args.query_file, sep="\t")

    # strip gene IDs of version number
    if args.strip:
        query_df["group_id"] = query_df["group_id"].str.replace(
            "\\.\\w+", "", regex=True
        )

    adata = prepare_anndata(
        args.anndata_file, args.gene_id_col, args.fetch_ids, args.make_unique
    )

    report_zero_frac(adata)

    # Match the sequence embeddings to the single cell data
    adata = add_seq_embeddings_to_anndata(
        adata, query_df, args.embeddings_file, reset_var_index_name="barcode"
    )

    adata = label_train_test_val(args.test_ids, args.valid_ids, adata)

    adata.write(args.out_name + ".h5ad", compression="gzip")
