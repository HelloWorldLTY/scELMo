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
..Input::
    A single cell (RNA) AnnData object with .obs being the genes and
    .var being the individual cells [gene x cell].
    Expects a .var column matching the cell_type_col_name argument.
    Expects .obs to be indexed of gene ID or symbols.
    Expects a .obs column matching the gene_col_name argument if one was
    provided that is not 'index'. If index is provided will use the gene ID
    index instead of a gene name.
    Expects a layer matching the layer argument to be present if specified.

..Arguments::
  -h, --help            Show this help message and exit
  --in IN_FILE          Input file is an anndata object saved as h5ad file.
  --genes GENES
                        List gene ids or symbols to compute the pseudobulk
                        aggregate for. Must match the entries in gene_col_name
                        of the anndata object. Default = ''.
  --gene_col_name GENE_COL_NAME
                        Name of .obs column where gene names can be found
                        that should be used for the aggregation. If set to
                        'index' will use the .obs index instead.
                        Default='index"
  --cell_type_col_name CELL_TYPE_COL
                        Name of the .var column that indicates the cell types
                        that will be used for the pseudobulking.
                        Default='cell types
  --method METHOD
                        Method to use for pseudobulking, supports:
                        'mean' - take the mean of the reads per gene
                        per cell type
                        'sum' - take the sum of the reads per gene per cell
                        type
                        'count_exp' - count the cells that express the gene at
                        or above an expression threshold provided per gene and
                        cell type
                        'perc_exp' - calculate the fraction of cells that
                        express
                        the gene at or above an expression threshold provided per
                        gene and cell type.
                        Default = 'mean'
  --expr_threshold EXP_THRESHOLD
                        Threshold at or above which a gene should be
                        considered as expressed. Matching the observed counts
                        in the anndata object.
  --layer LATER
                        If provided will use the anndata layer instead of the .X
                        counts.

..Usage::
python ./pseudobulk_anndata.py \
    --in my_anndata.h5ad \
    --out my_pseudobulked_anndata.h5ad  \
    --gene_col_name 'index' \
    --cell_type_col_name 'cell types'\
    --method 'mean'

..Output:: Output is AnnData object stored as .h5ad file under the --out
    location, with .obs being the genes and .var being the individual
    cell types [gene x cell types]. Where observed  counts were aggregated
    according to the chosen method.
"""
import argparse
import logging

import scanpy as sc

from seq2cells.utils.anndata_utils import pseudo_bulk

parser = argparse.ArgumentParser(
    description="Pseudobulk an AnnData object by cell type."
)
parser.add_argument(
    "--in",
    dest="in_file",
    type=str,
    required=True,
    help="Input anndata file in .h5ad format.",
)
parser.add_argument(
    "--genes",
    dest="genes",
    nargs="+",
    default="",
    required=False,
    help="List gene ids or symbols to compute the pseudobulk aggregate for. "
    "Must match the entries in gene_col_name of the anndata object.",
)
parser.add_argument(
    "--out",
    dest="out_file",
    default="./query_file_seq_model.tsv",
    type=str,
    required=True,
    help="Path and name for storing the pseudobulked anndata .h5ad",
)
parser.add_argument(
    "--gene_col_name",
    dest="gene_col_name",
    default="index",
    type=str,
    required=False,
    help="Name of .obs column where gene names can be found that should be "
    "used for the aggregation. If set to 'index' will use the .obs "
    "index instead. Default='index",
)
parser.add_argument(
    "--cell_type_col_name",
    dest="cell_type_col_name",
    default="cell types",
    type=str,
    required=False,
    help="Name of the .var column that indicates the cell types "
    "that will be used for the pseudobulking. "
    "Default='cell types",
)
parser.add_argument(
    "--method",
    dest="method",
    default="mean",
    type=str,
    required=False,
    help="Method to use for pseudobulking, supports:"
    "'mean' - take the mean of the reads per gene per cell type"
    "'sum' - take the sum of the reads per gene per cell type"
    "'count_exp' - count the cells that express the gene at or above an "
    "expression threshold provided per gene and cell type"
    "'perc_exp' - calculate the fraction of cells that express the "
    "gene at or above an expression threshold provided per gene and cell type. "
    "Default = 'mean",
)
parser.add_argument(
    "--expr_threshold",
    dest="expr_threshold",
    default=0.5,
    type=float,
    required=False,
    help="Threshold at or above which a gene should be considered as expressed. "
    "Matching the observed counts in the anndata object. Default = 0.5",
)
parser.add_argument(
    "--layer",
    dest="layer",
    default=None,
    type=str,
    required=False,
    help="If provided will use the anndata layer instead of the .X counts.",
)
parser.add_argument(
    "--mem_friendly",
    dest="mem_friendly",
    action="store_true",
    help="Flag to run in memory friendly mode. Takes oj the order of 10 times longer.",
)
parser.set_defaults(mem_friendly=False)
parser.add_argument(
    "--debug", dest="debug", action="store_true", help="Flag switch on debugging mode."
)
parser.set_defaults(debug=False)


if __name__ == "__main__":
    # fetch arguments
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # set scanpy verbosity
    # verbosity: errors (0), warnings (1), info (2), hints (3)
    if args.debug:
        sc.settings.verbosity = 3
    else:
        sc.settings.verbosity = 1

    # assert valid aggregation method selected
    assert args.method in [
        "mean",
        "sum",
        "perc_exp",
        "count_exp",
    ], "Invalid aggregation method selected!"

    # read anndata
    adata = sc.read_h5ad(args.in_file)

    # check selected genes
    if args.genes == "":
        genes = []
        num_genes = "all"
    else:
        genes = args.genes
        num_genes = len(genes)

    # run pseudobulking
    logger.info(f"Pseudo bulking {num_genes} genes ...")

    if args.mem_friendly:
        pseudo_adata = pseudo_bulk(
            adata,
            genes=genes,
            cell_type_col=args.cell_type_col_name,
            gene_col=args.gene_col_name,
            mode=args.method,
            expr_threshold=args.expr_threshold,
            mem_efficient_mode=True,
            layer=args.layer,
        )
    else:
        pseudo_adata = pseudo_bulk(
            adata,
            genes=genes,
            cell_type_col=args.cell_type_col_name,
            gene_col=args.gene_col_name,
            mode=args.method,
            expr_threshold=args.expr_threshold,
            mem_efficient_mode=False,
            layer=args.layer,
        )

    logger.info("Writting results to " + args.out_file)
    pseudo_adata.write(args.out_file)
