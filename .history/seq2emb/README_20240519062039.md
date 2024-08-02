# Preprocessing sequence embeddings

### 0) Intro - Workflow

Given some regions of interest (ROI), e.g. transcription start sites (TSS) the
aim of the sequence pre-processing
is to obtain:

1) A **query file** that specifies for each ROI: the DNA sequence
   window surrounding it and the location of the region of interest within this
   window.
2) Pre-computed DNA sequence **embeddings** for each ROI computed with the
   Enformer trunk
3) Gene (region) IDs that specifiy their intersection with Enformer training,
   test and validation sequences for splitting the dataset.

### 1) Query file

Enformer embeds and predicts over 896 bins of 128 bp covering the central
114,688 bp of the sequence queries of length 196,608 bp.
To extract embeddings of genomic ROIs, we construct sequence
queries of length 196,608 bp and identify the corresponding Enformer output
window
within which the ROI lies so the correct embedding can be extracted.

Using `create_seq_window_queries.py`

This script will take regions of interest, stitch them into patches if desired
and
construct sequence windows adhering to chromosome boundaries and create queries
for the sequence model.
Genomic position and the index (bin_id) of the prediction bin with which the
rois are intersecting are listed: 0-based!
The subsequent script calculating DNA sequence embeddings can then use the
bin_ids to extract the embeddings of interest.

Stitching: if enabled will group the rois based on a supplied grouping
variable (e.g. a gene name or id)
ROIs with the same grouping id will be grouped into patches. Patches are
split if they stretch over more than the supplied threshold (50kb default) and
sequence windows are constructed over the center of patches.
The position and bin id of the rois are listed in a comma separated string.

Notes:

* The stitching functionality is implemented but we do not use it for single
  cell expression predictions so far. To replicate the manuscript work run
  without stitching.

* ROIs only accept a single position, if larger regions of interests should be
  supplied then please center the coordinate first.

* By default, this script will center the sequence windows on ROIs or at the
  center of a stitched patch. Thus allowing predictions with a maximal
  sequence context reach for every roi.

* If the number of prediction bins is even, such as with the default
  Enformer setting, then the center of the sequence window is covered
  by the border of two bins. In that case the sequence window is shifted by
  minus half a bin size to center the ROI within a single bin.

#### Inputs:

1) A plain text file of ROI, where every line specifies a
   ROI supplied via the `--in` argument. Common formats are bed
   files or vcf file without header. Important, the genomic coodinates may be
   provided in bed-like (0-based, half open format) or as single column
   (1-based) vcf-like format.
   The coordinate handeling is controlled by the
   `position_col` and `position_base` arguments (see `--help`)

```angular2html
chr1    65418    65419    ENST00000641515.2    .    +    ENSG00000186092.7    OR4F5
chr1    451677    451678    ENST00000426406.4    .    -    ENSG00000284733.2    OR4F29
chr1    686653    686654    ENST00000332831.5    .    -    ENSG00000284662.2    OR4F16
chr1    923922    923923    ENST00000616016.5    .    +    ENSG00000187634.13    SAMD11
```

2) Reference genome in fasta file with .fai index present in same directory.

```angular2html
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

#### Usage:

```bash
python create_seq_window_queries.py \
    --in ./preprocessing_example_files/gencode.v41.basic.annotation.protein.coding.ensembl_canonical.tss.hg38.h10.bed \
    --ref_genome ./hg38.fa \
    --out ./query_tss_example.tsv \
    --chromosome_col 1\
    --position_col 3\
    --position_base 1 \
    --strand_col 6 \
    --group_id_col 7 \
    --additional_id_col 8 \
    --no-stitch
```

#### Output

Output is a tab-separated query file that lists the chrom start end strand of
the sequence window the ids of the stitched patch and the grouping and  
additional_id, the center of the sequence window the number of regions of
interest within the distance between multiple rois in the sequence and
the strands, position and bin id of the rois, comma separated if multiple
ones are available.

```angular2html
chr    seq_start    seq_end    seq_strand    patch_id    group_id    add_id    center    num_roi    stretch    strands_roi    positions_roi    bins_roi
chr1    1    196608    +    ENSG00000186092.7_0    ENSG00000186092.7    OR4F5    98369    1    0    ['+']    65419    191
chr1    353310    549917    -    ENSG00000284733.2_0    ENSG00000284733.2    OR4F29    451678    1    0    ['-']    451678    448
chr1    588286    784893    -    ENSG00000284662.2_0    ENSG00000284662.2    OR4F16    686654    1    0    ['-']    686654    448
chr1    825555    1022162    +    ENSG00000187634.13_0    ENSG00000187634.13    SAMD11    923923    1    0    ['+']    923923    448
chr1    860888    1057495    -    ENSG00000188976.11_0    ENSG00000188976.11    NOC2L    959256    1    0    ['-']    959256    448
chr1    862216    1058823    +    ENSG00000187961.15_0    ENSG00000187961.15    KLHL17    960584    1    0    ['+']    960584    448
chr1    868114    1064721    +    ENSG00000187583.11_0    ENSG00000187583.11    PLEKHN1    966482    1    0    ['+']    966482    448
chr1    883725    1080332    -    ENSG00000187642.10_0    ENSG00000187642.10    PERM1    982093    1    0    ['-']    982093    448
chr1    901729    1098336    -    ENSG00000188290.11_0    ENSG00000188290.11    HES4    1000097    1    0    ['-']    1000097    448
```

## 2) Sequence embeddings

The next step is to pre-compute the sequence embeddings over the ROIs now
specified in the query file.

Using `calc_embeddings_and_targets.py`

This script will take a query file as produced by
`create_seq_window_queries.py` and compute embeddings and optionally predicted
Enformer targets over the ROI.

Main idea here is that ROI are always centered on the
sequence model query window as much as possible to allow a balanced, maximal
sequence context for each prediction.

Ideally only a single region of interest or regions very close together are
supplied per query. Larger sets should be split in the prior pre-processing
step. E.g. split multiple clusters of TSS more than ~ 50 kb apart into
separate entities for summary later.

Embeddings and targets from multiple ROIs or with adjacent bins specified are
aggregated according to the specified methods. Default: Embeddings - mean,
Targets - sum.

#### Notes

* If the ROI / patch is located on the minus strand the reverse
  complement of the plus strand will be used as sequence input.
* If the reverse_complement is forced via `--rc_force` the reverse_complement
  is applied to plus strand patches and minus strand patches are processed
  from the plus strand. The position of ROIs are always
  mirrored where necessary to ensure the correct targets/embeddings are
  extracted.
* If the reverse complement augmentation is toggled on via `--rc_aug` then
  the reverse complement is applied randomly in 50 % of instances.
* `--rc_force` overwrites `--rc_aug`
* Shift augmentations are chosen randomly from the selected range of bp shifts
  selected a single bp shift if wanting to precisely control for that.
* Note: preprocessing with multiple ROIs per query is supported but all
  single cell work carried out by us was using a single ROI (TSS of
  canonical transcript).

#### Input

1) Query file as produced by `create_seq_window_queries.py` which is
   a raw text file
   including a header column that specifies the sequence windows to be
   processed by the seq model and the positions of the regions of interest
   within that sequence to be extracted (roi). Positions and bins of
   multiple ROI per query are comma separated in one string. Example format:

```angular2html
chr    seq_start    seq_end    seq_strand    patch_id    group_id    add_id    center    num_roi    stretch    strands_roi    positions_roi    bins_roi
chr1    1    196608    +    ENSG00000186092.7_0    ENSG00000186092.7    OR4F5    98369    1    0    ['+']    65419    191
chr1    353310    549917    -    ENSG00000284733.2_0    ENSG00000284733.2    OR4F29    451678    1    0    ['-']    451678    448
chr1    588286    784893    -    ENSG00000284662.2_0    ENSG00000284662.2    OR4F16    686654    1    0    ['-']    686654    448
chr1    825555    1022162    +    ENSG00000187634.13_0    ENSG00000187634.13    SAMD11    923923    1    0    ['+']    923923    448
chr1    860888    1057495    -    ENSG00000188976.11_0    ENSG00000188976.11    NOC2L    959256    1    0    ['-']    959256    448
chr1    862216    1058823    +    ENSG00000187961.15_0    ENSG00000187961.15    KLHL17    960584    1    0    ['+']    960584    448
chr1    868114    1064721    +    ENSG00000187583.11_0    ENSG00000187583.11    PLEKHN1    966482    1    0    ['+']    966482    448
chr1    883725    1080332    -    ENSG00000187642.10_0    ENSG00000187642.10    PERM1    982093    1    0    ['-']    982093    448
chr1    901729    1098336    -    ENSG00000188290.11_0    ENSG00000188290.11    HES4    1000097    1    0    ['-']    1000097    448
```

2) Reference genome in fasta format. Needs to be indexed (same name file
   with .fa.fai ending present)

```angular2html
>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

#### Usage

```bash
python calc_embeddings_and_targets.py \
--in_query ./preprocessing_example_files/query_tss_example.tsv \
--ref_genome hg38.fa \
--out_name enformer_out \
--position_base 1 \
--add_bins 0 \
--store_text \
--store_h5 \
--targets '4675:5312'  # for all Enformer cage-seq targets
```

#### Output

Output are one or two tab separated
text files storing the embeddings and optionally targets and/or an hdf5 file
storing the
embedding and target as pandas data frames under the 'emb' and 'tar' handle
respectively.
The header columns in the embedding file indicate the embedding dimensions.
The header columns in the target text file / data frame
correspond to the selected target ids (0-based) of Enformer targets
(see the
[published Basenji2 targets](https://github.com/calico/basenji/tree/master/manuscripts/cross2020)
).
Targets are subset to the selected targets, the indices of the selected are
stored in the header of the target output file (0-based)

Example raw text outputs:

```bash
head -n 3 enformer_out*tsv | cut -f 1,2,3
==> enformer_out_emb.tsv <==
0	1	2
-0.11201313883066177	-0.0001226698950631544	-0.10420460253953934
-0.1380479633808136	-8.836987944960129e-06	-0.14271216094493866

==> enformer_out_tar.tsv <==
4675	4676	4677
0.021540187299251556	0.012503976002335548	0.012968547642230988
0.01947534829378128	0.007085299119353294	0.007071667350828648
```

## 3) Intersect regions of interest with Enformer train / test / valid regions

For splitting genes into training, test and validation set we intersect the
position of their TSS with the regions over which Enformer is trained to
predict chromatin features and CAGE-seq coverage. See
[Kelley 2020](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008050#sec010)
For a description of the train, test, valid split region construction. The genes
whose TSS intersect with test and validation regions are extracted as test and
validation set for the single cell work. Where a TSS intersect with multiple
Enformer regions we select the one where the TSS is most central.

### Notes
 By default the Enformer input sequences are of length 196,608 bp.
 These regions were taken from the Basenji2 work with regions of length
 131,072 bp and extended by 32,768 bp to each side.
 The 131,072 bp sequences were shared by the authors.
 By default we trim the shared sequences to the central
 114,688 bp, because Enformer is only trained to predict over
 those  896 * 128 bp bins of each sequence window.
 The pruning can be disabled via the `--no_prune` flag. This will intersect
 the TSS with the 131,072 bp sequences.
 Alternatively, using `--extend` flag the sequence windows can be extended to
 the full 196,608 bp.

#### Input

1) Query file as produced by `create_seq_window_queries.py` which is
   a raw text file 
   including a header column that specifies the sequence windows to be
   processed by the seq model and the positions of the regions of interest
   within that sequence to be extracted (roi). Positions and bins of
   multiple ROI per query are comma separated in one string. 
   The 'patch_id' column is used for unique RSS/ROI identification 
    Example format:

```angular2html
chr    seq_start    seq_end    seq_strand    patch_id    group_id    add_id    center    num_roi    stretch    strands_roi    positions_roi    bins_roi
chr1    1    196608    +    ENSG00000186092.7_0    ENSG00000186092.7    OR4F5    98369    1    0    ['+']    65419    191
chr1    353310    549917    -    ENSG00000284733.2_0    ENSG00000284733.2    OR4F29    451678    1    0    ['-']    451678    448
chr1    588286    784893    -    ENSG00000284662.2_0    ENSG00000284662.2    OR4F16    686654    1    0    ['-']    686654    448
chr1    825555    1022162    +    ENSG00000187634.13_0    ENSG00000187634.13    SAMD11    923923    1    0    ['+']    923923    448
chr1    860888    1057495    -    ENSG00000188976.11_0    ENSG00000188976.11    NOC2L    959256    1    0    ['-']    959256    448
chr1    862216    1058823    +    ENSG00000187961.15_0    ENSG00000187961.15    KLHL17    960584    1    0    ['+']    960584    448
chr1    868114    1064721    +    ENSG00000187583.11_0    ENSG00000187583.11    PLEKHN1    966482    1    0    ['+']    966482    448
chr1    883725    1080332    -    ENSG00000187642.10_0    ENSG00000187642.10    PERM1    982093    1    0    ['-']    982093    448
chr1    901729    1098336    -    ENSG00000188290.11_0    ENSG00000188290.11    HES4    1000097    1    0    ['-']    1000097    448
```

2) Enformer sequences with train, test, validation assignment. The regions
   were [shared](https://console.cloud.google.com/storage/browser/basenji_barnyard/data)
   by the Basenji2/Enformer authors. And are also stored with thre files 
   required for pre-processing here ... #TODO

```angular2html
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
```

Using `intersect_queries_with_enformer_regions.py`

Run as

```bash
python intersect_queries_with_enformer_regions.py \
--query query_gencode_v41_protein_coding_canonical_tss_hg38_nostitch.tsv \
--enf_seqs sequences.bed \
--strip
```

#### Output

Three raw text files with the gene IDs belonging to train, test and 
validation set respectively. Those are used for 
tagging the genes in `add_embeddings_to_anndata.py`.

```bash
head -n 3 query_enf_intersect_*.txt
==> query_enf_intersect_test.txt <==
ENSG00000003096
ENSG00000004776
ENSG00000004777

==> query_enf_intersect_train.txt <==
ENSG00000000457
ENSG00000000460
ENSG00000000938

==> query_enf_intersect_valid.txt <==
ENSG00000000003
ENSG00000000005
ENSG00000000419
```
