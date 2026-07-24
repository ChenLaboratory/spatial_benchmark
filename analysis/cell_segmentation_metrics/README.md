# Cell Segmentation Metrics

Compares the three cell segmentation methods (Default, Baysor, Proseg v3 — see
[`../../run_cell_segmentation`](../../run_cell_segmentation)) across QC metrics, marker gene
co-expression, cell type composition and cluster separability. Feeds `figures/Fig4.R` (panels
b, c, f, g; panel e reads the annotated Seurat objects directly, see below).

```
01_qc_metrics.R              total/assigned transcripts, cell counts, cell area, density (Fig4b)
02_mecr.R                    mutually exclusive co-expression rate, after QC (Fig4c)
03_cell_type_composition.R   proportion of cells per annotated cell type (Fig4f)
04_silhouette_width.R        average silhouette width, annotated cell type grouping (Fig4g)
```

All four scripts read the same set of per-(sample, platform, segmentation-method) Seurat objects:
Default from `../../preprocessing`, Baysor/Proseg from `../../run_cell_segmentation` (each
script's own working directory — fill in the actual output paths for your run). Cell type labels
are the `SingleR_labels202605` column, produced for Default by
[`../iST_cell_type_annotation`](../iST_cell_type_annotation) and for Baysor/Proseg by
[`../../run_cell_segmentation/preliminary_analysis_seu/check_seu.R`](../../run_cell_segmentation/preliminary_analysis_seu/check_seu.R)
(QC filter → normalise/scale/PCA/UMAP/cluster → SingleR annotation against the matched sc/snRNA-seq
reference).

## 01 — QC metrics

Total transcripts (from the raw transcript file, shared by Default/Baysor/Proseg), transcripts
assigned to cells, total cell count, median cell area, and median transcripts per cell -- the four
metrics Fig4b plots (total cell count, median cell area, transcript assignment rate = assigned /
total transcripts, transcripts per cell). MERSCOPE's native cell area is a 3D volume; both
MERSCOPE Default and Proseg here use a 2D XY projected area instead (from
`cell_boundaries.parquet` / Proseg's own boundary polygons respectively). Baysor's cell area comes
from its own output metadata. All metrics are computed on the raw (pre-QC-filter) object -- this
script characterizes the raw segmentation output, it doesn't filter it (QC filtering happens
upstream, see below).

**Output**: `qc_df.rds` — one row per (sample, platform, segmentation).

## 02 — MECR (mutually exclusive co-expression rate)

For a panel of 23 canonical cell-type marker genes (restricted to the shared gene panel), MECR is
the Jaccard-style co-detection rate for each pair of markers assigned to *different* cell types —
elevated MECR indicates transcripts leaking between neighbouring cells (poor segmentation), though
a few pairs (e.g. myoepithelial or NK-cell markers) can reflect genuine co-expression biology
rather than segmentation error. Computed on QC-filtered (after-QC) datasets only, plus matched
sc/snRNA-seq as a reference baseline (already QC-processed). Summarised as median/mean MECR across
all gene pairs per (sample, platform, segmentation).

**Output**: `mecr_summary.rds`.

## 03 — cell type composition

Proportion of cells per annotated `SingleR_labels202605` type, per (sample, platform,
segmentation-method) dataset.

**Output**: `props_all.rds`.

## 04 — silhouette width (cell type)

Mean silhouette width (`bluster::approxSilhouette`, PCA space) using annotated cell type as the
grouping variable — higher values mean cells of the same annotated type are more transcriptionally
coherent and separable from other types, a direct readout of segmentation quality against expected
biology. Rare types (<10 cells) excluded; each retained type downsampled to at most 10,000 cells.

**Output**: `sil_ct_results.rds`.

## Not included

The source analyses also computed: cells-passing-QC count, median transcript density and median
genes/cell (01 — not among Fig4b's four plotted metrics); before-QC MECR and cluster-based
(unsupervised) silhouette width (using Seurat graph clustering instead of annotated cell type);
QC-metric boxplots/barplots and MECR boxplots/scatter-vs-%-assigned; per-cell-type cell area
distributions; composition heatmaps and by-method/by-platform (not by-sample) barplots. None of
these feed Fig4 panels b, c, e, f, g.
