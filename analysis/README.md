# Analysis

Platform-agnostic (or cross-platform) analysis that consumes the objects produced by
[`../preprocessing`](../preprocessing) and writes derived result tables used by
[`../figures`](../figures).

## [`FICTURE/`](FICTURE) — reference-guided spatial factor analysis (Xenium / MERSCOPE)

Segmentation-free spatial factor decomposition with [FICTURE](https://seqscope.github.io/ficture/),
guided by a fixed reference model built from the matched snRNA-seq annotation (see
[`scRNA_cell_type_annottaion/03_scrna_annotate_cell_types.R`](scRNA_cell_type_annottaion/03_scrna_annotate_cell_types.R)),
rather than the default unsupervised LDA factorization. Snakemake workflows are in
`FICTURE/Xenium/` and `FICTURE/MERSCOPE/`, shared utility scripts in `FICTURE/ficture_scripts/`.

## [`ST_alignment/`](ST_alignment) — Xenium/MERSCOPE to Visium spatial alignment

Registers Xenium/MERSCOPE cell coordinates onto the matched Visium H&E image with STalign
(landmark-based LDDMM), then bins the aligned cells into Visium spot-sized hexagons so all three
platforms can be compared at matched spatial units. See [`ST_alignment/README.md`](ST_alignment/README.md).

## [`scRNA_cell_type_annottaion/`](scRNA_cell_type_annottaion) — sc/snRNA-seq reference

Preprocessing and manual cluster annotation of the matched sc/snRNA-seq reference samples into the
final `cell_type202605` labels used for reference-based annotation throughout this repo. See
[`scRNA_cell_type_annottaion/README.md`](scRNA_cell_type_annottaion/README.md).

## [`iST_cell_type_annotation/`](iST_cell_type_annotation) — reference-based cell type annotation

Transfers `cell_type202605` reference labels onto Xenium/MERSCOPE data by SingleR
(`02_annotate_iST_cell_types.R`). See [`iST_cell_type_annotation/README.md`](iST_cell_type_annotation/README.md).

## [`run_cell_segmentation/`](run_cell_segmentation) — segmentation method comparison (Baysor vs proseg_v3)

Runs Baysor and proseg (v3) segmentation on raw Xenium/MERSCOPE transcripts, converts each to a
Seurat object, then applies the same QC/clustering/SingleR-annotation step used in
`iST_cell_type_annotation/` so the two methods' outputs are directly comparable. See
[`run_cell_segmentation/README.md`](run_cell_segmentation/README.md).

## [`process_aligned_data/`](process_aligned_data) — cross-platform data integration

Combines Visium with Xenium/MERSCOPE (via `ST_alignment`'s hexbin alignment) into per-sample,
cross-platform Seurat object lists with alignment QC metrics, then subsets to the spatial region
where every available platform has cells — the shared data structure used by downstream
cross-platform comparisons. See [`process_aligned_data/README.md`](process_aligned_data/README.md).

## [`sensitivity_metrics/`](sensitivity_metrics) — bin-level detection sensitivity

Bin-level transcript/gene detection metrics per platform (Xenium/MERSCOPE hex-binned at 110um,
Visium's native spots), before and after cross-platform alignment, for all genes and the shared
gene panel — feeds [`../figures/Fig2.R`](../figures/Fig2.R). See
[`sensitivity_metrics/README.md`](sensitivity_metrics/README.md).

## [`specificity_metrics/`](specificity_metrics) — negative-control-based detection specificity

Negative-control-probe and negative-control-codeword signal (off-target binding / decoding error)
per platform, at the 110um Visium-matched hexbin level and per individual probe/codeword, plus
Global Moran's I spatial-randomness testing for both negative controls and non-control genes —
feeds [`../figures/Fig3.R`](../figures/Fig3.R). See
[`specificity_metrics/README.md`](specificity_metrics/README.md).

## [`cell_segmentation_metrics/`](cell_segmentation_metrics) — segmentation method comparison metrics

QC metrics, marker-gene mutually exclusive co-expression rate, cell type composition and
silhouette-width separability, comparing Default/Baysor/Proseg segmentation across samples and
platforms — feeds [`../figures/Fig4.R`](../figures/Fig4.R). See
[`cell_segmentation_metrics/README.md`](cell_segmentation_metrics/README.md).

## [`cell_type_composition/`](cell_type_composition) — cross-platform cell type proportions

Cell type proportions per sample and platform (sc/snRNA-seq reference vs Xenium vs MERSCOPE), after
cross-platform alignment — feeds [`../figures/Fig5.R`](../figures/Fig5.R) panel b. See
[`cell_type_composition/README.md`](cell_type_composition/README.md).

## [`cell_cell_communication/`](cell_cell_communication) — ligand-receptor interaction analysis

Bin-level and cell-type-level ligand-receptor interaction testing (BLISA/CCI) for one sample
(MH0026 / "TNBC-02"), with pathologist-annotated tissue regions overlaid — feeds
[`../figures/Fig5.R`](../figures/Fig5.R) panels a, d, e, f. See
[`cell_cell_communication/README.md`](cell_cell_communication/README.md).

## Planned / not yet added

Further platform-agnostic result-table scripts (e.g. spatial statistics) will be added here as
they're written.
