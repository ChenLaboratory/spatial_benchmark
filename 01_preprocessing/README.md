# Preprocessing

Stage 1: raw vendor output → annotated Seurat objects. Every script here operates on **one sample
of one platform at a time** — there is no batch loop. Each is a generic template; the
sample-specific parameters (QC thresholds, clustering resolution) are tabulated below rather than
hardcoded, since they differ by sample.

Set `raw_data_dir` and `processed_data_dir` in the path block at the top of each script — see the
[root README](../README.md#running-it). Raw data comes from GEO GSE341352, laid out as
`<platform>/<sample>/`.

Sample IDs in the tables below are the **manuscript** IDs, which are also the GEO IDs and the IDs
in the `sample` column of every object. They are listed in the
[root README](../README.md#samples).

## Order of work

The sc/snRNA-seq experiment was performed **first**, and its manual cluster annotation produces the
`cell_type202605` labels that every other platform is annotated against. So it is documented first
here, and must be run first:

```
scRNA/  ──►  cell_type202605 reference labels
              │
              ├──►  iST/02_annotate_iST.R      (SingleR label transfer onto Xenium/MERSCOPE)
              ├──►  ../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R
              │                                (same transfer, Baysor/Proseg objects)
              └──►  ../02_analysis/01_FICTURE/ (reference model for guided factorization)

visium/  and  iST/  are independent of each other and can be run in either order.
```

---

# `scRNA/` — sc/snRNA-seq reference

Preprocessing and manual annotation of the matched sc/snRNA-seq reference samples. This is the
reference used for cell type annotation and deconvolution throughout the repository.

## `01_reprocess.R`

Applied to each sample.

**Input:** 10X data directory, passed to `Read10X()`.

**Output:** `seu.rds` — Seurat object with low quality cells filtered out.

**QC filtering thresholds per sample:**

| Sample | Type | min_nFeature | max_nCount | max_percent_mt |
|--------|------|--------------|------------|----------------|
| ER_01   | scRNA | 150 | 12,000 | 20 |
| ER_02   | scRNA | 200 | 10,000 | 10 |
| TNBC_01 | snRNA | 200 | 10,000 | 5  |
| TNBC_02 | snRNA | 200 | 17,000 | 10 |
| TNBC_03 | scRNA | 200 | 10,000 | 16 |
| TNBC_04 | scRNA | 100 | 16,000 | 30 |

## `02_inferCNV.R` — infer copy number variation

Infers copy number variation from expression data to identify cancer cells.

```bash
Rscript 02_inferCNV.R <seu.rds> <sample_name>
```

| Argument | Description |
|----------|-------------|
| `seu.rds` | Seurat object with `seurat_clusters` metadata |
| `sample_name` | Sample label used in plot titles |

Also requires two external files not covered by the GEO accession:

- `seu_normal_human_breast.rds` — normal human breast scRNA-seq reference, used as the reference
  group for CNV inference
- `hg38_gencode_v27.txt` — gene order file (chromosome positions, GENCODE v27 / hg38)

**Output:** written to `inferCNV_res/` in the working directory:

| File | Description |
|------|-------------|
| `ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png` | CNV heatmap per cell group |
| `Boxplot_infercnv_instability_score.png` | Per-group CNV instability score boxplot |
| `instability_score_from_infercnv.rds` | Per-cell instability scores |
| `heatmap_list_combined.rds` | Combined heatmap object |
| `run.final.infercnv_obj` | Final inferCNV object |

## `03_annotate_cell_types.R` — cell type annotation

Manual, iterative cluster annotation (per sample, at increasing clustering resolution as needed)
into the final `cell_type202605` labels. Sources [`../R/plot_dotplot.R`](../R/plot_dotplot.R)
for marker dot plots.

**Output:** `seu_<sample>.rds` — reference object carrying `cell_type202605`.

Reference cells labelled `"mixture"` (unresolved sub-clustering) are excluded wherever this
reference is used for label transfer.

---

# `visium/` — Visium

## `01_preprocess_visium.R`

Applied to each sample, up to and including spot clustering (no SVG detection, marker/deconvolution
analysis or inferCNV).

**Input:** Space Ranger `outs/` directory, passed to `Load10X_Spatial()`.

**Output:** `seu_ST_Filtered.rds` — Seurat object with low-quality spots and lowly expressed genes
filtered out, normalized, and clustered (`seurat_clusters` in metadata).

**Gene-level filtering (fixed across all samples):** genes expressed in fewer than 10 spots, or with
total detected expression < 30, were removed prior to normalization.

**QC filtering thresholds and clustering resolution per sample:**

| Sample | min_nCount_Spatial | min_nFeature_Spatial | max_percent_mt | resolution | Notes |
|--------|--------------------|----------------------|----------------|------------|-------|
| ER_01   | 1,500 | 1,000 | 7    | 0.4 | |
| ER_02   | 500   | 500   | 7    | 0.4 | |
| TNBC_01 | 1,000 | 1,600 | 3.5  | 0.3 | |
| TNBC_02 | 1,000 | 2,500 | 10   | 0.4 | also filtered `nCount_Spatial` < 100,000 |
| TNBC_03 | 1,500 | 1,000 | 11   | 0.4 | |
| TNBC_04 | 1,500 | 1,000 | 2.5  | 0.4 | |

PCA/UMAP dimensionality (`dims = 1:30`) and the number of HVGs (2,000) were the same for every
sample.

---

# `iST/` — Xenium / MERSCOPE

This is the **default (vendor) segmentation** arm. The two alternative segmentation methods produce
the same `so.rds` artifact from the same raw transcripts — see
[`../00_segmentation_baysor`](../00_segmentation_baysor) and
[`../00_segmentation_proseg_v3`](../00_segmentation_proseg_v3).

## `01_preprocess_iST.R`

Applied to each sample of both platforms, up to and including cell clustering.

**Input:**

| Argument | Description |
|----------|-------------|
| `platform` | `"Xenium"` or `"MERSCOPE"` — selects the raw data format read in |
| Xenium output directory | Must contain `cell_feature_matrix.h5` and `cells.csv.gz` |
| MERSCOPE output directory | Must contain `cell_by_gene.csv` and `cell_metadata.csv` |

**Output:**

- `so_raw.rds` — unfiltered Seurat object with a `spatial` DimReduc (cell centroids) and a
  `density` metadata column (`nCount_RNA` / cell area or volume).
- `so.rds` — low-quality cells filtered out, normalized, and clustered (`seurat_clusters`).

**QC filtering thresholds (fixed across all samples and both platforms):** cells kept if
`nCount_RNA >= 10` and `nFeature_RNA >= 5`.

**Dimensionality and clustering resolution per sample:**

| Platform | Sample | dimUsed | resolution |
|----------|--------|---------|------------|
| Xenium   | ER_01   | 30 | 0.3 |
| Xenium   | ER_02   | 30 | 0.3 |
| Xenium   | TNBC_01 | 30 | 0.3 |
| Xenium   | TNBC_02 | 30 | 0.3 |
| Xenium   | TNBC_03 | 30 | 0.3 |
| Xenium   | TNBC_04 | 30 | 0.3 |
| MERSCOPE | ER_01   | 50 | 0.1 |
| MERSCOPE | ER_02   | 50 | 0.1 |
| MERSCOPE | TNBC_01 | 30 | 0.3 |
| MERSCOPE | TNBC_02 | 30 | 0.3 |
| MERSCOPE | TNBC_03 | 30 | 0.3 |
| MERSCOPE | TNBC_04 | —  | — (preprocessing not completed for this sample) |

## `02_annotate_iST.R` — reference-based cell type annotation

Transfers the `cell_type202605` reference labels from `scRNA/` onto each Xenium/MERSCOPE sample by
SingleR (Wilcoxon DE-based label transfer), producing `SingleR_labels202605`. Applied per sample and
platform.

**Input:**

| Argument | Description |
|----------|-------------|
| `sample_id` | Internal sample ID, e.g. `TNBC_01` |
| `platform` | `"Xenium"` or `"MERSCOPE"` |
| `ref_type` | `"scRNA"` or `"snRNA"` — which reference to use (sample-specific, see below) |
| `so.rds` | iST Seurat object from `01_preprocess_iST.R` |
| `seu_<sample_id>.rds` | Matched reference from `scRNA/03_annotate_cell_types.R` |

**Output:** `so.rds` (overwritten in place) — same object with an added `SingleR_labels202605`
metadata column.

**Reference type per sample:**

| Sample | Platforms | ref_type |
|--------|-----------|----------|
| ER_01   | Xenium, MERSCOPE | scRNA |
| ER_02   | Xenium, MERSCOPE | scRNA |
| TNBC_01 | Xenium, MERSCOPE | snRNA |
| TNBC_02 | Xenium, MERSCOPE | snRNA |
| TNBC_03 | Xenium, MERSCOPE | scRNA |
| TNBC_04 | Xenium only      | scRNA |

TNBC_04 has no MERSCOPE row — no completed MERSCOPE preprocessing for this sample, consistent with
the clustering table above.
