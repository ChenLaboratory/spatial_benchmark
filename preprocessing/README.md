# Preprocessing

Run-once, platform-specific loading, QC filtering, normalization and clustering. Each script is a
generic template applied to every sample of its platform; sample-specific parameters (QC thresholds,
clustering resolution) are documented per section below rather than hardcoded, since they differ by
sample. Shared plotting/marker-dotplot helpers are sourced from [`../R`](../R).

## `00_visium.R` — Visium

Applied to each sample, up to and including spot clustering (no SVG detection, marker/deconvolution
analysis or inferCNV).

**Input:** Space Ranger `outs/` directory, passed to `Load10X_Spatial()`.

**Output:** `seu_ST_Filtered.rds` — Seurat object with low-quality spots and lowly expressed genes
filtered out, normalized, and clustered (`seurat_clusters` in metadata).

**Gene-level filtering (fixed across all samples):** genes expressed in fewer than 10 spots, or with
total detected expression < 30, were removed prior to normalization.

**QC filtering thresholds per sample:** spot-level QC thresholds and the clustering resolution were
sample-specific:

| Sample | min_nCount_Spatial | min_nFeature_Spatial | max_percent_mt | resolution | Notes |
|--------|---------------------|-----------------------|-----------------|------------|-------|
| ER_0114_T3 | 500 | 500 | 7 | 0.4 | |
| ER_0360 | 1,500 | 1,000 | 7 | 0.4 | |
| MH0007 | 1,000 | 1,600 | 3.5 | 0.3 | |
| MH0026 | 1,000 | 2,500 | 10 | 0.4 | also filtered nCount_Spatial < 100,000 |
| TN_B1_0177 | 1,500 | 1,000 | 11 | 0.4 | |
| TN_B1_0554 | 1,500 | 1,000 | 2.5 | 0.4 | |

PCA/UMAP dimensionality (`dims = 1:30`) and the number of HVGs (2,000) were the same for every sample.

## `01_iST.R` — Xenium / MERSCOPE

Applied to each sample of both platforms, up to and including cell clustering (no cell type
annotation — annotation uses the matched sc/snRNA-seq reference; see
[`../analysis/iST_cell_type_annotation`](../analysis/iST_cell_type_annotation)).

**Input:**

| Argument | Description |
|----------|-------------|
| `platform` | `"Xenium"` or `"MERSCOPE"` — selects the raw data format read in |
| Xenium output directory | Must contain `cell_feature_matrix.h5` and `cells.csv.gz` (10x Xenium output) |
| MERSCOPE output directory | Must contain `cell_by_gene.csv` and `cell_metadata.csv` (Vizgen output) |

**Output:**

`so_raw.rds` — unfiltered Seurat object with a `spatial` DimReduc (cell centroids) and a `density`
metadata column (`nCount_RNA` / cell area or volume).

`so.rds` — Seurat object with low-quality cells filtered out, normalized, and clustered
(`seurat_clusters` in metadata).

**QC filtering thresholds (fixed across all samples and both platforms):** cells were kept if
`nCount_RNA >= 10` and `nFeature_RNA >= 5`.

**Dimensionality and clustering resolution per sample:** `dimUsed` (PCA/UMAP/neighbor dimensions) and
`resolution` (`FindClusters`) were platform/sample-specific:

| Platform | Sample | dimUsed | resolution |
|----------|--------|---------|------------|
| Xenium | ER_0114 | 30 | 0.3 |
| Xenium | ER_0360 | 30 | 0.3 |
| Xenium | MH0007 | 30 | 0.3 |
| Xenium | MH0026 | 30 | 0.3 |
| Xenium | TN_0177 | 30 | 0.3 |
| Xenium | TN_0554 | 30 | 0.3 |
| MERSCOPE | ER_0114 | 50 | 0.1 |
| MERSCOPE | ER_0360 | 50 | 0.1 |
| MERSCOPE | MH0007 | 30 | 0.3 |
| MERSCOPE | MH0026 | 30 | 0.3 |
| MERSCOPE | TN_0177 | 30 | 0.3 |
| MERSCOPE | TN_0554 | — | — (preprocessing not completed for this sample) |

## sc/snRNA-seq reference preprocessing and annotation

Moved to [`../analysis/scRNA_cell_type_annottaion`](../analysis/scRNA_cell_type_annottaion)
(`01_scrna_reprocessing.R`, `02_scrna_inferCNV.R`, `03_scrna_annotate_cell_types.R`), since the
sc/snRNA-seq reference exists specifically to support iST (and Visium) cell type annotation/deconvolution
rather than as a platform preprocessed independently of that downstream use.
