# Preprocessing

Raw vendor output → QC-filtered, clustered, cell-type-annotated Seurat object, one platform at a
time. Every script here operates on **one sample of one platform** — there is no batch loop. 

---

## `00_scRNA/` — sc/snRNA-seq reference

Preprocessing and manual annotation of the matched sc/snRNA-seq samples: the reference used for
cell type annotation throughout the repository.

| Script | Description |
|---|---|
| `01_reprocess.R` | Cell Ranger matrix → QC-filtered, normalised, embedded and clustered object. Prints a top-10-marker heatmap for manual cluster inspection. Output: `seu.rds` |
| `02_inferCNV.R` | Infers copy number variation from expression to identify malignant cells. Takes `<seu.rds> <sample_name>`; needs `seurat_clusters` metadata. Writes to `inferCNV_res/` |
| `03_annotate_cell_types.R` | Manual, iterative cluster → cell type annotation producing the final `cell_type202605` labels. Interactive: run one sample's block at a time and inspect the plots. Output: `seu_<sample>.rds` |

**Parameters — `01_reprocess.R` QC thresholds**

| Sample | Type | min_nFeature | max_nCount | max_percent_mt |
|--------|------|--------------|------------|----------------|
| TNBC_01 | snRNA | 200 | 10,000 | 5  |
| TNBC_02 | snRNA | 200 | 17,000 | 10 |
| TNBC_03 | scRNA | 200 | 10,000 | 16 |
| TNBC_04 | scRNA | 100 | 16,000 | 30 |
| ER_01   | scRNA | 150 | 12,000 | 20 |
| ER_02   | scRNA | 200 | 10,000 | 10 |

**`02_inferCNV.R` outputs**, written to `inferCNV_res/`:

| File | Description |
|------|-------------|
| `ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png` | CNV heatmap per cell group |
| `Boxplot_infercnv_instability_score.png` | Per-group CNV instability score boxplot |
| `instability_score_from_infercnv.rds` | Per-cell instability scores |
| `heatmap_list_combined.rds` | Combined heatmap object |
| `run.final.infercnv_obj` | Final inferCNV object |

---

## `visium/` — Visium

| Script | Description |
|---|---|
| `01_preprocess_visium.R` | Space Ranger `outs/` → QC-filtered, normalised and clustered object, up to and including spot clustering. Output: `so_raw.rds` (unfiltered, with `percent.mt`) and `seu_ST_Filtered.rds` |

**Parameters — spot QC thresholds and clustering resolution**

| Sample | min_nCount_Spatial | min_nFeature_Spatial | max_percent_mt | resolution | Notes |
|--------|--------------------|----------------------|----------------|------------|-------|
| TNBC_01 | 1,000 | 1,600 | 3.5  | 0.3 | |
| TNBC_02 | 1,000 | 2,500 | 10   | 0.4 | also filtered `nCount_Spatial` < 100,000 |
| TNBC_03 | 1,500 | 1,000 | 11   | 0.4 | |
| TNBC_04 | 1,500 | 1,000 | 2.5  | 0.4 | |
| ER_01   | 1,500 | 1,000 | 7    | 0.4 | |
| ER_02   | 500   | 500   | 7    | 0.4 | |

**Fixed across all samples:** genes expressed in fewer than 10 spots, or with total detected
expression < 30, are removed before normalisation; 2,000 HVGs; `dims = 1:30` for PCA/UMAP.

---

## `iST/` — Xenium / MERSCOPE

The **default (vendor) segmentation** arm. The two alternative segmentation methods produce the
same `so.rds` artifact from the same raw transcripts — see
[`../00_segmentation_baysor`](../00_segmentation_baysor) and
[`../00_segmentation_proseg_v3`](../00_segmentation_proseg_v3).

| Script | Description |
|---|---|
| `01_preprocess_iST.R` | Vendor cell-by-gene matrix and cell metadata → QC-filtered, normalised and clustered object, up to and including cell clustering. Set `platform` to `"Xenium"` or `"MERSCOPE"`; it selects the raw file format read. Output: `so_raw.rds` (unfiltered, with a `spatial` DimReduc of cell centroids and a `density` column) and `so.rds` |
| `02_annotate_iST.R` | Transfers `cell_type202605` from the matched reference by SingleR (Wilcoxon DE-based), giving `SingleR_labels202605`. Reference cells labelled `mixture` are dropped first. Output: `so.rds`, updated in place |

**Parameters — dimensionality and clustering resolution**

| Platform | Sample | dimUsed | resolution |
|----------|--------|---------|------------|
| Xenium   | TNBC_01 | 30 | 0.3 |
| Xenium   | TNBC_02 | 30 | 0.3 |
| Xenium   | TNBC_03 | 30 | 0.3 |
| Xenium   | TNBC_04 | 30 | 0.3 |
| Xenium   | ER_01   | 30 | 0.3 |
| Xenium   | ER_02   | 30 | 0.3 |
| MERSCOPE | TNBC_01 | 30 | 0.3 |
| MERSCOPE | TNBC_02 | 30 | 0.3 |
| MERSCOPE | TNBC_03 | 30 | 0.3 |
| MERSCOPE | ER_01   | 50 | 0.1 |
| MERSCOPE | ER_02   | 50 | 0.1 |

**Fixed across all samples and both platforms:** cells kept if `nCount_RNA >= 10` and
`nFeature_RNA >= 5`.

**Raw files expected** in each sample's directory (set in [`../config.R`](../config.R)):

| Platform | Files |
|----------|-------|
| Xenium | `cell_feature_matrix.h5`, `cells.csv.gz` |
| MERSCOPE | `cell_by_gene.csv`, `cell_metadata.csv` |

**Reference used by `02_annotate_iST.R`:** each sample's own matched reference from `00_scRNA/`,
read from `<processed_data_dir>/scRNA/<sample>/seu_<sample>.rds` regardless of modality — snRNA for
TNBC_01 and TNBC_02, scRNA for the rest.

TNBC_04 has no MERSCOPE row: its MERSCOPE run is excluded for poor quality, so that sample is
processed for Xenium alone.
