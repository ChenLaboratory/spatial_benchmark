# sc/snRNA-seq Reference

Preprocessing and annotation of the sc/snRNA-seq reference samples used for iST (and Visium) cell
type annotation/deconvolution elsewhere in this repo.

## `01_scrna_reprocessing.R`

Applied to each sample.

**Input:** 10X data directory, passed to `Read10X()`.

**Output:** `seu.rds` — Seurat object with low quality cells filtered out.

**QC filtering thresholds per sample:**

| Sample | Type | min_nFeature | max_nCount | max_percent_mt |
|--------|------|-------------|------------|------------|
| ER_0114 | scRNA | 200 | 10,000 | 10 |
| ER_0360 | scRNA | 150 | 12,000 | 20 |
| TN_0177 | scRNA | 200 | 10,000 | 16 |
| TN_0554 | scRNA | 100 | 16,000 | 30 |
| MH0007  | snRNA | 200 | 10,000 | 5  |
| MH0026  | snRNA | 200 | 17,000 | 10 |

## `02_scrna_inferCNV.R` — infer copy number variation

Infers copy number variation from expression data to identify cancer cells.

```bash
Rscript 02_scrna_inferCNV.R <seu.rds> <sample_name>
```

| Argument | Description |
|----------|-------------|
| `seu.rds` | Seurat object with `seurat_clusters` metadata |
| `sample_name` | Sample label used in plot titles |

The script also requires:
- A normal human breast scRNA-seq reference Seurat object (`seu_normal_human_breast.rds`), used as the reference group for CNV inference
- A gene order file (`hg38_gencode_v27.txt`), containing chromosome positions of genes based on GENCODE v27 / hg38

**Output:** results are written to `inferCNV_res/` in the working directory:

| File | Description |
|------|-------------|
| `ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png` | CNV heatmap per cell group |
| `Boxplot_infercnv_instability_score.png` | Per-group CNV instability score boxplot |
| `instability_score_from_infercnv.rds` | Per-cell instability scores |
| `heatmap_list_combined.rds` | Combined heatmap object |
| `run.final.infercnv_obj` | Final inferCNV object |

## `03_scrna_annotate_cell_types.R` — cell type annotation

Manual, iterative cluster annotation (per sample, at increasing clustering resolution as needed) into
the final `cell_type202605` labels used as the reference cell type annotation throughout the rest of
the manuscript (Visium/iST deconvolution, FICTURE reference model, `../02_annotate_iST_cell_types.R`).
Sources [`../../../R/plot_dotplot.R`](../../../R/plot_dotplot.R) for marker dot plots.
