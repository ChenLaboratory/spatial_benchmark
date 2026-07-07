# scRNA/snRNA-seq Analysis

## Preprocessing

The script `01_preprocessing.Rmd` was applied to each sample.

### Input

| Argument | Description |
|----------|-------------|
| 10X data directory | Cell Ranger output directory passed to `Read10X()` |

### Output

`seu.rds` — Seurat object with low quality cells filtered out.

### QC filtering thresholds per sample

The following sample-specific thresholds were used:

| Sample | Type | nFeature_RNA | nCount_RNA | percent.mt |
|--------|------|-------------|------------|------------|
| ER_0114 | scRNA | > 200 | < 10,000 | < 10 |
| ER_0360 | scRNA | > 150 | < 12,000 | < 20 |
| TN_0177 | scRNA | > 200 | < 10,000 | < 16 |
| TN_0554 | scRNA | > 100 | < 16,000 | < 30 |
| MH0007  | snRNA | > 200 | < 10,000 | < 5  |
| MH0026  | snRNA | > 200 | < 17,000 | < 10 |


## inferCNV

Infers copy number variation from expression data to identify cancer cells.

### Input

```bash
Rscript run_inferCNV.R <seu.rds> <sample_name>
```

| Argument | Description |
|----------|-------------|
| `seu.rds` | Seurat object with `seurat_clusters` metadata |
| `sample_name` | Sample label used in plot titles |

The script also requires:
- A normal human breast scRNA-seq reference Seurat object (`seu_normal_human_breast.rds`), used as the reference group for CNV inference
- A gene order file (`hg38_gencode_v27.txt`), containing chromosome positions of genes based on GENCODE v27 / hg38

### Output

Results are written to `inferCNV_res/` in the working directory:

| File | Description |
|------|-------------|
| `ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png` | CNV heatmap per cell group |
| `Boxplot_infercnv_instability_score.png` | Per-group CNV instability score boxplot |
| `instability_score_from_infercnv.rds` | Per-cell instability scores |
| `heatmap_list_combined.rds` | Combined heatmap object |
| `run.final.infercnv_obj` | Final inferCNV object |

Results are visualized in `../scRNAseq_analysis_template.Rmd` under the `## inferCNV` section.
