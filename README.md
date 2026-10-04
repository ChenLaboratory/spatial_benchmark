<!-- Optional: uncomment and add a logo to img/logo.png
<p align="center">
  <img src="img/logo.png" width="400" alt="Project logo">
</p>
-->

<h1 align="center">Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tumours </h1>

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Introduction](#introduction)
- [Data Availability](#data-availability)
- [Workflows](#workflows)
- [Repository layout](#repository-layout)
- [Citation](#citation)
- [Acknowledgements](#acknowledgements)

## Overview

This repository contains the code to reproduce the analysis in:

Qin L*, Yip RKH*, et al. (2026). **Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tumours**.


## Introduction

The field of spatial transcriptomics has dramatically expanded over the past few years, with the popularisation of a plethora of commercial platforms driving the growth. However, to date, no independent studies have comparatively analysed these techniques to benchmark their performance. Here, we present a well-controlled direct comparison of three spatial transcriptomics platforms. We evaluate these methods with multiple quality metrics and provide recommendations for data analysis and methods selection. This resource will help researchers to choose the optimal multiplexed in situ imaging technologies and is valuable for the development of new analysis tools.  


## Data Availability

The raw sequencing and ST data have been deposited in the Gene Expression Omnibus (GEO) under the
accession number [GSE344378].


| Manuscript ID | sc/snRNA | Visium | Xenium | MERSCOPE V1 | MERSCOPE V2 |
|---------------|----------|--------|--------|-------------|-------------|
| TNBC_01 | snRNA | ✓ | ✓ | ✓ | — |
| TNBC_02 | snRNA | ✓ | ✓ | ✓ | — |
| TNBC_03 | scRNA | ✓ | ✓ | ✓ | ✓ |
| TNBC_04 | scRNA | ✓ | ✓ | ✓ | — |
| TNBC_05 | —     | — | — | ✓ | ✓ |
| ER_01   | scRNA | ✓ | ✓ | ✓ | ✓ |
| ER_02   | scRNA | ✓ | ✓ | ✓ | ✓ |

## Workflows

The following analysis workflows are available in this repository, numbered in the order they run.
Each folder contains its own README with full usage instructions.

- [**00_segmentation_baysor/**](00_segmentation_baysor/) — Baysor resegmentation of the raw
  Xenium/MERSCOPE transcripts. 

- [**00_segmentation_proseg_v3/**](00_segmentation_proseg_v3/) — Proseg v3 resegmentation of the
  same transcripts, plus an optional 2D cell-area metric computed from Proseg's boundary polygons.

- [**00_ST_alignment/**](00_ST_alignment/) — Covers landmark-based STalign registration of iST cells onto the matched Visium H&E image,
  assignment of aligned cells to Visium spot-orientated 110 µm hexagons, and construction of the
  cross-platform Seurat object lists — including the overlapped region every comparison downstream
  depends on.

- [**01_preprocessing/**](01_preprocessing/) — Covers sc/snRNA-seq reference processing and cell type
  annotation, Visium spot QC
  and clustering, and Xenium/MERSCOPE processing under the default vendor segmentation, with SingleR label transfer from the matched reference.

- [**02_analysis/**](02_analysis/) — Cross-platform comparison. Covers detection sensitivity at matched spatial units, negative-control specificity and
  its spatial structure, segmentation-method benchmarking (QC metrics, marker co-expression,
  silhouette width), cell type composition, ligand-receptor interaction analysis, and guided
  FICTURE factor decomposition.

- [**03_figures/**](03_figures/) — One script per manuscript figure.



## Repository layout

```
config.R                             Data locations, sample lists and palettes — the only file you edit
R/
  blisa.R                            BLISA / CCI ligand-receptor method library
  plot_dotplot.R                     Marker dot plot helper

01_preprocessing/
  00_scRNA/
    01_reprocess.R                   QC, normalise, embed and cluster one sc/snRNA-seq sample
    02_inferCNV.R                    inferCNV instability score, used to call tumour clusters
    03_annotate_cell_types.R         Manual cluster annotation -> cell_type202605 labels
  visium/
    01_preprocess_visium.R           Space Ranger output -> QC-filtered, clustered object
  iST/
    01_preprocess_iST.R              Xenium/MERSCOPE vendor output -> Seurat object
    02_annotate_iST.R                SingleR label transfer from the matched reference

00_segmentation_baysor/
  01_run_baysor_xenium.sh            Segmentation (Xenium)
  01_run_baysor_merscope.sh          Segmentation (MERSCOPE)
  02_run_baysor_to_seu.sh            Segmentation output -> Seurat object
  generate_seu.R                     Seurat-conversion logic, shared by both platforms
  map_transcripts_baysor_*.py        Transcript mapping, called by 02
  xenium.toml, merscope.toml         Baysor config, one per platform

00_segmentation_proseg_v3/
  01_run_proseg_xenium.sh            Segmentation (Xenium)
  01_run_proseg_merscope.sh          Segmentation (MERSCOPE)
  02_run_proseg_to_seu.sh            Segmentation output -> Seurat object
  03_run_cell_area_2d.sh             Optional 2D cell area from Proseg boundary polygons
  generate_seu.R                     Seurat-conversion logic
  zarr_to_h5ad.py                    Format conversion, called by 02
  zarr_cell_area_2d.py               Cell-area computation, called by 03

00_ST_alignment/
  01_STalign_iST_to_Visium.py        Landmark-based registration onto the Visium H&E image
  02_hexbin.py                       Assign aligned cells to Visium spot-sized 110 um hexagons
  03_build_aligned_objects.R         Cross-platform object lists and the overlapped region

02_analysis/
  01_FICTURE/
    Xenium/, MERSCOPE/
      01_run_ficture_<platform>.smk  Snakemake workflow: prepare, fit, decode, plot
      01_submit.sh                   SLURM submission wrapper
      create_model_matrix.R          Guided reference model from annotated cell types
      parquet_to_csv.R               Raw transcripts -> FICTURE input CSV (Xenium)
    ficture_scripts/                 Helper scripts called by the workflow
  02_sensitivity_metrics/
    01_bin_metrics_before_alignment.R  Per-platform bin metrics, unaligned
    02_bin_metrics_after_alignment.R   Shared-panel bin metrics, overlapped region
    03_gene_count_scatter_data.R       Per-gene totals per sample and platform
  03_specificity_metrics/
    01_summed_control_counts.R       Summed negative-control signal per hexbin
    02_individual_control_moran.R    Moran's I per individual control feature
    03_individual_control_spatial.R  Per-hexbin control counts at 50 um, all samples
    04_individual_gene_moran.R       Moran's I and totals for non-control genes
  04_segmentation_metrics/
    00_preliminary_analysis_seu.R    QC, clustering and annotation of Baysor/Proseg objects
    01_qc_metrics.R                  Transcripts, cell counts, cell area
    02_mecr.R                        Mutually exclusive co-expression rate
    03_cell_type_composition.R       Composition per segmentation method
    04_silhouette_width.R            Silhouette width by annotated cell type
  05_cell_type_composition/
    01_proportion_after_alignment.R  Composition per sample x platform, after alignment
  06_cell_cell_communication/
    01_prepare_annotated_data.R      Cross-platform objects + pathology annotation
    02_run_blisa.R                   Bin-level ligand-receptor interaction testing
    03_run_cci.R                     Cell-type-level interaction scores
    04_region_lr_proportion.R        Significant-bin proportion per region x LR pair

03_figures/
  Fig2_sensitivity.R                 Figure 2: detection sensitivity across platforms
  Fig3_specificity.R                 Figure 3: negative-control specificity
  Fig4_segmentation.R                Figure 4: segmentation method comparison
  Fig5_celltype_cci.R                Figure 5: cell type composition and cell-cell communication
  Ext_Fig5_gene_scatter.R            Extended Data Figure 5: per-gene count scatter
  Ext_Fig6_specificity.R             Extended Data Figure 6: supplementary specificity panels
```


## Citation

If you use this dataset or these workflows, please cite:

Qin L*, Yip RKH*, et al. (2026). **Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tumours**.

> Full citation — complete author list, journal and DOI — to be added on publication.


## Acknowledgements

We thank the WEHI Advanced Genomics Facility, Center for Dynamic Imaging and Advanced Histotechnology facility for supporting the generation and analysis of single-cell and spatial transcriptomics data for this project.

