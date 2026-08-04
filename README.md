<!-- Optional: uncomment and add a logo to img/logo.png
<p align="center">
  <img src="img/logo.png" width="400" alt="Project logo">
</p>
-->

<h1 align="center">Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tissues </h1>

# Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tumours

## Contents
- [Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tissues](#cross-platform-spatial-transcriptomics-profiles-of-human-breast-tissues)
  - [Contents](#contents)
  - [Overview](#overview)
  - [Introduction](#introduction)
  - [Data Availability](#data-availability)
  - [Repository structure](#repository-structure)
  - [Citation](#citation)
  - [Acknowledgements](#acknowledgements)

## Overview

This repository contains the code to reproduce the analysis in:

Qin L*, Yip RKH*, et al. (2026). **Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tumours**.


## Introduction

<!-- TODO: 1–2 short paragraphs.
     What the study is, what data/platforms it covers, and what this repo provides.
     Optionally add a study-design figure:
     <p align="center"><img src="img/study_design.png" width="700" alt="Study design"></p>
-->

The field of spatial transcriptomics has dramatically expanded over the past few years, with the popularisation of a plethora of commercial platforms driving the growth. However, to date, no independent studies have comparatively analysed these techniques to benchmark their performance. Here, we present a well-controlled direct comparison of three spatial transcriptomics platforms. We evaluate these methods with multiple quality metrics and provide recommendations for data analysis and methods selection. This resource will help researchers to choose the optimal multiplexed in situ imaging technologies and is valuable for the development of new analysis tools.  


## Data Availability

<!-- TODO: Where the data is deposited + accession links.
     Keep or delete the table skeletons below as needed. -->

_Add a sentence on where the dataset is deposited (e.g. GEO / ArrayExpress / BioStudies) with accession links._

<!--
| Batch | Technology | Condition | Sample ID |
|-------|------------|-----------|-----------|
|       |            |           | [ID](link) |
-->


## Repository structure

```
R/              shared helper functions, sourced at the top of scripts
preprocessing/  run-once, platform-specific loading/QC/clustering (Visium, Xenium/MERSCOPE, sc/snRNA-seq)
analysis/       platform-agnostic/cross-platform analysis, writes result tables (FICTURE, ST alignment, ...)
figures/        one folder per manuscript figure, reads from results/ and data/processed/
data/           raw/ and processed/ data (gitignored — not tracked in this repo)
results/        tables/ and figures/ written by analysis/ and figures/ (gitignored)
```

See the `README.md` inside `preprocessing/` and `analysis/` for what each script does, its inputs/outputs,
and any sample-specific parameters.


## Citation

<!-- TODO: full citation once available (authors, title, journal/preprint, DOI). -->

If you use this dataset or these workflows, please cite:

Qin L*, Yip RKH*, et al. (2026). **Cross-Platform Spatial Transcriptomics Profiles of Human Breast Tissues**.


## Acknowledgements

We thank the WEHI Advanced Genomics Facility, Center for Dynamic Imaging and Advanced Histotechnology facility for supporting the generation and analysis of single-cell and spatial transcriptomics data for this project.

