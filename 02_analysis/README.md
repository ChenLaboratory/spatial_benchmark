# Cross-platform benchmarking

## Analysis pipeline

### 01 — FICTURE

Guided factor decomposition of the Xenium and MERSCOPE transcripts, run for TNBC_01 and TNBC_02.
A fixed reference model built from the annotated cell types replaces FICTURE's default
unsupervised LDA, so the resulting factors correspond to known cell types rather than arbitrary
components.

| Script | Description |
|---|---|
| `<platform>/create_model_matrix.R` | Pseudobulk each annotated cell type (`SingleR_labels202605`) with edgeR, write log-CPM model matrix + colour scheme in FICTURE's format |
| `Xenium/parquet_to_csv.R` | Convert raw `transcripts.parquet` to the CSV form the workflow expects; drops control probes and incomplete rows |
| `<platform>/01_run_ficture_<platform>.smk` | Snakemake workflow: prepare input, fit the guided model, decode pixels, plot |

---

### 02 — Sensitivity metrics

| Script | Description |
|---|---|
| `01_bin_metrics_before_alignment.R` | Per-platform bin metrics computed independently, no cross-platform join |
| `02_bin_metrics_after_alignment.R` | Shared-panel bin metrics on the overlapped region |
| `03_gene_count_scatter_data.R` | Per-gene total counts per sample and platform, on the overlapped region |

**01 is the exception to the overlapped-region rule.** It reads raw, unaligned `so_raw.rds`
objects, so it is the only stage covering every dataset that exists — including the MERSCOPE-only
pilot TNBC_05 and the secondary MERSCOPE acquisitions, both filtered out downstream. 

| Output | Contents |
|---|---|
| `df_bins_summary_110um_all_ST_all_genes_before_alignment.rds` | `sample, platform, n_cells, n_transcripts, n_genes` — all genes |
| `df_bins_summary_110um_all_ST_213_shared_genes_before_alignment.rds` | Same, shared panel only |
| `bin_gene_counts_shared_genes_after_alignment.rds` | `Unique_Genes, Platform, sample` |
| `bin_nCount_shared_genes_unscaled_after_alignment.rds` | Long format; each platform's own unscaled per-hexbin totals |
| `gene_counts_per_sample.rds` | `Xenium, MERSCOPE, Visium, Visium_scaled, genes, sample` |

---

### 03 — Specificity metrics

Negative-control-based specificity for Xenium and MERSCOPE: how much signal falls on off-target
probe binding (control probes) or decoding error (control codewords), and whether that signal is
spatially random — as a true technical artefact should be — or spatially structured, which is a red
flag. 

| Script | Description |
|---|---|
| `01_summed_control_counts.R` | Per-hexbin control signal summed across all probes/codewords (Fig 3b–h) |
| `02_individual_control_moran.R` | Global Moran's I per individual control probe/codeword (Fig 3i, k; Ext Fig 6b) |
| `03_individual_control_spatial.R` | Per-hexbin counts of individual control features, all samples (Fig 3j; Ext Fig 6a) |
| `04_individual_gene_moran.R` | Global Moran's I + total counts for non-control genes (Fig 3i, k; Ext Fig 6b) |


| Output | Contents |
|---|---|
| `visium_summed_control.rds` | `list(sample -> Visium object)` with the summed-control columns added |
| `summed_control_counts_long.rds` | `Sample, Platform, NC, metric, Counts` across all samples |
| `all_moran_individual_NC.rds` | Per (sample, platform, control feature): `lisa`, `pval`, `type`, `platform2` |
| `individual_control_spatial.rds` | Nested list of per-platform 50 µm hexbin objects + control features |
| `all_moran_genes.rds` | Per (sample, platform, gene): `lisa`, `pval`, `platform2` |
| `total_counts_df.rds` | `feature_clean, total_counts, sample, platform2` — genes and controls |

---

### 04 — Segmentation metrics

Compares the three segmentation methods — Default (vendor),
[Baysor](../00_segmentation_baysor) and [Proseg v3](../00_segmentation_proseg_v3) — on QC metrics,
marker co-expression, composition and cluster separability. 

| Script | Description |
|---|---|
| `00_preliminary_analysis_seu.R` | Shared QC → cluster → SingleR annotation for the Baysor and Proseg objects; run once after segmenting |
| `01_qc_metrics.R` | Total/assigned transcripts, cell counts, cell area, transcripts per cell (Fig 4b) |
| `02_mecr.R` | Mutually exclusive co-expression rate, after QC (Fig 4c) |
| `03_cell_type_composition.R` | Proportion of cells per annotated cell type (Fig 4f) |
| `04_silhouette_width.R` | Average silhouette width by annotated cell type (Fig 4g) |


| Output | Contents |
|---|---|
| `qc_df.rds` | One row per (sample, platform, segmentation) |
| `mecr_summary.rds` | Median/mean MECR across marker pairs per dataset |
| `props_all.rds` | Cell type proportions per dataset |
| `sil_ct_results.rds` | Mean silhouette width per dataset |

---

### 05 — Cell type composition

Cell type proportions across samples and platforms after alignment — sc/snRNA-seq reference vs
Xenium vs MERSCOPE. 

| Script | Description |
|---|---|
| `01_proportion_after_alignment.R` | Cell type proportions per sample × platform, on the overlapped region |


**Output:** `comp_df.rds` — one row per (sample, platform, cell type): `n`, `prop`.

---

### 06 — Cell-cell communication

Ligand-receptor interaction analysis for **TNBC_02**, with pathologist-annotated tissue regions
overlaid. 

| Script | Description |
|---|---|
| `01_prepare_annotated_data.R` | Cross-platform objects + pathology region/type annotation |
| `02_run_blisa.R` | Bin-level LR interaction testing (BLISA), per platform |
| `03_run_cci.R` | Cell-type-level LR interaction scores (CCI) |
| `04_region_lr_proportion.R` | Proportion of significant bins per pathology region × LR pair |


| Output | Contents |
|---|---|
| `seu_{x,m,v}_annotated.rds`, `annotation_df.rds` | Annotated per-platform objects; cell-level QuPath table for region ordering |
| `BLISA_{v,x,m}.rds`, `spot_sf.rds`, `d_theoretical.rds` | Bin-level results; Visium spot geometry; theoretical spot spacing |
| `CCI_x.rds`, `CCI_m.rds` | Cell-type-pair interaction scores |
| `prop_all.rds` | Significant-bin proportion per region × LR pair |
