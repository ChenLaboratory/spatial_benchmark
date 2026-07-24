# Specificity Metrics

Negative-control-based measures of gene detection specificity for Xenium and MERSCOPE — how much
signal falls on off-target probe binding (negative control probes) or decoding error (negative
control codewords), and whether that signal is spatially random (as it should be for a true
technical artefact) or spatially structured (a red flag). Feeds
[`../../figures/Fig3`](../../figures/Fig3) panels b-k and
[`../../figures/figS6_gene_specificity`](../../figures/figS6_gene_specificity).

```
01_summed_control_counts.R        per-hexbin NC signal, summed across all probes/codewords (Fig3 b-h)
02_individual_control_moran.R     Global Moran's I per individual NC probe/codeword (Fig3 i, k; figS6 b)
03_individual_control_spatial.R   per-hexbin counts of individual NC features, all samples (Fig3 j; figS6 a)
04_individual_gene_moran.R        Global Moran's I + total counts for non-control genes (Fig3 i, k; figS6 b)
```

All scripts use `overlapped_data_list_AllSample.rds` from
[`../process_aligned_data`](../process_aligned_data) (the QC-filtered, cross-platform overlapped
region, 110um Visium-matched hexbins). "TNBC-02" in the manuscript = sample `MH0026` — the only
sample with both Xenium and MERSCOPE V1 in the aligned data (used for Fig3 panels b, j, k; the
other 5 samples are shown in figS6 panel a).

Xenium's negative control probes/codewords are excluded from the Xenium Seurat object's expression
matrix in this pipeline; wherever their counts are needed, they're pulled directly from the raw
`transcripts.parquet` files. MERSCOPE's negative controls remain as rows in its count matrix,
identified by name pattern (`Blank` = codeword; `Sabrina`/`FC[0-9]+` = probe, depending on
MERSCOPE hardware version). Xenium's per-cell `control_probe_counts`/`control_codeword_counts` are
native 10x Xenium QC metadata columns, expected already present from preprocessing.

## 01 — summed control counts

For each sample, sums all NC probe counts and all NC codeword counts per cell, then per Visium
hexbin, for each platform. Derives per-hexbin `total`, `PerCell` (per cell in that hexbin),
`percent` (of platform's total hexbin counts), and `FDR` (false discovery rate: percent scaled by
the ratio of gene panel size to number of NC probes/codewords).

**Output**:
| File | Contents |
|------|----------|
| `visium_summed_control.rds` | `list(sample -> Visium Seurat object)` with the above columns added; feeds Fig3b spatial plots |
| `summed_control_counts_long.rds` | Long-format table (`Sample, Platform, NC, metric, Counts`) across all samples/platforms/params; feeds Fig3c-h |

## 02 — individual control Moran's I

Adds per-hexbin counts for each *individual* NC probe/codeword (not summed), then computes Global
Moran's I (`scider::globalMoran`, k=10 spatial neighbors, 9999 permutations) for each one, per
sample and platform.

**Output**: `all_moran_individual_NC.rds` — one row per (sample, platform, NC feature): `lisa`,
`pval`, `type` (probe/codeword), `platform2` (MERSCOPE_V1/V2/Xenium).

## 03 — individual control spatial (all samples)

Same idea as 02 but at each platform's own finer 50um hexbins (not the 110um Visium-matched bins),
since this feeds spatial visualizations (Fig3j, figS6a) rather than a cross-sample summary
statistic. Run for all 6 samples (Xenium) / 5 samples with MERSCOPE.

**Output**: `individual_control_spatial.rds` — `list(sample -> list(Xenium = list(spe_hex,
neg_controls), MERSCOPE = list(spe_hex, neg_controls)))` (MERSCOPE absent for samples without it).

## 04 — individual gene Moran's I (non-control genes)

Same Moran's I method as 02, applied to every non-control gene instead of negative controls, plus
total transcript counts per gene/probe/codeword per sample and platform (the x-axis for Fig3k /
figS6b).

**Output**:
| File | Contents |
|------|----------|
| `all_moran_genes.rds` | One row per (sample, platform, gene): `lisa`, `pval`, `platform2` |
| `total_counts_df.rds` | `feature_clean, total_counts, sample, platform2` — genes + NC probes/codewords |

## Not included

The source analyses also computed: Global Moran's I on the *summed* (not individual) NC signal (a
barplot variant not used in the manuscript panels); cell-level (not bin-level) NC violin plots with
matched scRNA-seq; extreme/outlier Moran's I gene tables and labeled scatter variants; low-count
random-gene and platform-discordant-gene analyses; and the all-negative-controls (not
shared-gene-restricted) and n-labeled variants of the count-range threshold scatter. None of these
feed Fig3 panels b-k or figS6.
