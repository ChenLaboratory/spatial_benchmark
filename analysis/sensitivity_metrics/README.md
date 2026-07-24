# Sensitivity Metrics

Per-sample, per-platform bin-level transcript and gene detection metrics used to compare
sequencing depth and gene detection sensitivity across Xenium, MERSCOPE and Visium (the source
data for [`../../figures/Fig2`](../../figures/Fig2)). Xenium/MERSCOPE cells are binned into 110um
hexagons to match the Visium spot centre-to-centre distance; Visium's own native spots are used
directly as bins.

```
01_bin_metrics_before_alignment.R   per-platform bin metrics, computed independently (no cross-platform join)
02_bin_metrics_after_alignment.R    shared-panel bin metrics on the cross-platform overlapped region
03_gene_count_scatter_data.R        per-gene total counts per sample/platform, on the overlapped region
```

## 01 — before alignment

For each sample/platform, loads the raw (unfiltered, unaligned) Seurat object, bins iST cells into
110um hexagons with `scider::gridDensity`/`gridSPE`, and computes per-bin cell count, transcript
count and gene count. Run twice: once on all genes, once restricted to the shared gene panel
(`shared_genes_withVisium.rds` from [`../process_aligned_data`](../process_aligned_data)).

MERSCOPE control probes (`Blank`, `Sabrina`, `FC[0-9]+`) are removed before computing any metric.
Covers the 6-sample benchmarking cohort plus the pilot sample 24MH0042 (present in the source data
for a few platforms via extra `v1_so_raw.rds`/`v2_so_raw.rds` runs); 24MH0042 is not part of the
benchmarking cohort and is filtered out downstream in the figure scripts.

**Input**: per-sample `so_raw.rds` (and, where available, extra `v1_so_raw.rds`/`v2_so_raw.rds`
MERSCOPE runs) from `../../preprocessing`; `shared_genes_withVisium.rds` from
`../process_aligned_data`.

**Output**:
| File | Contents |
|------|----------|
| `df_bins_summary_110um_all_ST_all_genes_before_alignment.rds` | `sample, platform, n_cells, n_transcripts, n_genes` — all genes |
| `df_bins_summary_110um_all_ST_213_shared_genes_before_alignment.rds` | Same columns, restricted to the shared gene panel |

## 02 — after alignment

Uses `overlapped_data_list_AllSample.rds` from `../process_aligned_data` (the QC-filtered,
cross-platform overlapped region, benchmarking cohort only) to compute, for the shared gene panel
only:
- detected genes per hexbin (iST cells aggregated by assigned Visium hexbin, then counted for
  non-zero shared-panel genes; Visium counted directly from its own spot)
- total shared-panel transcript counts per hexbin, using each platform's own (unscaled) values —
  i.e. no area-based rescaling of the Visium counts.

**Input**: `overlapped_data_list_AllSample.rds`, `shared_genes_withVisium.rds` from
`../process_aligned_data`.

**Output**:
| File | Contents |
|------|----------|
| `bin_gene_counts_shared_genes_after_alignment.rds` | `Unique_Genes, Platform, sample` |
| `bin_nCount_shared_genes_unscaled_after_alignment.rds` | `orig.ident, variable, value, sample, platform` (long format) |

## 03 — gene count scatter data

For each of the 6 benchmarking-cohort samples, sums transcript counts per gene across all
cells/spots (Xenium, MERSCOPE, Visium raw and area-scaled to hexbin size, `* 402/181`), on the
same cross-platform overlapped region as 02, restricted to the shared gene panel present in that
sample's data. Feeds the gene count scatter panels in
[`../../figures/Fig2`](../../figures/Fig2) (f2h) and
[`../../figures/figS5_gene_scatter`](../../figures/figS5_gene_scatter) (Extended Data Fig. 5).

**Input**: `overlapped_data_list_AllSample.rds`, `shared_genes_withVisium.rds` from
`../process_aligned_data`.

**Output**: `gene_counts_per_sample.rds` — `Xenium, MERSCOPE, Visium, Visium_scaled, genes, sample`.

## Not included

The source sensitivity analysis also computed cell-level (not bin-level) QC comparisons, all-genes
bin metrics after alignment, MERSCOPE-vs-Xenium per-gene fold-sensitivity, and outlier-gene
inspection (cell-/bin-level spatial plots, CCC) for the pilot sample 24MH0042 — none of these feed
the manuscript figures and are not reproduced here.
