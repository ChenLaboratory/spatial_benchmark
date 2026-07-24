# Cell Type Composition

Cross-sample, cross-platform cell type proportions (sc/snRNA-seq reference vs Xenium vs MERSCOPE),
after cross-platform alignment. Feeds [`../../figures/Fig5.R`](../../figures/Fig5.R) panel b.

```
01_proportion_after_alignment.R   cell type proportions per sample x platform, after alignment
```

For each of the 6 benchmarking samples: sc/snRNA-seq reference cell type proportions (labels
already assigned, "mixture" cells excluded) shown as-is; Xenium/MERSCOPE proportions restricted to
cells retained in the cross-platform overlapped hexbins (`../process_aligned_data`) — i.e. the same
cells used throughout the rest of the cross-platform comparisons in this repo, so the composition
comparison is on a matched spatial footprint across platforms.

**Output**: `comp_df.rds` — one row per (sample, platform, cell type): `n`, `prop`.

## Not included

The source analysis (`cell_type_identification.Rmd`) also computes the equivalent "before
alignment" composition (all QC-passed cells, no cross-platform overlap restriction) — not
reproduced here since Fig5b only shows the after-alignment comparison.
