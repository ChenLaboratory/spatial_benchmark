# Cell-Cell Communication

Ligand-receptor (LR) interaction analysis for sample MH0026 ("TNBC-02" in the manuscript), across
its full Visium tissue region (not restricted to the cross-platform-overlapped hexbins used
elsewhere in this repo), with pathologist-annotated tissue regions overlaid. Feeds `figures/Fig5.R`
panels a, d, e, f.

```
01_prepare_annotated_data.R   cross-platform data + pathology region/type annotation
02_run_blisa.R                bin-level LR interaction testing (BLISA), per platform
03_run_cci.R                  cell-type-level LR interaction scores (CCI)
04_region_lr_proportion.R     proportion of significant bins per pathology region x LR pair
runBLISA.R                    BLISA/CCI method library, sourced by 02 and 03 (not run directly)
```

## BLISA/CCI method library

`runBLISA.R` is the underlying statistical method, `source()`-d by `02_run_blisa.R` and
`03_run_cci.R`: `filterLRpairs` (LR pair QC), `runBLISA.default.isolates.removed` (bin-level
spatial LR interaction testing — a bivariate local Moran's I on ligand/receptor expression between
neighbouring bins, with isolated/low-cell bins excluded from the statistic rather than given
arbitrary neighbours), and `runCCI_SpotAligned_MatchingWeights` (aggregates bin-level hotspots to
cell-type-pair interaction scores). Also included: `hex_binning_cells`/`aggregate_cells_to_spots`
(cell-to-bin aggregation), `visium_spot_distance`, and the `CCIspatial`/`CCIheatmap*`/`plotLRIsum`/
`plot_LRI` visualization helpers used elsewhere in the source analysis (not all of which are used
by Fig5 — see "Not included" below).

## 01 — prepare annotated data

Loads raw (unfiltered) cross-platform data with alignment metadata already attached (from
`../process_aligned_data`), intersects with each platform's independent QC-filtered object, then
restricts Xenium/MERSCOPE to cells whose assigned Visium spot passed Visium QC — giving every
iST cell in the Visium tissue footprint, not just the benchmarking cross-platform overlap.
Cell type labels (`SingleR_labels202605`) are transferred from the QC-filtered objects (not
present on the raw ones). Pathology annotation is read from per-region QuPath cell-export CSVs:
Xenium cells are matched by cell ID directly; Visium spots take the majority-vote annotation of
their assigned (annotated) Xenium cells; MERSCOPE cells inherit their assigned Visium spot's
annotation. `pathology_region` = specific annotated region (e.g. one replicate of a tissue type);
`pathology_type` = region with the replicate suffix stripped.

**Output**: `seu_x_annotated.rds`, `seu_m_annotated.rds`, `seu_v_annotated.rds`,
`annotation_df.rds` (cell-level QuPath table, used for region/type ordering in figures).

## 02 — BLISA

Bin-level (Visium spot / hexbin) LR interaction testing, one run per platform: Visium spots used
directly; Xenium/MERSCOPE cells aggregated to their assigned spot first. LR pairs restricted to the
shared gene panel and filtered for minimum ligand/receptor expression (CellChatDB.human).

**Output**: `BLISA_v.rds`, `BLISA_x.rds`, `BLISA_m.rds`, `spot_sf.rds` (Visium spot geometry +
per-spot iST cell counts), `d_theoretical.rds` (theoretical Visium spot spacing, used as the
neighbourhood distance).

## 03 — CCI

Aggregates each platform's bin-level BLISA result to (sender cell type -> receiver cell type)
interaction scores per LR pair.

**Output**: `CCI_x.rds`, `CCI_m.rds`.

## 04 — region x LR-pair proportion

For each (platform, LR pair, pathology region), the proportion of that region's bins that were
BLISA-significant for that LR pair.

**Output**: `prop_all.rds`.

## Not included

The source analysis also computed: MERSCOPE-annotated / Xenium-annotated overview panels per
pathology sub-type (not used in the manuscript panels); several alternate spatial p-value plot
styles (unfilled, all-bins, region-bordered instead of type-bordered); the full genome-wide LRI
spatial PDF and per-LR CCI heatmap PDF loops (commented out / disabled in the source itself);
type x LR-pair (rather than region x LR-pair) proportion heatmaps; and the "all CellChat LR pairs
in shared genes" extended heatmap variants (with/without Xenium-order clustering). None of these
feed Fig5 panels a, d, e, f.
