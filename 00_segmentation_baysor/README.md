# Baysor segmentation

Alternative cell segmentation for Xenium and MERSCOPE. Takes the same raw transcripts as
[`../01_preprocessing/iST`](../01_preprocessing/iST) and produces the same artifact, so the three
segmentations are directly comparable:

| Method | Where it lives |
|--------|----------------|
| **Default** (vendor) | [`../01_preprocessing/iST/01_preprocess_iST.R`](../01_preprocessing/iST/01_preprocess_iST.R) |
| **Baysor** | this folder |
| **Proseg v3** | [`../00_segmentation_proseg_v3`](../00_segmentation_proseg_v3) |

The metrics comparing them are in
[`../02_analysis/04_segmentation_metrics`](../02_analysis/04_segmentation_metrics), feeding
[`../03_figures/Fig4_segmentation.R`](../03_figures/Fig4_segmentation.R).

```
01_run_baysor_xenium.sh       segmentation (Xenium)
01_run_baysor_merscope.sh     segmentation (MERSCOPE)
02_run_baysor_to_seu.sh       segmentation output -> Seurat object (calls generate_seu.R)
generate_seu.R                Seurat-conversion logic, shared by both platforms
map_transcripts_baysor_xenium.py / map_transcripts_baysor_merscope.py   utility, called by 02
xenium.toml / merscope.toml   Baysor config, one per platform
```

Numbered scripts are the ones you run; the `.py` and `generate_seu.R` files are called by them,
not run directly.

## Pipeline

Run **once per sample per platform** — there is no batch loop. Fill in the placeholder variables
at the top of each script (or copy it per sample) and run it from that sample's own working
directory.

1. **Segmentation** — `01_run_baysor_<platform>.sh` on that sample's raw transcripts.
2. **Seurat conversion** — `02_run_baysor_to_seu.sh` in the same working directory, giving
   `<sample>_so_raw.rds`.
3. **QC / clustering / annotation** —
   [`../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R`](../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R),
   shared with Proseg. It loops over both methods and all samples, so run it once after
   segmenting everything rather than per sample.

**Output:** `<processed_data_dir>/segmentation/baysor/<platform>/<sample>_so.rds`, with
`processed_data_dir` set in [`../config.R`](../config.R).

## Notes

- **Sample coverage:** Xenium for TNBC_01–04 and ER_01/ER_02; MERSCOPE for the same set minus
  TNBC_04, whose MERSCOPE run is not part of the benchmarking cohort. Sample IDs are the
  manuscript IDs — see the [root README](../README.md#samples).
- Step 3 duplicates the SingleR transfer in
  [`../01_preprocessing/iST/02_annotate_iST.R`](../01_preprocessing/iST/02_annotate_iST.R),
  deliberately, so each arm runs end to end independently. Change the annotation procedure in
  both or neither.
- **Viewer-import steps are intentionally excluded.** The original workflow also re-imported the
  segmentation into each platform's native viewer format (XeniumRanger for Xenium, VPT/`.vzg` for
  MERSCOPE). Neither feeds the analysis, so only the steps producing the analysed Seurat object
  are included here.
