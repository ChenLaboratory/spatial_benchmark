# Proseg (v3) segmentation pipeline

Run **once per sample per platform** — there is no batch loop. Fill in the placeholder variables
at the top of each script (or copy it per sample) and run it from that sample's own working
directory.

```
01_run_proseg_xenium.sh       segmentation (Xenium)
01_run_proseg_merscope.sh     segmentation (MERSCOPE)
02_run_proseg_to_seu.sh       segmentation output -> Seurat object (calls generate_seu.R)
03_run_cell_area_2d.sh        optional: 2D cell-area metric (calls zarr_cell_area_2d.py)
generate_seu.R                Seurat-conversion logic
zarr_to_h5ad.py / zarr_cell_area_2d.py   utility, called by 02 / 03
```

Numbered scripts are the ones to run.

1. **Segmentation** — `01_run_proseg_<platform>.sh` on that sample's raw transcripts.
2. **Seurat conversion** — `02_run_proseg_to_seu.sh` in the same working directory, giving
   `<sample>_so_raw.rds`.
3. **2D cell area** *(optional)* — `03_run_cell_area_2d.sh` computes a 2D cell-area metric from
   Proseg's boundary polygons, since Proseg's default output is a 3D volume/surface area. It
   writes `<sample>_proseg_xy_area.csv`, which
   [`../02_analysis/04_segmentation_metrics/01_qc_metrics.R`](../02_analysis/04_segmentation_metrics/01_qc_metrics.R)
   reads for the median-cell-area panel; skip it and that panel has no Proseg values.
4. **QC / clustering / annotation** —
   [`../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R`](../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R),
   shared with Baysor. It loops over both methods and all samples, so run it once after
   segmenting everything rather than per sample.

**Output:** `<processed_data_dir>/segmentation/proseg_v3/<platform>/<sample>_so.rds`, with
`processed_data_dir` set in [`../config.R`](../config.R).

**Sample coverage:** Xenium for TNBC_01–04 and ER_01/ER_02. MERSCOPE for the five
benchmarking-cohort datasets only — TNBC_01 and TNBC_02 on the V1 instrument, TNBC_03, ER_01 and
ER_02 on V2. TNBC_04's MERSCOPE run is excluded for poor quality.