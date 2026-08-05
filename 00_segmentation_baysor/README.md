# Baysor segmentation


## Pipeline

Run **once per sample per platform** — there is no batch loop. Fill in the placeholder variables
at the top of each script (or copy it per sample) and run it from that sample's own working
directory.

```
01_run_baysor_xenium.sh       segmentation (Xenium)
01_run_baysor_merscope.sh     segmentation (MERSCOPE)
02_run_baysor_to_seu.sh       segmentation output -> Seurat object (calls generate_seu.R)
generate_seu.R                Seurat-conversion logic, shared by both platforms
map_transcripts_baysor_xenium.py / map_transcripts_baysor_merscope.py   utility, called by 02
xenium.toml / merscope.toml   Baysor config, one per platform
```

Numbered scripts are the ones to run.

1. **Segmentation** — `01_run_baysor_<platform>.sh` on that sample's raw transcripts.
2. **Seurat conversion** — `02_run_baysor_to_seu.sh` in the same working directory, giving
   `<sample>_so_raw.rds`.
3. **QC / clustering / annotation** —
   [`../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R`](../02_analysis/04_segmentation_metrics/00_preliminary_analysis_seu.R),
   shared with Proseg. It loops over both methods and all samples, so run it once after
   segmenting everything rather than per sample.

**Output:** `<processed_data_dir>/segmentation/baysor/<platform>/<sample>_so.rds`.

**Sample coverage:** Xenium for TNBC_01–04 and ER_01/ER_02. MERSCOPE for the five
benchmarking-cohort datasets only — TNBC_01 and TNBC_02 on the V1 instrument, TNBC_03, ER_01 and
ER_02 on V2. TNBC_04's MERSCOPE run is excluded for poor quality.
