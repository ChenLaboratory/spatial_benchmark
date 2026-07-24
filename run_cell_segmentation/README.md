# Cell Segmentation

Compares two cell segmentation methods for Xenium/MERSCOPE data — Baysor and proseg (v3). Each
method's scripts are run **once per sample** (same steps, repeated individually for every
sample/platform combination) — there is no batch loop; fill in the placeholder variables at the
top of each script (or copy it per sample) and run it from within that sample's own working
directory. Both methods end at a Seurat object, which then goes through a shared QC/clustering/
SingleR-annotation step so results are comparable across methods.

```
baysor/
  01_run_baysor_xenium.sh     segmentation (Xenium)
  01_run_baysor_merscope.sh   segmentation (MERSCOPE)
  02_run_baysor_to_seu.sh     segmentation output -> Seurat object (calls generate_seu.R)
  generate_seu.R              Seurat-conversion logic, shared by both platforms
  map_transcripts_baysor_xenium.py / map_transcripts_baysor_merscope.py   utility, called by 02
  xenium.toml / merscope.toml Baysor config, one per platform

proseg_v3/
  01_run_proseg_xenium.sh     segmentation (Xenium)
  01_run_proseg_merscope.sh   segmentation (MERSCOPE)
  02_run_proseg_to_seu.sh     segmentation output -> Seurat object (calls generate_seu.R)
  generate_seu.R              Seurat-conversion logic
  03_run_cell_area_2d.sh      optional: 2D cell-area metric (calls zarr_cell_area_2d.py)
  zarr_to_h5ad.py / zarr_cell_area_2d.py   utility scripts, called by 02 / 03

preliminary_analysis_seu/
  check_seu.R                 shared QC + clustering + SingleR annotation, one sample at a time
```

Numbered scripts are the ones you actually run; the unnumbered `.py`/`generate_seu.R` files are
called by them, not run directly.

## Pipeline

For each sample, for each method:

1. **Segmentation** — run the platform-specific `01_run_*.sh` on that sample's raw transcripts.
2. **Seurat conversion** — run `02_run_*_to_seu.sh` in the same working directory; it builds a
   Seurat object from the segmentation output.
3. **QC / clustering / annotation** — run `preliminary_analysis_seu/check_seu.R` on that Seurat
   object: removes control features, applies a fixed QC cutoff, runs the standard Seurat pipeline
   (normalize → scale → PCA → UMAP → cluster), then transfers cell type labels from the matched
   sc/snRNA-seq reference via SingleR.

proseg_v3 additionally has an optional step 3, `03_run_cell_area_2d.sh`, which computes a 2D
cell-area metric from proseg's boundary polygons (proseg's default is a 3D volume/surface area).

## Notes

- **Sample coverage**: MERSCOPE was not run for one sample (TN0554) with either method; both
  methods otherwise cover the same sample set.
- **Segmentation-import steps are intentionally excluded.** Both methods originally also
  re-imported their segmentation into the platform's native viewer format (XeniumRanger for
  Xenium, VPT/`.vzg` for MERSCOPE) — confirmed unused by the actual downstream analysis, so not
  included here; only the steps that feed the analyzed Seurat object are.
- A legacy, superseded version of proseg (v1/v2) and a one-off debugging notebook from the source
  material are likewise not included, for the same reason: neither fed the analysis actually used.
