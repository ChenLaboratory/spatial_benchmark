# ST Alignment (iST to Visium)

Aligns Xenium/MERSCOPE (iST) cell coordinates onto the matched Visium H&E image, then bins the
aligned cells into Visium spot-sized hexagons for cross-platform comparison at matched spatial units.

Run in the `stalign` conda/micromamba environment (`STalign`, `torch`, `opencv-python`; see
`/vast/projects/Spatial/lei/stalign_env`).

## 1. `01_STalign_iST_to_Visium.py` — landmark-based image registration

### Input

| Argument | Description |
|----------|-------------|
| Xenium `cells.csv.gz` / MERSCOPE `cell_metadata.csv.gz` | iST cell centroid coordinates |
| Visium `spatial/tissue_hires_image.png` | Registration target |
| `<sample>_<iST_type>_points.npy`, `<sample>_Visium_points.npy` | Manually annotated landmark pairs |

### Output

| File | Description |
|------|-------------|
| `<sample>_<iST_type>_STalign_to_Visium_tissue_hires_image_with_point_annotator.csv.gz` | Original iST cell metadata plus `aligned_x`/`aligned_y` in Visium hires-image pixel coordinates |

### Parameters (fixed across all samples and both platforms)

Rasterization at `dx=2` µm; LDDMM run with `niter=200`, `sigmaP=2e-1`, `sigmaM=sigmaB=sigmaA=0.18`,
`diffeo_start=100`, `epL=5e-11`, `epT=5e-4`, `epV=5e1` (all default/unchanged for every sample).

## 2. `02_Hexbin.py` — hexagonal binning to Visium spots

Computes hexagon-shaped bins matching the Visium spot grid geometry (from `tissue_positions.csv` /
`scalefactors_json.json`), then assigns each STalign-aligned iST cell to the Visium spot barcode
whose hexagon it falls inside.

```bash
python 02_Hexbin.py <sample> <iST_type> <aligned_fname>
```

| Argument | Description |
|----------|-------------|
| `sample` | Sample name, e.g. `ER_0114` |
| `iST_type` | `"Xenium"` or `"MERSCOPE"` |
| `aligned_fname` | Output `.csv.gz` from `01_STalign_iST_to_Visium.py` |

### Output

| File | Description |
|------|-------------|
| `<sample>_<iST_type>_hexbin.csv.gz` | iST cells with an added `Visium_spot_id` column (empty if no enclosing hexbin) |
| `<sample>_<iST_type>_hexbin_results_overview.pdf` | Diagnostic plots: hexbin centres, hexbin pattern, aligned overlay, binning result |
