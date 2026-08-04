# Integration

Stage 2: bring the platforms into a **shared coordinate space**. Unlike
[`../01_preprocessing`](../01_preprocessing), which handles one platform at a time, everything here
operates on several platforms at once, and everything in [`../02_analysis`](../02_analysis) that compares
platforms depends on its output.

```
01_STalign_iST_to_Visium.py   register Xenium/MERSCOPE cells onto the matched Visium H&E image
02_hexbin.py                  assign aligned cells to Visium spot-sized hexagons
03_build_aligned_objects.R    join the alignment onto the Seurat objects; build the overlapped region
```

Run 01 and 02 in the `stalign` conda/micromamba environment (`STalign`, `torch`,
`opencv-python`). Set the paths in the block at the top of each script — see the
[root README](../README.md#workflows).

Sample IDs below are the manuscript IDs, which are also the GEO IDs and the IDs used in filenames;
they are listed in the [root README](../README.md#samples).

---

## `01_STalign_iST_to_Visium.py` — landmark-based image registration

Registers iST cell coordinates onto the matched Visium H&E image with STalign (landmark-based
LDDMM).

**Input**

| Argument | Description |
|----------|-------------|
| Xenium `cells.csv.gz` / MERSCOPE `cell_metadata.csv.gz` | iST cell centroid coordinates |
| Visium `spatial/tissue_hires_image.png` | Registration target |
| `<sample>_<iST_type>_points.npy`, `<sample>_Visium_points.npy` | Manually annotated landmark pairs |

**Output**

| File | Description |
|------|-------------|
| `<sample>_<iST_type>_STalign_to_Visium_tissue_hires_image_with_point_annotator.csv.gz` | Original iST cell metadata plus `aligned_x`/`aligned_y` in Visium hires-image pixel coordinates |

**Parameters (fixed across all samples and both platforms):** rasterization at `dx=2` µm; LDDMM with
`niter=200`, `sigmaP=2e-1`, `sigmaM=sigmaB=sigmaA=0.18`, `diffeo_start=100`, `epL=5e-11`,
`epT=5e-4`, `epV=5e1` — all default/unchanged for every sample.

---

## `02_hexbin.py` — hexagonal binning to Visium spots

Computes hexagon-shaped bins matching the Visium spot grid geometry (from `tissue_positions.csv` /
`scalefactors_json.json`), then assigns each aligned iST cell to the Visium spot barcode whose
hexagon it falls inside. The 110 µm hexagon matches the Visium spot centre-to-centre distance, so
all three platforms can be compared at matched spatial units.

```bash
python 02_hexbin.py <sample> <iST_type> <aligned_fname>
```

| Argument | Description |
|----------|-------------|
| `sample` | Internal sample name, e.g. `ER_02` |
| `iST_type` | `"Xenium"` or `"MERSCOPE"` |
| `aligned_fname` | Output `.csv.gz` from `01_STalign_iST_to_Visium.py` |

**Output**

| File | Description |
|------|-------------|
| `<sample>_<iST_type>_hexbin.csv.gz` | iST cells with an added `Visium_spot_id` column (empty if no enclosing hexbin) |
| `<sample>_<iST_type>_hexbin_results_overview.pdf` | Diagnostic plots: hexbin centres, hexbin pattern, aligned overlay, binning result |

---

## `03_build_aligned_objects.R` — cross-platform object lists

Combines Visium with Xenium/MERSCOPE into per-sample, cross-platform Seurat object lists, with
alignment QC metrics attached to each Visium spot and iST cell.

**Input**

| Source | What |
|--------|------|
| [`../01_preprocessing/visium`](../01_preprocessing/visium) output | `<sample>/so_raw.rds`, `<sample>/seu_ST_Filtered.rds` |
| [`../01_preprocessing/iST`](../01_preprocessing/iST) output | `<platform>/<sample>/so_raw.rds`, `so.rds` |
| `02_hexbin.py` output | `<sample>_<platform>_hexbin.csv.gz` |
| Xenium/MERSCOPE gene panel design files | used to compute the shared gene set across all three platforms |

**What it does**

1. **Raw data** — for each sample and platform, joins the hexbin alignment onto the *unfiltered*
   Seurat objects. Produces `raw_data_list_AllSample.rds`.
2. **Processed data** — the same joins on the *QC-filtered* objects, re-counting hexbin membership
   against the filtered cell/spot sets. Produces intermediate `processed_data_list_AllSample`.
3. **Overlapped region** — subsets each sample's Visium spots to only the hexbins containing cells
   from *every* platform available for that sample, then subsets each iST object to cells in those
   hexbins. This is the common spatial region used for direct cross-platform comparison downstream.
   Produces `overlapped_data_list_AllSample.rds`.

**QC metrics added to Visium metadata** (both raw and processed; `<platform>` = `Xenium`/`MERSCOPE`):

| Column | Meaning |
|--------|---------|
| `<platform>_cell_count` (raw) / `<platform>_cell_count_filtered` (processed) | Number of iST cells assigned to this hexbin/spot |
| `<platform>_total_nCount_hexbin` | Summed `nCount_RNA` across all iST cells in this hexbin (full gene panel) |
| `<platform>_avg_nCount_hexbin` | The above / cell count — average transcripts per iST cell in this hexbin |
| `<platform>_total_nCount_hexbin_shared_genes` | Same sum, restricted to the cross-platform shared gene panel |
| `<platform>_avg_nCount_hexbin_shared_genes` | Same average, restricted to the shared gene panel |
| `nCount_Spatial_shared_genes` | The Visium spot's own total count, restricted to the shared gene panel |

**QC metrics added to iST (Xenium/MERSCOPE) metadata:**

| Column | Meaning |
|--------|---------|
| `Visium_spot_id` | Hexbin/Visium spot this cell was assigned to (empty if unassigned) |
| `aligned_spatial` DimReduc | Cell coordinates in the aligned Visium image space |
| `aligned_x`, `aligned_y`, `aligned_x_in_fullres`, `aligned_y_in_fullres`, `aligned` DimReduc | Same, recomputed for the overlapped-region subset (overlapped data only) |

**Output**

| File | Contents |
|------|----------|
| `raw_data_list_AllSample.rds` | `list(sample -> list(Visium, Xenium, MERSCOPE))`, unfiltered objects + QC metrics above |
| `overlapped_data_list_AllSample.rds` | Same structure, built from the QC-filtered join, subset to the common overlapped spatial region |
| `shared_genes_iST.rds` | Genes shared by the Xenium and MERSCOPE panels |
| `shared_genes_withVisium.rds` | The 213-gene panel shared across all three platforms |
