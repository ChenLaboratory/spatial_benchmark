# Process Aligned Data

Combines Visium with Xenium/MERSCOPE (via the STalign+hexbin alignment in
[`../ST_alignment`](../ST_alignment)) into per-sample, cross-platform Seurat object lists, with
alignment QC metrics attached to each Visium spot and iST cell — the shared cross-platform data
structure used by later comparison analyses.

## Input

| Source | What |
|--------|------|
| `../../preprocessing/00_visium.R` output | `<sample>/so_raw.rds`, `<sample>/seu_ST_Filtered.rds` |
| `../../preprocessing/01_iST.R` output | `<platform>/<sample>/so_raw.rds`, `so.rds` |
| `../ST_alignment` output | `<sample>_<platform>_hexbin.csv.gz` (cell → Visium-spot assignment + aligned coordinates) |
| Xenium/MERSCOPE gene panel design file | used to compute the shared gene set across all three platforms |

## What it does

1. **Raw data** — for each sample and platform, joins the hexbin alignment onto the *unfiltered*
   Seurat objects (see QC metrics below). Produces `raw_data_list_AllSample.rds`.
2. **Processed data** — same joins, but on the *QC-filtered* Seurat objects, re-counting hexbin membership against the filtered cell/spot sets. Produces intermediate `processed_data_list_AllSample`.
3. **Overlapped region** — built from `processed_data_list_AllSample`: subsets each sample's Visium spots to only the hexbins that contain cells from *every* platform available for that sample, then subsets each iST object to cells falling in those overlapped hexbins. This is the common spatial region used for direct cross-platform comparison downstream. Produces `overlapped_data_list_AllSample.rds`.

## QC metrics added

Added to **Visium** metadata (both raw and processed data, `<platform>` = `Xenium`/`MERSCOPE`):

| Column | Meaning |
|--------|---------|
| `<platform>_cell_count` (raw) / `<platform>_cell_count_filtered` (processed) | Number of iST cells assigned to this hexbin/spot |
| `<platform>_total_nCount_hexbin` | Summed `nCount_RNA` across all iST cells in this hexbin (full gene panel) |
| `<platform>_avg_nCount_hexbin` | `<platform>_total_nCount_hexbin` / cell count — average transcripts per iST cell in this hexbin |
| `<platform>_total_nCount_hexbin_shared_genes` | Same sum, restricted to the cross-platform shared gene panel |
| `<platform>_avg_nCount_hexbin_shared_genes` | Same average, restricted to the shared gene panel |
| `nCount_Spatial_shared_genes` | The Visium spot's own total count, restricted to the shared gene panel |

Added to **iST** (Xenium/MERSCOPE) metadata:

| Column | Meaning |
|--------|---------|
| `Visium_spot_id` | Hexbin/Visium-spot this cell was assigned to by alignment (empty if unassigned) |
| `aligned_spatial` DimReduc | Cell coordinates in the aligned Visium image space |
| `aligned_x`, `aligned_y`, `aligned_x_in_fullres`, `aligned_y_in_fullres`, `aligned` DimReduc (overlapped data only) | Same, recomputed for the overlapped-region subset |

## Output

| File | Contents |
|------|----------|
| `raw_data_list_AllSample.rds` | `list(sample -> list(Visium, Xenium, MERSCOPE))`, unfiltered objects + QC metrics above |
| `overlapped_data_list_AllSample.rds` | Same structure, built from the QC-filtered/processed join, subset to the common overlapped spatial region |

