# Bin-level ligand-receptor (LR) interaction testing (BLISA) at the Visium
# spot/hexbin level, for each platform independently (Visium native spots;
# Xenium/MERSCOPE cells aggregated to their assigned Visium spot). Feeds
# 03_run_cci.R, 04_region_lr_proportion.R and Fig5 panels d/e. See README.md.

library(Seurat)
library(sf)
library(FNN)
library(CellChat)

sample <- "MH0026"

shared_genes_file <- "/path/to/analysis/process_aligned_data/shared_genes_withVisium.rds"

seu_x <- readRDS("seu_x_annotated.rds")
seu_m <- readRDS("seu_m_annotated.rds")
seu_v <- readRDS("seu_v_annotated.rds")
shared_genes <- readRDS(shared_genes_file)

source("runBLISA.R")  # runBLISA.default.isolates.removed, filterLRpairs, aggregate_cells_to_spots, visium_spot_distance

## ---- Visium spot geometry ---------------------------------------------------

coords <- GetTissueCoordinates(seu_v)
coords_df <- data.frame(cell_id = rownames(coords), x_centroid = coords[, "imagecol"],
                         y_centroid = coords[, "imagerow"], row.names = rownames(coords))
spot_ids <- coords_df$cell_id
spot_sf  <- sf::st_as_sf(coords, coords = c("imagerow", "imagecol"), crs = NA)
rownames(spot_sf) <- spot_ids

xen_cnt <- table(seu_x$Visium_spot_id)
mer_cnt <- table(seu_m$Visium_spot_id)
spot_sf$Xenium_cell_count_filtered   <- as.integer(xen_cnt[rownames(spot_sf)])
spot_sf$MERSCOPE_cell_count_filtered <- as.integer(mer_cnt[rownames(spot_sf)])
spot_sf$Xenium_cell_count_filtered[is.na(spot_sf$Xenium_cell_count_filtered)]     <- 0L
spot_sf$MERSCOPE_cell_count_filtered[is.na(spot_sf$MERSCOPE_cell_count_filtered)] <- 0L

d_theoretical <- visium_spot_distance(seu_v, "lowres")

## ---- BLISA, one run per platform -------------------------------------------

counts_v_spot <- GetAssayData(seu_v, assay = "Spatial", slot = "counts")
counts_v_spot <- counts_v_spot[intersect(rownames(counts_v_spot), shared_genes), , drop = FALSE]
LR_v <- filterLRpairs(counts = counts_v_spot, min_ligand = 1, min_receptor = 1, LR_df = CellChatDB.human$interaction)
BLISA_v <- runBLISA.default.isolates.removed(
  counts_matrix = counts_v_spot, bin_sf = spot_sf, LR_df = LR_v,
  dmax = d_theoretical * 2, hex_size = d_theoretical, n_cells_col = NA
)
saveRDS(BLISA_v, "BLISA_v.rds")

counts_x_spot <- aggregate_cells_to_spots(seu_x, shared_genes, spot_ids_keep = spot_ids)
LR_x <- filterLRpairs(counts = counts_x_spot, min_ligand = 1, min_receptor = 1, LR_df = CellChatDB.human$interaction)
BLISA_x <- runBLISA.default.isolates.removed(
  counts_matrix = counts_x_spot, bin_sf = spot_sf, LR_df = LR_x,
  dmax = d_theoretical * 2, hex_size = d_theoretical, n_cells_col = "Xenium_cell_count_filtered"
)
saveRDS(BLISA_x, "BLISA_x.rds")

counts_m_spot <- aggregate_cells_to_spots(seu_m, shared_genes, spot_ids_keep = spot_ids)
LR_m <- filterLRpairs(counts = counts_m_spot, min_ligand = 1, min_receptor = 1, LR_df = CellChatDB.human$interaction)
BLISA_m <- runBLISA.default.isolates.removed(
  counts_matrix = counts_m_spot, bin_sf = spot_sf, LR_df = LR_m,
  dmax = d_theoretical * 2, hex_size = d_theoretical, n_cells_col = "MERSCOPE_cell_count_filtered"
)
saveRDS(BLISA_m, "BLISA_m.rds")

saveRDS(spot_sf, "spot_sf.rds")
saveRDS(d_theoretical, "d_theoretical.rds")
