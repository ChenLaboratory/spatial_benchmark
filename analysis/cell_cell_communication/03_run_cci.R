# Cell-type-level ligand-receptor interaction scores (CCI): aggregates the
# bin-level BLISA result up to (sender cell type -> receiver cell type) pairs
# for each LR pair, per platform. Feeds Fig5 panel f. See README.md.

library(Seurat)

source("runBLISA.R")  # runCCI_SpotAligned_MatchingWeights

ct_group <- "SingleR_labels202605"

seu_x <- readRDS("seu_x_annotated.rds")
seu_m <- readRDS("seu_m_annotated.rds")
BLISA_x <- readRDS("BLISA_x.rds")
BLISA_m <- readRDS("BLISA_m.rds")
spot_sf <- readRDS("spot_sf.rds")
spot_ids <- rownames(spot_sf)

CCI_x <- runCCI_SpotAligned_MatchingWeights(seu = seu_x, BLISA_res = BLISA_x, spot_ids = spot_ids, ct_group = ct_group)
saveRDS(CCI_x, "CCI_x.rds")

CCI_m <- runCCI_SpotAligned_MatchingWeights(seu = seu_m, BLISA_res = BLISA_m, spot_ids = spot_ids, ct_group = ct_group)
saveRDS(CCI_m, "CCI_m.rds")
