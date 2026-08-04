# Purpose:  Cell-type-level ligand-receptor interaction scores (CCI): aggregate the
#           bin-level BLISA result up to (sender cell type -> receiver cell type) pairs
#           for each LR pair, per platform. Xenium and MERSCOPE only — Visium spots are
#           not single cells, so they carry no sender/receiver cell type.
#           See README.md.
# Inputs:   <results_table_dir>/seu_x_annotated.rds, seu_m_annotated.rds   from 01
#           <results_table_dir>/BLISA_x.rds, BLISA_m.rds, spot_sf.rds     from 02
#           ../../R/blisa.R                                               CCI implementation
# Outputs:  <results_table_dir>/CCI_x.rds, CCI_m.rds
#                     feed figures/Fig5_celltype_cci.R panel f

library(Seurat)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

source("../../R/blisa.R")  # runCCI_SpotAligned_MatchingWeights

ct_group <- "SingleR_labels202605"

seu_x <- readRDS(file.path(results_table_dir, "seu_x_annotated.rds"))
seu_m <- readRDS(file.path(results_table_dir, "seu_m_annotated.rds"))
BLISA_x <- readRDS(file.path(results_table_dir, "BLISA_x.rds"))
BLISA_m <- readRDS(file.path(results_table_dir, "BLISA_m.rds"))
spot_sf <- readRDS(file.path(results_table_dir, "spot_sf.rds"))
spot_ids <- rownames(spot_sf)

CCI_x <- runCCI_SpotAligned_MatchingWeights(seu = seu_x, BLISA_res = BLISA_x, spot_ids = spot_ids, ct_group = ct_group)
saveRDS(CCI_x, file.path(results_table_dir, "CCI_x.rds"))

CCI_m <- runCCI_SpotAligned_MatchingWeights(seu = seu_m, BLISA_res = BLISA_m, spot_ids = spot_ids, ct_group = ct_group)
saveRDS(CCI_m, file.path(results_table_dir, "CCI_m.rds"))
