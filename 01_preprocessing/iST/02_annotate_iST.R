# Purpose:  Annotate cell types in one Xenium or MERSCOPE sample by SingleR label
#           transfer from the matched sc/snRNA-seq reference, using the manual
#           `cell_type202605` labels from ../scRNA/03_annotate_cell_types.R. Reference
#           cells labelled "mixture" (unresolved sub-clustering) are dropped first.
#           Run once per sample per platform; `sample_id`, `platform` and `ref_type`
#           are set by hand at the top — see the table in README.md.
# Inputs:   <processed_data_dir>/<platform>/<sample_id>/so.rds   clustered object from 01
#           <processed_data_dir>/scRNA/<sample_id>/seu_<sample_id>.rds
#                     annotated sc/snRNA reference (snRNA for TNBC_01/TNBC_02,
#                     scRNA otherwise -- see `sc_modality` in config.R)
# Outputs:  <processed_data_dir>/<platform>/<sample_id>/so.rds   same object, updated in
#                     place with the SingleR_labels202605 metadata column added

library(Seurat)
library(SingleR)
library(SingleCellExperiment)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

sample_id <- "TNBC_02"   # any ID in `xenium_samples` / `merscope_samples`; see config.R
platform  <- "Xenium"    # or "MERSCOPE"

so_path  <- file.path(processed_data_dir, platform, sample_id, "so.rds")
ref_file <- file.path(processed_data_dir, "scRNA", sample_id, paste0("seu_", sample_id, ".rds"))

## SingleR annotation
# Reference cells labelled "mixture" (unresolved sub-clustering, see single_cell_reference/
# 03_scrna_annotate_cell_types.R) are excluded from the reference before label transfer.

ref.seu <- readRDS(ref_file)
ref.sce <- as.SingleCellExperiment(ref.seu)
ref.sce <- ref.sce[, ref.sce$cell_type202605 != "mixture"]

so <- readRDS(so_path)
sce <- as.SingleCellExperiment(so)
pred <- SingleR::SingleR(
  test      = sce,
  ref       = ref.sce,
  labels    = ref.sce$cell_type202605,
  de.method = "wilcox"
)
so$SingleR_labels202605 <- pred$labels

## Save

saveRDS(so, so_path)
