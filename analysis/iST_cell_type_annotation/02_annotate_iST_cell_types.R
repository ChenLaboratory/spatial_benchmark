# Cell type annotation of iST (Xenium/MERSCOPE) data by SingleR label transfer from the
# matched sc/snRNA-seq reference (`cell_type202605`, from single_cell_reference/03_scrna_annotate_cell_types.R).
# Applied per sample and platform; see README.md for the sample/platform/reference-type table.

library(Seurat)
library(SingleR)
library(SingleCellExperiment)

sample_id <- "sample_id"
platform  <- "Xenium"  # or "MERSCOPE"
ref_type  <- "scRNA"   # or "snRNA" -- see README.md

so_path  <- file.path(base_dir, platform, sample_id, "so.rds")
ref_file <- file.path(base_dir, ref_type, sample_id, paste0("seu_", sample_id, ".rds"))

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
