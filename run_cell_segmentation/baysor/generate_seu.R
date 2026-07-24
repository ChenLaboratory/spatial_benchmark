# Baysor segmentation.csv (via MTX, see 02_run_baysor_to_seu.sh) -> Seurat object.
# Shared across platforms; the only real differences are cell-ID matching and
# whether density needs to be derived.
# Usage: Rscript generate_seu.R <xenium|merscope> <conf_cutoff>

library(Seurat)
library(dplyr)

args        <- commandArgs(trailingOnly = TRUE)
platform    <- args[1]   # "xenium" or "merscope"
conf_cutoff <- args[2]   # e.g. 0 or 0.9 — used only in output filenames

PLATFORM_DIR <- file.path("/path/to/analysis/cell_segmentation/baysor", platform)
MTX_DIR    <- file.path(PLATFORM_DIR, "results/mtx")
SEU_DIR    <- file.path(PLATFORM_DIR, "results/seu")
BAYSOR_DIR <- PLATFORM_DIR

samples <- if (platform == "xenium") {
  c("ER0114", "ER0360", "MH0007", "MH0026", "TN0177", "TN0554")
} else {
  c("ER0114", "ER0360", "MH0007", "MH0026", "TN0177")  # no TN0554 for MERSCOPE
}

for (sample in samples) {
  message("Processing ", sample, "...")

  mtx_dir    <- file.path(MTX_DIR, paste0(sample, "_conf", conf_cutoff))
  stats_path <- file.path(BAYSOR_DIR, paste0("output_", platform, "_", sample), "segmentation_cell_stats.csv")

  if (!dir.exists(mtx_dir))     { warning("MTX not found, skipping: ", sample); next }
  if (!file.exists(stats_path)) { warning("cell_stats not found, skipping: ", sample); next }

  # counts
  expression_matrix <- Read10X(data.dir = mtx_dir)

  seu <- CreateSeuratObject(
    counts       = expression_matrix,
    assay        = "RNA",
    min.cells    = 0,
    min.features = 0
  )

  # cell stats
  cell_stats <- read.csv(stats_path)

  if (platform == "xenium") {
    # Xenium cell IDs are integers; barcodes.tsv was written with 'cell_' prefix
    cell_stats$cell_id_formatted <- paste0("cell_", cell_stats$cell)
    cell_stats <- cell_stats %>% filter(cell_id_formatted %in% colnames(seu))
    rownames(cell_stats) <- cell_stats$cell_id_formatted
  } else {
    # MERSCOPE cell IDs are alphanumeric strings (e.g. CR6d88aca73-1)
    # barcodes.tsv was written without any prefix, so colnames(seu) match cell_stats$cell directly
    cell_stats <- cell_stats %>% filter(cell %in% colnames(seu))
    rownames(cell_stats) <- cell_stats$cell
  }

  seu <- AddMetaData(seu, metadata = cell_stats)

  if (platform == "xenium") {
    # density not pre-computed for Xenium — derive from counts and area
    seu$density <- seu$nCount_RNA / seu$area
  }
  # density is already in segmentation_cell_stats.csv for MERSCOPE — no need to recompute

  # spatial coordinates
  sp_embeddings <- seu@meta.data[, c("x", "y")]
  colnames(sp_embeddings) <- c("Sp_1", "Sp_2")
  seu[["spatial"]] <- CreateDimReducObject(
    embeddings = as.matrix(sp_embeddings),
    assay      = "RNA",
    key        = "Sp_"
  )

  out_path <- file.path(SEU_DIR, paste0(sample, "_so_raw.rds"))
  saveRDS(seu, out_path)
  message("Saved ", out_path)
}
