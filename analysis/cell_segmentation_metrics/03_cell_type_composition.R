# Cell type composition (proportion of cells per annotated type) for each
# (sample, platform, segmentation-method) dataset, read from the
# SingleR-annotated `SingleR_labels202605` column. Feeds Fig4f (MH0007 only).
# See README.md.

library(Seurat)
library(dplyr)

samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")

# Default: annotated so.rds from preprocessing/iST_cell_type_annotation
default_dir <- "/path/to/preprocessing"  # <platform>/<sample_id>/so.rds

# Baysor / Proseg: annotated so.rds from run_cell_segmentation, sample IDs
# without underscore (e.g. ER0114)
seg_meta <- tibble::tribble(
  ~sample_id,  ~sample_file, ~platform,
  "MH0007",    "MH0007",     "Xenium",
  "MH0026",    "MH0026",     "Xenium",
  "ER_0114",   "ER0114",     "Xenium",
  "ER_0360",   "ER0360",     "Xenium",
  "TN_0177",   "TN0177",     "Xenium",
  "TN_0554",   "TN0554",     "Xenium",
  "MH0007",    "MH0007",     "MERSCOPE",
  "MH0026",    "MH0026",     "MERSCOPE",
  "ER_0114",   "ER0114",     "MERSCOPE",
  "ER_0360",   "ER0360",     "MERSCOPE",
  "TN_0177",   "TN0177",     "MERSCOPE"
)

baysor_xenium_dir    <- "/path/to/run_cell_segmentation/baysor/results/xenium"      # <sample_file>_so.rds
baysor_merscope_dir  <- "/path/to/run_cell_segmentation/baysor/results/merscope"    # <sample_file>_so.rds
proseg_xenium_dir     <- "/path/to/run_cell_segmentation/proseg_v3/results/xenium"    # <sample_file>_so.rds
proseg_merscope_dir   <- "/path/to/run_cell_segmentation/proseg_v3/results/merscope"  # <sample_file>_so.rds

# proportion of cells per annotated cell type, for one dataset
compute_props <- function(labels, sample_id, platform, method) {
  tab <- table(labels)
  df  <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(df) <- c("cell_type", "count")
  df$proportion <- df$count / sum(df$count)
  df$sample_id <- sample_id; df$platform <- platform; df$method <- method
  df[, c("sample_id", "platform", "method", "cell_type", "count", "proportion")]
}

## ---- Default -----------------------------------------------------------

default_meta <- tibble::tribble(
  ~sample_id,  ~platform,
  "MH0007",  "Xenium",   "MH0026",  "Xenium",   "ER_0114", "Xenium",
  "ER_0360", "Xenium",   "TN_0177", "Xenium",   "TN_0554", "Xenium",
  "MH0007",  "MERSCOPE", "MH0026",  "MERSCOPE", "ER_0114", "MERSCOPE",
  "ER_0360", "MERSCOPE", "TN_0177", "MERSCOPE"
)

props_default <- lapply(seq_len(nrow(default_meta)), function(i) {
  sid <- default_meta$sample_id[i]; plt <- default_meta$platform[i]
  so <- readRDS(file.path(default_dir, plt, sid, "so.rds"))
  compute_props(so$SingleR_labels202605, sid, plt, "default")
})
props_default <- do.call(rbind, props_default)

## ---- Baysor --------------------------------------------------------------

props_baysor <- lapply(seq_len(nrow(seg_meta)), function(i) {
  sid <- seg_meta$sample_id[i]; sfile <- seg_meta$sample_file[i]; plt <- seg_meta$platform[i]
  dir <- if (plt == "Xenium") baysor_xenium_dir else baysor_merscope_dir
  so <- readRDS(file.path(dir, paste0(sfile, "_so.rds")))
  compute_props(so$SingleR_labels202605, sid, plt, "baysor")
})
props_baysor <- do.call(rbind, props_baysor)

## ---- Proseg ---------------------------------------------------------------

props_proseg <- lapply(seq_len(nrow(seg_meta)), function(i) {
  sid <- seg_meta$sample_id[i]; sfile <- seg_meta$sample_file[i]; plt <- seg_meta$platform[i]
  dir <- if (plt == "Xenium") proseg_xenium_dir else proseg_merscope_dir
  so <- readRDS(file.path(dir, paste0(sfile, "_so.rds")))
  compute_props(so$SingleR_labels202605, sid, plt, "proseg_v3")
})
props_proseg <- do.call(rbind, props_proseg)

## ---- combine and save ------------------------------------------------------

props_all <- rbind(props_default, props_baysor, props_proseg)
props_all$method     <- factor(props_all$method, levels = c("default", "baysor", "proseg_v3"))
props_all$platform   <- factor(props_all$platform, levels = c("Xenium", "MERSCOPE"))
props_all$sample_id  <- factor(props_all$sample_id, levels = samples)

saveRDS(props_all, "props_all.rds")
