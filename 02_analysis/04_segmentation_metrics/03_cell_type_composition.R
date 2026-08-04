# Purpose:  Cell type composition — proportion of cells per annotated type — for each
#           (sample, platform, segmentation-method) dataset, taken from the
#           SingleR_labels202605 column. Shows whether the choice of segmentation
#           changes the apparent cell type makeup of a sample. See README.md.
# Inputs:   <processed_data_dir>/<platform>/<sample>/so.rds          default segmentation,
#                     annotated by preprocessing/iST/02_annotate_iST.R
#           <processed_data_dir>/segmentation/<method>/<platform>/<sample>_so.rds
#                     Baysor and Proseg, annotated by segmentation/03_postprocess_seu.R
# Outputs:  <results_table_dir>/props_all.rds
#                     feeds figures/Fig4_segmentation.R panel f (TNBC_01 only)

library(Seurat)
library(dplyr)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

samples <- bench_samples

# Default: annotated so.rds from preprocessing/iST/02_annotate_iST.R
default_dir <- processed_data_dir  # <platform>/<sample_id>/so.rds

# Baysor / Proseg: annotated so.rds from segmentation/, keyed by sample ID
# (e.g. ER_02)
seg_meta <- tibble::tribble(
  ~sample_id,  ~sample_file, ~platform,
  "TNBC_01",    "TNBC_01",     "Xenium",
  "TNBC_02",    "TNBC_02",     "Xenium",
  "ER_02",   "ER_02",     "Xenium",
  "ER_01",   "ER_01",     "Xenium",
  "TNBC_03",   "TNBC_03",     "Xenium",
  "TNBC_04",   "TNBC_04",     "Xenium",
  "TNBC_01",    "TNBC_01",     "MERSCOPE",
  "TNBC_02",    "TNBC_02",     "MERSCOPE",
  "ER_02",   "ER_02",     "MERSCOPE",
  "ER_01",   "ER_01",     "MERSCOPE",
  "TNBC_03",   "TNBC_03",     "MERSCOPE"
)

baysor_xenium_dir    <- file.path(processed_data_dir, "segmentation", "baysor", "Xenium")      # <sample_file>_so.rds
baysor_merscope_dir  <- file.path(processed_data_dir, "segmentation", "baysor", "MERSCOPE")    # <sample_file>_so.rds
proseg_xenium_dir     <- file.path(processed_data_dir, "segmentation", "proseg_v3", "Xenium")    # <sample_file>_so.rds
proseg_merscope_dir   <- file.path(processed_data_dir, "segmentation", "proseg_v3", "MERSCOPE")  # <sample_file>_so.rds

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
  "TNBC_01",  "Xenium",   "TNBC_02",  "Xenium",   "ER_02", "Xenium",
  "ER_01", "Xenium",   "TNBC_03", "Xenium",   "TNBC_04", "Xenium",
  "TNBC_01",  "MERSCOPE", "TNBC_02",  "MERSCOPE", "ER_02", "MERSCOPE",
  "ER_01", "MERSCOPE", "TNBC_03", "MERSCOPE"
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

saveRDS(props_all, file.path(results_table_dir, "props_all.rds"))
