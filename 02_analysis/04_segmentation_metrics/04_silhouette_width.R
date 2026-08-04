# Purpose:  Silhouette width in PCA space, grouped by the annotated cell type labels
#           (SingleR_labels202605). Higher values mean cells of the same annotated type
#           form tighter, better-separated transcriptional groups — a readout of
#           segmentation quality relative to expected biology. See README.md.
# Inputs:   <processed_data_dir>/<platform>/<sample>/so.rds          default segmentation
#           <processed_data_dir>/segmentation/<method>/<platform>/<sample>_so.rds
#                     Baysor and Proseg segmentations
# Outputs:  <results_table_dir>/sil_ct_results.rds
#                     feeds figures/Fig4_segmentation.R panel g
#
# Rare cell types (<10 cells) are excluded, and each retained type is
# downsampled to at most 10,000 cells before computing silhouette width, to
# keep computation tractable and stop large dominant types (e.g. Tumor) from
# disproportionately driving the mean.

library(Seurat)
library(dplyr)
library(bluster)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

samples <- bench_samples

default_dir <- processed_data_dir  # <platform>/<sample>/so.rds
baysor_xenium_dir     <- file.path(processed_data_dir, "segmentation", "baysor", "Xenium")      # <sample_file>_so.rds
baysor_merscope_dir   <- file.path(processed_data_dir, "segmentation", "baysor", "MERSCOPE")    # <sample_file>_so.rds
proseg_xenium_dir      <- file.path(processed_data_dir, "segmentation", "proseg_v3", "Xenium")    # <sample_file>_so.rds
proseg_merscope_dir    <- file.path(processed_data_dir, "segmentation", "proseg_v3", "MERSCOPE")  # <sample_file>_so.rds

seg_meta <- tibble::tribble(
  ~sample,     ~sample_file, ~platform,
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

sil_manifest <- bind_rows(
  seg_meta %>% mutate(method = "default", path = file.path(default_dir, platform, sample, "so.rds")),
  seg_meta %>% mutate(method = "baysor",
                       path = ifelse(platform == "Xenium",
                                      file.path(baysor_xenium_dir, paste0(sample_file, "_so.rds")),
                                      file.path(baysor_merscope_dir, paste0(sample_file, "_so.rds")))),
  seg_meta %>% mutate(method = "proseg_v3",
                       path = ifelse(platform == "Xenium",
                                      file.path(proseg_xenium_dir, paste0(sample_file, "_so.rds")),
                                      file.path(proseg_merscope_dir, paste0(sample_file, "_so.rds"))))
)

# Mean silhouette width across cells, grouped by annotated cell type in PCA
# space (rare types <10 cells excluded; each retained type downsampled to
# at most `downsample` cells)
getSilhouetteWidth_ct <- function(seu_obj, downsample = 10000) {
  labels <- seu_obj$SingleR_labels202605

  ct_counts  <- table(labels)
  keep_ct    <- names(ct_counts)[ct_counts >= 10]
  seu_sub    <- seu_obj[, labels %in% keep_ct]

  Idents(seu_sub) <- seu_sub$SingleR_labels202605
  seu_sub <- subset(seu_sub, downsample = downsample)

  silhouette <- as.data.frame(bluster::approxSilhouette(
    Embeddings(seu_sub, "pca")[, 1:10], clusters = seu_sub$SingleR_labels202605
  ))

  data.frame(sample_id = unique(seu_sub$sample_id), platform = unique(seu_sub$platform),
             value = round(mean(silhouette$width), digits = 3))
}

compute_silhouette_ct <- function(path, sample_id, platform, method) {
  so <- readRDS(path)
  so$sample_id <- sample_id
  so$platform  <- platform
  res <- getSilhouetteWidth_ct(so)
  res$method <- method
  res
}

sil_ct_results <- do.call(rbind, lapply(seq_len(nrow(sil_manifest)), function(i) {
  m <- sil_manifest[i, ]
  compute_silhouette_ct(m$path, m$sample, m$platform, m$method)
}))

sil_ct_results$method    <- factor(sil_ct_results$method, levels = c("default", "baysor", "proseg_v3"))
sil_ct_results$platform  <- factor(sil_ct_results$platform, levels = c("Xenium", "MERSCOPE"))
sil_ct_results$sample_id <- factor(sil_ct_results$sample_id, levels = samples)

saveRDS(sil_ct_results, file.path(results_table_dir, "sil_ct_results.rds"))
