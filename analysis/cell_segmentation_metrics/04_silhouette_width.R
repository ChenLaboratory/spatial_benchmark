# Silhouette width computed with annotated cell type labels (SingleR_labels202605)
# as the grouping variable in PCA space: higher values indicate cells of the
# same annotated type form tighter, better-separated transcriptional groups --
# a direct readout of segmentation quality relative to expected biology. Feeds
# Fig4g. See README.md.
#
# Rare cell types (<10 cells) are excluded, and each retained type is
# downsampled to at most 10,000 cells before computing silhouette width, to
# keep computation tractable and stop large dominant types (e.g. Tumor) from
# disproportionately driving the mean.

library(Seurat)
library(dplyr)
library(bluster)

samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")

default_dir <- "/path/to/preprocessing"  # <platform>/<sample>/so.rds
baysor_xenium_dir     <- "/path/to/run_cell_segmentation/baysor/results/xenium"      # <sample_file>_so.rds
baysor_merscope_dir   <- "/path/to/run_cell_segmentation/baysor/results/merscope"    # <sample_file>_so.rds
proseg_xenium_dir      <- "/path/to/run_cell_segmentation/proseg_v3/results/xenium"    # <sample_file>_so.rds
proseg_merscope_dir    <- "/path/to/run_cell_segmentation/proseg_v3/results/merscope"  # <sample_file>_so.rds

seg_meta <- tibble::tribble(
  ~sample,     ~sample_file, ~platform,
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

saveRDS(sil_ct_results, "sil_ct_results.rds")
