# Purpose:  Mutually Exclusive Co-expression Rate (MECR): for each pair of marker genes
#           belonging to different cell types, the fraction of expressing cells in which
#           both genes are co-detected. High MECR means transcripts are leaking between
#           neighbouring cells, i.e. poor segmentation. Computed per QC-filtered
#           (sample, platform, segmentation) dataset, with matched sc/snRNA-seq as the
#           reference baseline. See README.md.
# Inputs:   <processed_data_dir>/<platform>/<sample>/so.rds          default segmentation
#           <processed_data_dir>/segmentation/<method>/<platform>/<sample>_so.rds
#                     Baysor and Proseg segmentations
#           <processed_data_dir>/scRNA/<sample>/seu_<sample>.rds     reference baseline
#           <processed_data_dir>/shared_genes_withVisium.rds         from ../../integration
# Outputs:  <results_table_dir>/mecr_summary.rds
#                     feeds figures/Fig4_segmentation.R panel c (after-QC median MECR)
#
# Only the after-QC computation is retained here (source also computes a
# before-QC version, used elsewhere in the source Rmd but not by Fig4c).

library(Seurat)
library(tidyverse)
library(Matrix)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

samples_xenium   <- bench_samples_xenium
samples_merscope <- bench_samples_merscope
seg_levels   <- c("Default", "Baysor", "Proseg", "scRNA")
plt_levels   <- c("Xenium", "MERSCOPE", "scRNA")
group_levels <- c("scRNA", "MERSCOPE_Proseg", "MERSCOPE_Baysor", "MERSCOPE_Default",
                   "Xenium_Proseg", "Xenium_Baysor", "Xenium_Default")

shared_genes_file <- file.path(processed_data_dir, "shared_genes_withVisium.rds")

xenium_after_paths   <- setNames(file.path(file.path(processed_data_dir, "Xenium"), samples_xenium, "so.rds"), samples_xenium)
merscope_after_paths <- setNames(file.path(file.path(processed_data_dir, "MERSCOPE"), samples_merscope, "so.rds"), samples_merscope)
baysor_xenium_dir     <- file.path(processed_data_dir, "segmentation", "baysor", "Xenium")      # <sample>_so.rds
baysor_merscope_dir   <- file.path(processed_data_dir, "segmentation", "baysor", "MERSCOPE")    # <sample>_so.rds
proseg_xenium_dir      <- file.path(processed_data_dir, "segmentation", "proseg_v3", "Xenium")    # <sample>_so.rds
proseg_merscope_dir    <- file.path(processed_data_dir, "segmentation", "proseg_v3", "MERSCOPE")  # <sample>_so.rds
sc_dir <- file.path(processed_data_dir, "scRNA")  # <sample>/seu_<sample>.rds (snRNA for TNBC_01/TNBC_02, scRNA otherwise)

get_sc_seu <- function(sample) {
  readRDS(file.path(sc_dir, sample, paste0("seu_", sample, ".rds")))
}

# ---- marker gene panel ------------------------------------------------------

shared_genes <- readRDS(shared_genes_file)

marker_df <- data.frame(
  gene = c("EPCAM", "KRT19", "KRT8", "CD3E", "CD3D", "CD8A", "NKG7", "MS4A1", "CD79A",
           "PECAM1", "CLDN5", "VWF", "C1QA", "C1QB", "CD14", "FCGR3A", "ITGAX", "ITGAM",
           "PDGFRA", "DPT", "COL1A1", "MYH11", "ACTG2"),
  cell_type = c("Epithelial", "Epithelial", "Epithelial", "T", "T", "T", "T", "B", "B",
                "Endo", "Endo", "Endo", "Macro", "Macro", "Macro", "Macro", "Macro", "Macro",
                "Fibro", "Fibro", "Fibro", "Muscle", "Muscle"),
  stringsAsFactors = FALSE
)
rownames(marker_df) <- marker_df$gene
marker_df <- marker_df[intersect(shared_genes, marker_df$gene), ]

# ---- MECR computation -------------------------------------------------------

get.coexpression.rate <- function(obj, marker_df) {
  coexp.rates <- c()
  gene_pairs  <- c()
  genes <- intersect(rownames(obj), rownames(marker_df))
  if (length(genes) > 25) genes <- sample(genes, 25)
  mtx <- as.matrix(GetAssayData(obj, slot = "counts")[genes, ])
  for (g1 in genes) {
    for (g2 in genes) {
      if ((g1 != g2) && (g1 > g2) && (marker_df[g1, "cell_type"] != marker_df[g2, "cell_type"])) {
        c1 <- mtx[g1, ]; c2 <- mtx[g2, ]
        coexp.rates <- c(coexp.rates, sum(c1 > 0 & c2 > 0) / sum(c1 > 0 | c2 > 0))
        gene_pairs  <- c(gene_pairs, paste(g1, g2, sep = "_"))
      }
    }
  }
  data.frame(gene_pair = gene_pairs, rate = coexp.rates, stringsAsFactors = FALSE)
}

compute_mecr <- function(obj, sample, platform, segmentation, marker_df) {
  message(paste(sample, platform, segmentation, sep = " | "))
  mecr <- get.coexpression.rate(obj, marker_df)
  mecr$sample <- sample; mecr$platform <- platform; mecr$segmentation <- segmentation
  mecr %>% dplyr::select(sample, platform, segmentation, gene_pair, rate)
}

# Baysor objects carry a multi-layer RNA assay; flatten the "Gene Expression"
# layer into a plain single-layer object so GetAssayData(slot="counts") works
load_seurat <- function(path) {
  obj <- readRDS(path)
  if ("RNA" %in% names(obj@assays)) {
    layers <- tryCatch(Layers(obj[["RNA"]]), error = function(e) NULL)
    if (!is.null(layers) && "counts.Gene Expression" %in% layers) {
      obj <- CreateSeuratObject(counts = GetAssayData(obj, assay = "RNA", layer = "counts.Gene Expression"))
    }
  }
  obj
}

prep_mecr <- function(mecr_df) {
  mecr_df %>%
    mutate(
      segmentation = factor(segmentation, levels = seg_levels),
      platform     = factor(platform,     levels = plt_levels),
      group = case_when(platform == "scRNA" ~ "scRNA", TRUE ~ paste(platform, segmentation, sep = "_"))
    ) %>%
    mutate(
      segmentation = as.character(segmentation),
      segmentation = ifelse(platform == "scRNA", "scRNA", segmentation),
      segmentation = factor(segmentation, levels = seg_levels)
    ) %>%
    mutate(group = factor(group, levels = group_levels))
}

# ---- After-QC MECR, all (sample, platform, segmentation) combinations -----

res_after <- list()
for (s in samples_xenium)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(xenium_after_paths[s]), s, "Xenium", "Default", marker_df)
for (s in samples_merscope)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(merscope_after_paths[s]), s, "MERSCOPE", "Default", marker_df)
for (s in samples_xenium)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(file.path(baysor_xenium_dir, paste0(s, "_so.rds"))), s, "Xenium", "Baysor", marker_df)
for (s in samples_merscope)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(file.path(baysor_merscope_dir, paste0(s, "_so.rds"))), s, "MERSCOPE", "Baysor", marker_df)
for (s in samples_xenium)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(file.path(proseg_xenium_dir, paste0(s, "_so.rds"))), s, "Xenium", "Proseg", marker_df)
for (s in samples_merscope)
  res_after[[length(res_after) + 1]] <- compute_mecr(load_seurat(file.path(proseg_merscope_dir, paste0(s, "_so.rds"))), s, "MERSCOPE", "Proseg", marker_df)

mecr_df_after <- do.call(rbind, res_after)
rownames(mecr_df_after) <- NULL

# scRNA reference (already QC-processed; one dataset per sample, no per-platform/segmentation split)
res_scrna <- lapply(samples_xenium, function(s)
  compute_mecr(get_sc_seu(s), s, "scRNA", "scRNA", marker_df))
mecr_df_scrna <- do.call(rbind, res_scrna)
rownames(mecr_df_scrna) <- NULL

mecr_after <- bind_rows(prep_mecr(mecr_df_after), prep_mecr(mecr_df_scrna)) %>%
  mutate(qc_status = "After QC")

# ---- aggregate across gene pairs: one row per (sample, platform, segmentation) ----
# median/mean MECR across all mutually-exclusive marker gene pairs (all pairs
# included, including biologically ambiguous ones -- see MECR_compute.Rmd)

mecr_summary <- mecr_after %>%
  group_by(sample, platform, segmentation, qc_status) %>%
  summarise(median_mecr = median(rate, na.rm = TRUE), mean_mecr = mean(rate, na.rm = TRUE), .groups = "drop") %>%
  mutate(segmentation = factor(segmentation, levels = seg_levels), platform = factor(platform, levels = plt_levels))

saveRDS(mecr_summary, file.path(results_table_dir, "mecr_summary.rds"))
