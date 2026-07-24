# Bin-level transcript and gene detection metrics per sample/platform, computed
# independently per platform *before* cross-platform alignment. iST (Xenium/
# MERSCOPE) cells are binned into 110um hexagons (matching the Visium spot
# centre-to-centre distance); Visium itself is not re-binned, its native spots
# are used directly. Two gene panels: all genes, and the shared panel across
# Xenium/MERSCOPE/Visium. See README.md.

library(Seurat)
library(SpatialExperiment)
library(scider)
library(dplyr)

xenium_dir        <- "/path/to/preprocessing/iST/Xenium"     # <sample>/so_raw.rds
visium_dir        <- "/path/to/preprocessing/Visium"         # <sample>/so_raw.rds
merscope_dir       <- "/path/to/preprocessing/iST/MERSCOPE"   # <sample>/so_raw.rds, v1_so_raw.rds, v2_so_raw.rds
shared_genes_file <- "/path/to/analysis/process_aligned_data/shared_genes_withVisium.rds"

# all_samples includes the 6-sample benchmarking cohort plus one pilot sample
# (24MH0042, not part of the benchmarking cohort but included here to match
# the original per-platform sensitivity computation; figures filter it out)
all_samples <- c("MH0007", "MH0026", "ER_0360", "ER_0114", "TN_0177", "TN_0554", "24MH0042")
bench_samples_with_MERSCOPE_v1 <- c("MH0026", "MH0007")
bench_samples_with_MERSCOPE_v2 <- c("ER_0360", "ER_0114", "TN_0177")

shared_genes <- readRDS(shared_genes_file)

# ---- bin-level helpers -----------------------------------------------------

ST_so2spe <- function(so) {
  counts <- GetAssayData(so, slot = "counts")
  spatial_coords <- so@reductions$spatial@cell.embeddings
  colData <- so@meta.data

  SpatialExperiment(assays = list(counts = counts), colData = colData, spatialCoords = spatial_coords)
}

BinSPE <- function(spe, grid.type = "hex", grid.length.x = 110, ct_group = "orig.ident") {
  ct1 <- spe@colData[[ct_group]]
  coi <- as.character(unique(ct1))
  spe_hex_cell <- scider::gridDensity(spe, coi = coi, grid.type = grid.type, grid.length.x = grid.length.x, id = ct_group)
  spe_hex <- scider::gridSPE(spe_hex_cell, cell.count = TRUE, id = ct_group)
  spe_hex[, spe_hex$cell_counts_overall > 0]
}

process_seu <- function(seu, platform, sample_id) {
  spe <- ST_so2spe(seu)
  spe_hex <- BinSPE(spe)
  data.frame(
    sample = rep(sample_id, ncol(spe_hex)),
    platform = rep(platform, ncol(spe_hex)),
    n_cells = spe_hex$cell_counts_overall,
    n_transcripts = colSums(counts(spe_hex)),
    n_genes = colSums(counts(spe_hex) > 0)
  )
}

# MERSCOPE control probes removed before computing any metric
drop_merscope_controls <- function(seu) {
  seu[!(grepl("Blank", rownames(seu)) | grepl("Sabrina", rownames(seu)) | grepl("^FC[0-9]+$", rownames(seu))), ]
}

# ---- main loop --------------------------------------------------------------
# gene_subset = NULL -> all genes; gene_subset = shared_genes -> shared panel

compute_bin_summary <- function(gene_subset = NULL) {
  df_bins_summary <- data.frame()

  for (s in all_samples) {
    message("=== ", s, " ===")

    xenium_path <- file.path(xenium_dir, s, "so_raw.rds")
    visium_path <- file.path(visium_dir, s, "so_raw.rds")
    merscope_path <- file.path(merscope_dir, s, "so_raw.rds")
    merscope_v1_path <- file.path(merscope_dir, s, "v1_so_raw.rds")
    merscope_v2_path <- file.path(merscope_dir, s, "v2_so_raw.rds")

    if (file.exists(xenium_path)) {
      seu <- readRDS(xenium_path)
      if (!is.null(gene_subset)) seu <- seu[gene_subset, ]
      df_bins_summary <- rbind(df_bins_summary, process_seu(seu, "Xenium", s))
    }

    if (file.exists(visium_path)) {
      seu <- readRDS(visium_path)
      if (!is.null(gene_subset)) seu <- seu[gene_subset, ]
      counts <- as.data.frame(seu@assays$Spatial@layers$counts)
      df_bins_summary <- rbind(df_bins_summary, data.frame(
        sample = rep(s, ncol(seu)),
        platform = rep("Visium", ncol(seu)),
        n_cells = rep(NA, ncol(seu)),
        n_transcripts = colSums(counts),
        n_genes = colSums(counts > 0)
      ))
    }

    # MERSCOPE run used in the main benchmarking cohort
    if (file.exists(merscope_path)) {
      seu <- drop_merscope_controls(readRDS(merscope_path))
      if (!is.null(gene_subset)) seu <- seu[gene_subset, ]
      mer_col <- if (s %in% bench_samples_with_MERSCOPE_v1) "MERSCOPE_V1" else "MERSCOPE_V2"
      df_bins_summary <- rbind(df_bins_summary, process_seu(seu, mer_col, s))
    }

    # extra MERSCOPE v1/v2 acquisitions for samples not covered above
    # (e.g. TN_0554, which only has a v1 run outside the main cohort pipeline)
    if (file.exists(merscope_v1_path)) {
      seu <- drop_merscope_controls(readRDS(merscope_v1_path))
      if (!is.null(gene_subset)) seu <- seu[gene_subset, ]
      df_bins_summary <- rbind(df_bins_summary, process_seu(seu, "MERSCOPE_V1", s))
    }
    if (file.exists(merscope_v2_path)) {
      seu <- drop_merscope_controls(readRDS(merscope_v2_path))
      if (!is.null(gene_subset)) seu <- seu[gene_subset, ]
      df_bins_summary <- rbind(df_bins_summary, process_seu(seu, "MERSCOPE_V2", s))
    }
  }

  df_bins_summary
}

df_bins_summary_all_genes <- compute_bin_summary(gene_subset = NULL)
saveRDS(df_bins_summary_all_genes, "df_bins_summary_110um_all_ST_all_genes_before_alignment.rds")

df_bins_summary_shared_genes <- compute_bin_summary(gene_subset = shared_genes)
saveRDS(df_bins_summary_shared_genes, "df_bins_summary_110um_all_ST_213_shared_genes_before_alignment.rds")
