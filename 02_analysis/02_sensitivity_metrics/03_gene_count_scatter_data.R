# Purpose:  Total transcript count per gene, summed across all cells/spots, per sample
#           and platform, on the cross-platform overlapped region. See README.md.
# Inputs:   <processed_data_dir>/overlapped_data_list_AllSample.rds   from ../../integration
#           <processed_data_dir>/shared_genes_withVisium.rds          from ../../integration
# Outputs:  <results_table_dir>/gene_counts_per_sample.rds
#                     feeds the per-gene scatter panels in
#                     figures/Fig2_sensitivity.R (panel h) and
#                     figures/Ext_Fig5_gene_scatter.R

library(Seurat)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

bench_samples_with_MERSCOPE <- bench_samples_merscope

overlapped_data_file <- file.path(processed_data_dir, "overlapped_data_list_AllSample.rds")
shared_genes_file     <- file.path(processed_data_dir, "shared_genes_withVisium.rds")

determine_platforms <- function(x) {
  if (x %in% bench_samples_with_MERSCOPE) c("Xenium", "MERSCOPE") else "Xenium"
}

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)
shared_genes <- readRDS(shared_genes_file)

# per-gene total counts (Xenium, MERSCOPE, Visium raw + area-scaled) for one
# sample, restricted to the shared gene panel actually present in its data
gene_counts_sample <- function(sample, data_used = overlapped_data_list_AllSample) {
  platforms <- determine_platforms(sample)
  xenium <- data_used[[sample]]$Xenium
  visium <- data_used[[sample]]$Visium
  genes_in_data <- intersect(rownames(xenium), rownames(visium))
  sg <- intersect(shared_genes, genes_in_data)
  xenium_counts <- as.matrix(GetAssayData(xenium, layer = "counts", assay = "RNA"))
  visium_counts <- as.matrix(GetAssayData(visium, layer = "counts", assay = "Spatial"))
  df <- data.frame(
    Xenium   = sapply(sg, function(g) sum(xenium_counts[g, ])),
    MERSCOPE = NA_real_,
    Visium   = sapply(sg, function(g) sum(visium_counts[g, ])),
    genes    = sg,
    row.names = sg
  )
  if ("MERSCOPE" %in% platforms) {
    merscope <- data_used[[sample]]$MERSCOPE
    sg2 <- intersect(sg, rownames(merscope))
    merscope_counts <- as.matrix(GetAssayData(merscope, layer = "counts", assay = "RNA"))
    df[sg2, "MERSCOPE"] <- sapply(sg2, function(g) sum(merscope_counts[g, ]))
  }
  # Visium spot counts area-rescaled to hexbin size (402 / 181 = hexbin:spot area ratio)
  df$Visium_scaled <- df$Visium * 402 / 181
  df$sample <- sample
  df
}

gene_counts_per_sample <- do.call(rbind, lapply(bench_samples, gene_counts_sample))
rownames(gene_counts_per_sample) <- NULL

saveRDS(gene_counts_per_sample, file.path(results_table_dir, "gene_counts_per_sample.rds"))
