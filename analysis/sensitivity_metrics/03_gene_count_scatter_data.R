# Total transcript counts per gene, summed across all cells/spots, per sample
# and platform, on the cross-platform overlapped region -- feeds the gene
# count scatter panels in figures/Fig2 (f2h) and figures/figS5_gene_scatter.
# See README.md.

library(Seurat)

overlapped_data_file <- "/path/to/analysis/process_aligned_data/overlapped_data_list_AllSample.rds"
shared_genes_file     <- "/path/to/analysis/process_aligned_data/shared_genes_withVisium.rds"

bench_samples <- c("MH0007", "MH0026", "ER_0360", "ER_0114", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")

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

saveRDS(gene_counts_per_sample, "gene_counts_per_sample.rds")
