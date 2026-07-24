# After-alignment, bin-level gene and transcript counts, shared gene panel
# only, unscaled Visium values -- computed on the cross-platform overlapped
# region produced by ../process_aligned_data. See README.md.

library(Seurat)
library(dplyr)
library(reshape2)

overlapped_data_file <- "/path/to/analysis/process_aligned_data/overlapped_data_list_AllSample.rds"
shared_genes_file     <- "/path/to/analysis/process_aligned_data/shared_genes_withVisium.rds"

bench_samples <- c("MH0007", "MH0026", "ER_0360", "ER_0114", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")

determine_platforms <- function(x) {
  if (x %in% bench_samples_with_MERSCOPE) c("Xenium", "MERSCOPE") else "Xenium"
}

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)
shared_genes <- readRDS(shared_genes_file)

## ---- detected genes per hexbin, shared gene panel -------------------------
# sum iST cell counts within each Visium hexbin, then count genes with >0 reads
get_unique_gene_counts <- function(seurat_obj, platform_name, sample) {
  meta <- seurat_obj@meta.data
  counts <- as.matrix(GetAssayData(seurat_obj, layer = "counts"))

  aggregated_counts <- aggregate(t(counts), by = list(meta$Visium_spot_id), FUN = sum)
  rownames(aggregated_counts) <- aggregated_counts$Group.1
  aggregated_counts <- aggregated_counts[-1, ]
  aggregated_counts$Group.1 <- NULL
  aggregated_counts <- t(aggregated_counts)
  aggregated_counts <- aggregated_counts[rownames(aggregated_counts) %in% shared_genes, ]

  unique_genes <- colSums(aggregated_counts > 0)
  data.frame(Visium_spot_id = colnames(aggregated_counts), Unique_Genes = unique_genes,
             Platform = platform_name, sample = sample)
}

bin_gene_counts_shared_genes <- data.frame()

for (sample in names(overlapped_data_list_AllSample)) {
  data_used <- overlapped_data_list_AllSample[[sample]]

  visium_data <- data_used$Visium
  counts <- as.matrix(GetAssayData(visium_data, layer = "counts"))
  counts <- counts[rownames(counts) %in% shared_genes, ]
  unique_genes <- colSums(counts > 0)
  bin_gene_counts_shared_genes <- rbind(bin_gene_counts_shared_genes, data.frame(
    Visium_spot_id = colnames(counts), Unique_Genes = unique_genes, Platform = "Visium", sample = sample
  ))

  bin_gene_counts_shared_genes <- rbind(bin_gene_counts_shared_genes,
    get_unique_gene_counts(data_used$Xenium, "Xenium", sample))

  if (sample %in% bench_samples_with_MERSCOPE) {
    bin_gene_counts_shared_genes <- rbind(bin_gene_counts_shared_genes,
      get_unique_gene_counts(data_used$MERSCOPE, "MERSCOPE", sample))
  }
}
bin_gene_counts_shared_genes$Visium_spot_id <- NULL

saveRDS(bin_gene_counts_shared_genes, "bin_gene_counts_shared_genes_after_alignment.rds")

## ---- total transcript counts per hexbin, shared gene panel, unscaled ------
## Visium (each Visium hexbin's own shared-panel spot count, compared against
## the summed iST shared-panel counts of cells assigned to that hexbin -- no
## area-based rescaling of the Visium value)
bin_nCount_shared_genes_unscaled <- data.frame()

for (sample in bench_samples) {
  meta <- overlapped_data_list_AllSample[[sample]]$Visium@meta.data
  platforms <- determine_platforms(sample)
  vars <- c("nCount_Spatial_shared_genes", paste(platforms, "total_nCount_hexbin_shared_genes", sep = "_"))

  data_melted <- reshape2::melt(meta, id.vars = "orig.ident", measure.vars = vars) %>%
    mutate(sample = sample, platform = sapply(strsplit(as.character(variable), "_"), `[`, 1))
  bin_nCount_shared_genes_unscaled <- rbind(bin_nCount_shared_genes_unscaled, data_melted)
}
bin_nCount_shared_genes_unscaled$platform <- gsub("nCount", "Visium", bin_nCount_shared_genes_unscaled$platform)

saveRDS(bin_nCount_shared_genes_unscaled, "bin_nCount_shared_genes_unscaled_after_alignment.rds")
