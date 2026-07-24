# Global Moran's I for every non-control gene (bin-level, same method as the
# negative controls in 02_individual_control_moran.R), plus total transcript
# counts per gene/probe/codeword per sample and platform. Feeds Fig3 panel k
# (Moran's I vs. total signal, gene-targeting probes vs. negative controls).
# See README.md.

library(Seurat)
library(SpatialExperiment)
library(scider)
library(dplyr)
library(tibble)
library(arrow)

overlapped_data_file <- "/path/to/analysis/process_aligned_data/overlapped_data_list_AllSample.rds"
xenium_raw_transcripts_files <- c(
  "MH0007"  = "/path/to/raw/Xenium/MH0007/transcripts.parquet",
  "MH0026"  = "/path/to/raw/Xenium/MH0026/transcripts.parquet",
  "ER_0360" = "/path/to/raw/Xenium/ER_0360/transcripts.parquet",
  "ER_0114" = "/path/to/raw/Xenium/ER_0114/transcripts.parquet",
  "TN_0177" = "/path/to/raw/Xenium/TN_0177/transcripts.parquet",
  "TN_0554" = "/path/to/raw/Xenium/TN_0554/transcripts.parquet"
)

bench_samples <- c("MH0007", "MH0026", "ER_0360", "ER_0114", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")
bench_samples_with_MERSCOPE_v1 <- c("MH0026", "MH0007")

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)

visium_seu2spe_withImg <- function(seu, sample_id = "sample01", img_id) {
  sce <- Seurat::as.SingleCellExperiment(seu)
  spatialCoords <- as.matrix(GetTissueCoordinates(seu, scale = NULL)[, c(2, 1)])
  img <- SpatialExperiment::SpatialImage(x = as.raster(GetImage(seu, mode = "raster")))
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = names(Seurat::ScaleFactors(seu[[img_id]])),
    data = I(list(img)),
    scaleFactor = as.numeric(Seurat::ScaleFactors(seu[[img_id]]))
  )
  if ("lowres" %in% imgData$image_id) {
    lowres_row <- imgData[imgData$image_id == "lowres", , drop = FALSE]
    other_rows <- imgData[imgData$image_id != "lowres", , drop = FALSE]
    imgData <- rbind(lowres_row, other_rows)
  }
  spe <- SpatialExperiment(
    assays = assays(sce), rowData = rowData(sce), colData = colData(sce),
    metadata = metadata(sce), reducedDims = reducedDims(sce), altExps = altExps(sce),
    sample_id = sample_id, spatialCoords = spatialCoords, imgData = imgData
  )
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  spe
}

is_control_gene_xenium <- function(gene) grepl("NegControl", gene)
is_control_gene_merscope <- function(gene) grepl("Blank", gene) | grepl("Sabrina", gene) | grepl("^FC[0-9]+$", gene)

compute_moran <- function(spe, feature_name) {
  dat <- spe[[feature_name]]
  res <- globalMoran(spe, data1 = dat, at = "cell", permutations = 9999)
  tibble(feature = feature_name, lisa = res$lisa, pval = res$p.val)
}

## ---- Global Moran's I for every non-control gene --------------------------

all_results <- list()

for (sample in bench_samples) {
  message("==== ", sample, " ====")
  visium_seu <- overlapped_data_list_AllSample[[sample]]$Visium

  xenium_seu <- overlapped_data_list_AllSample[[sample]]$Xenium
  xenium_genes <- rownames(xenium_seu)[!is_control_gene_xenium(rownames(xenium_seu))]
  counts_x <- as.data.frame(t(as.matrix(GetAssayData(xenium_seu, layer = "counts"))))
  counts_x[, "Visium_spot_id"] <- xenium_seu$Visium_spot_id
  gene_hexbin_x <- counts_x %>%
    group_by(Visium_spot_id) %>%
    summarise(across(all_of(xenium_genes), sum, .names = "Xenium_{col}"), .groups = "drop") %>%
    filter(Visium_spot_id != "")

  visium_seu_x <- visium_seu
  visium_seu_x@meta.data <- visium_seu_x@meta.data %>%
    rownames_to_column(var = "Visium_spot_id") %>%
    left_join(gene_hexbin_x, by = "Visium_spot_id") %>%
    column_to_rownames("Visium_spot_id")

  spe_x <- visium_seu2spe_withImg(visium_seu_x, img_id = "slice1")
  spe_x <- findNbrsSpatial(spe_x, k = 10)

  xenium_feat_cols <- intersect(paste0("Xenium_", xenium_genes), colnames(spe_x@colData))
  for (feat in xenium_feat_cols) {
    tmp <- compute_moran(spe_x, feat)
    tmp$sample <- sample; tmp$platform <- "Xenium"; tmp$gene <- sub("^Xenium_", "", feat)
    all_results[[length(all_results) + 1]] <- tmp
  }

  if (sample %in% bench_samples_with_MERSCOPE) {
    merscope_seu <- overlapped_data_list_AllSample[[sample]]$MERSCOPE
    merscope_genes <- rownames(merscope_seu)[!is_control_gene_merscope(rownames(merscope_seu))]
    counts_m <- as.data.frame(t(as.matrix(GetAssayData(merscope_seu, layer = "counts"))))
    counts_m[, "Visium_spot_id"] <- merscope_seu$Visium_spot_id
    gene_hexbin_m <- counts_m %>%
      group_by(Visium_spot_id) %>%
      summarise(across(all_of(merscope_genes), sum, .names = "MERSCOPE_{col}"), .groups = "drop") %>%
      filter(Visium_spot_id != "")

    visium_seu_m <- visium_seu
    visium_seu_m@meta.data <- visium_seu_m@meta.data %>%
      rownames_to_column(var = "Visium_spot_id") %>%
      left_join(gene_hexbin_m, by = "Visium_spot_id") %>%
      column_to_rownames("Visium_spot_id")

    spe_m <- visium_seu2spe_withImg(visium_seu_m, img_id = "slice1")
    spe_m <- findNbrsSpatial(spe_m, k = 10)

    merscope_feat_cols <- intersect(paste0("MERSCOPE_", merscope_genes), colnames(spe_m@colData))
    for (feat in merscope_feat_cols) {
      tmp <- compute_moran(spe_m, feat)
      tmp$sample <- sample; tmp$platform <- "MERSCOPE"; tmp$gene <- sub("^MERSCOPE_", "", feat)
      all_results[[length(all_results) + 1]] <- tmp
    }
  }
}

all_moran_genes <- bind_rows(all_results)
all_moran_genes$platform2 <- dplyr::case_when(
  grepl("MERSCOPE", all_moran_genes$platform) & (all_moran_genes$sample %in% bench_samples_with_MERSCOPE_v1) ~ "MERSCOPE_V1",
  grepl("MERSCOPE", all_moran_genes$platform) & !(all_moran_genes$sample %in% bench_samples_with_MERSCOPE_v1) ~ "MERSCOPE_V2",
  grepl("Xenium", all_moran_genes$platform) ~ "Xenium",
  TRUE ~ "Unknown"
)
all_moran_genes$sample_factor <- factor(all_moran_genes$sample, levels = bench_samples)
all_moran_genes$platform2 <- factor(all_moran_genes$platform2, levels = c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium"))
all_moran_genes$sig <- all_moran_genes$pval < 0.05

saveRDS(all_moran_genes, "all_moran_genes.rds")

## ---- total transcript counts per gene/probe/codeword ----------------------
## (non-control genes + Xenium NegControls pulled from raw transcripts +
## MERSCOPE controls already in the count matrix)

get_xenium_neg_counts <- function(sample, xenium_seu) {
  df <- read_parquet(xenium_raw_transcripts_files[sample])
  neg_counts <- df %>%
    filter(grepl("NegControl", feature_name), cell_id != "UNASSIGNED") %>%
    group_by(feature_name, cell_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(feature_name, cell_id = colnames(xenium_seu), fill = list(count = 0))
  neg_wide <- neg_counts %>%
    pivot_wider(names_from = cell_id, values_from = count, values_fill = 0) %>%
    column_to_rownames("feature_name")
  neg_wide <- neg_wide[, intersect(colnames(neg_wide), colnames(xenium_seu)), drop = FALSE]
  data.frame(feature_clean = rownames(neg_wide), total_counts = rowSums(neg_wide), stringsAsFactors = FALSE)
}

scatter_counts_list <- list()

for (sample in bench_samples) {
  xenium_seu <- overlapped_data_list_AllSample[[sample]]$Xenium
  counts_x <- GetAssayData(xenium_seu, layer = "counts")

  scatter_counts_list[[paste(sample, "Xenium")]] <- bind_rows(
    data.frame(feature_clean = rownames(counts_x), total_counts = rowSums(counts_x), stringsAsFactors = FALSE),
    get_xenium_neg_counts(sample, xenium_seu)
  ) %>% mutate(sample = sample, platform2 = "Xenium")

  if (sample %in% bench_samples_with_MERSCOPE) {
    counts_m <- GetAssayData(overlapped_data_list_AllSample[[sample]]$MERSCOPE, layer = "counts")
    plat2 <- ifelse(sample %in% bench_samples_with_MERSCOPE_v1, "MERSCOPE_V1", "MERSCOPE_V2")
    scatter_counts_list[[paste(sample, "MERSCOPE")]] <- data.frame(
      feature_clean = rownames(counts_m), total_counts = rowSums(counts_m),
      sample = sample, platform2 = plat2, stringsAsFactors = FALSE
    )
  }
}

total_counts_df <- bind_rows(scatter_counts_list) %>%
  mutate(feature_clean = make.names(feature_clean))

saveRDS(total_counts_df, "total_counts_df.rds")
