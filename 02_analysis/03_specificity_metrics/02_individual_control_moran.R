# Purpose:  Global Moran's I (spatial autocorrelation) for each individual negative
#           control probe/codeword, computed on its per-hexbin count, testing whether
#           control signal is spatially random. See README.md.
# Inputs:   <processed_data_dir>/overlapped_data_list_AllSample.rds   from ../../integration
#           raw_xenium_dir[<sample>]/transcripts.parquet   (set in config.R)
#                     Xenium negative controls are absent from the Seurat object's
#                     expression matrix, so their per-cell counts come from here
# Outputs:  <results_table_dir>/all_moran_individual_NC.rds
#                     feeds figures/Fig3_specificity.R panel i, and panel k together
#                     with 04_individual_gene_moran.R

library(Seurat)
library(SpatialExperiment)
library(scider)
library(dplyr)
library(tibble)
library(arrow)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

bench_samples_with_MERSCOPE_v1 <- bench_samples_merscope_v1
bench_samples_with_MERSCOPE_v2 <- bench_samples_merscope_v2

overlapped_data_file <- file.path(processed_data_dir, "overlapped_data_list_AllSample.rds")
# raw Xenium transcripts.parquet per sample (negative control rows are
# excluded from the Xenium Seurat object, so pulled from here instead)
xenium_raw_transcripts_files <- setNames(
  file.path(raw_xenium_dir[bench_samples], "transcripts.parquet"), bench_samples)

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)

# convert a Visium Seurat object (with spatial image) to a SpatialExperiment;
# custom wrapper (not part of any CRAN/Bioconductor package)
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

get_xenium_neg_matrix <- function(sample) {
  df <- read_parquet(xenium_raw_transcripts_files[sample])
  all_cells <- df %>% filter(cell_id != "UNASSIGNED") %>% pull(cell_id) %>% unique()
  neg_df <- df %>% filter(grepl("NegControl", feature_name), cell_id != "UNASSIGNED")
  neg_counts <- neg_df %>%
    group_by(feature_name, cell_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(feature_name, cell_id = all_cells, fill = list(count = 0))
  neg_matrix <- neg_counts %>%
    pivot_wider(names_from = cell_id, values_from = count, values_fill = 0) %>%
    as.data.frame()
  rownames(neg_matrix) <- neg_matrix$feature_name
  neg_matrix$feature_name <- NULL
  neg_matrix
}

## ---- per-hexbin counts for each individual NC probe/codeword -------------

for (sample in bench_samples) {
  visium_seu <- overlapped_data_list_AllSample[[sample]]$Visium

  if (!is.null(overlapped_data_list_AllSample[[sample]]$MERSCOPE)) {
    iST_seu <- overlapped_data_list_AllSample[[sample]]$MERSCOPE
    counts <- as.data.frame(t(as.matrix(GetAssayData(iST_seu, layer = "counts"))))
    counts[, "Visium_spot_id"] <- iST_seu$Visium_spot_id
    NC <- rownames(iST_seu)[grepl("Sabrina", rownames(iST_seu)) |
                             grepl("^FC[0-9]+$", rownames(iST_seu)) |
                             grepl("Blank", rownames(iST_seu))]
    NC_summed <- counts %>%
      group_by(Visium_spot_id) %>%
      summarise(across(all_of(NC), sum, .names = "MERSCOPE_{col}")) %>%
      filter(Visium_spot_id != "")
    visium_seu@meta.data <- visium_seu@meta.data %>%
      rownames_to_column(var = "Visium_spot_id") %>%
      left_join(NC_summed, by = "Visium_spot_id") %>%
      column_to_rownames("Visium_spot_id")
  }

  iST_seu <- overlapped_data_list_AllSample[[sample]]$Xenium
  neg_matrix <- get_xenium_neg_matrix(sample)
  neg_matrix <- as.data.frame(t(neg_matrix[, colnames(iST_seu)]))
  NC <- colnames(neg_matrix)
  neg_matrix[, "Visium_spot_id"] <- iST_seu$Visium_spot_id
  NC_summed <- neg_matrix %>%
    group_by(Visium_spot_id) %>%
    summarise(across(all_of(NC), sum, .names = "Xenium_{col}")) %>%
    filter(Visium_spot_id != "")
  visium_seu@meta.data <- visium_seu@meta.data %>%
    rownames_to_column(var = "Visium_spot_id") %>%
    left_join(NC_summed, by = "Visium_spot_id") %>%
    column_to_rownames("Visium_spot_id")

  overlapped_data_list_AllSample[[sample]]$Visium <- visium_seu
}

## ---- Global Moran's I per individual NC probe/codeword -------------------

merscope_codeword_pattern <- "^MERSCOPE_Blank"
merscope_probe_pattern    <- "^(MERSCOPE_Sabrina|MERSCOPE_FC)"
xenium_codeword_pattern   <- "^Xenium_NegControlCodeword"
xenium_probe_pattern      <- "^Xenium_NegControlProbe"

compute_moran <- function(spe, feature_name) {
  dat <- spe[[feature_name]]
  res <- globalMoran(spe, data1 = dat, at = "cell", permutations = 9999)
  tibble(feature = feature_name, lisa = res$lisa, pval = res$p.val)
}

all_results <- list()

for (sample in names(overlapped_data_list_AllSample)) {
  visium_seu <- overlapped_data_list_AllSample[[sample]]$Visium
  spe <- visium_seu2spe_withImg(visium_seu, img_id = "slice1")
  spe <- findNbrsSpatial(spe, k = 10)

  feature_names <- colnames(spe@colData)
  feature_sets <- list(
    MERSCOPE_codeword = grep(merscope_codeword_pattern, feature_names, value = TRUE),
    MERSCOPE_probe    = grep(merscope_probe_pattern, feature_names, value = TRUE),
    Xenium_codeword   = grep(xenium_codeword_pattern, feature_names, value = TRUE),
    Xenium_probe      = grep(xenium_probe_pattern, feature_names, value = TRUE)
  )

  for (type in names(feature_sets)) {
    message("Processing: ", sample, " - ", type)
    for (feat in feature_sets[[type]]) {
      tmp <- compute_moran(spe, feat)
      tmp$sample <- sample
      tmp$platform <- sub("_.*", "", type)
      tmp$type <- sub(".*_", "", type)
      all_results[[length(all_results) + 1]] <- tmp
    }
  }
}

all_moran_individual_NC <- bind_rows(all_results)

all_moran_individual_NC$platform2 <- dplyr::case_when(
  grepl("MERSCOPE", all_moran_individual_NC$platform) & (all_moran_individual_NC$sample %in% bench_samples_with_MERSCOPE_v1) ~ "MERSCOPE_V1",
  grepl("MERSCOPE", all_moran_individual_NC$platform) & (all_moran_individual_NC$sample %in% bench_samples_with_MERSCOPE_v2) ~ "MERSCOPE_V2",
  grepl("Xenium", all_moran_individual_NC$platform) ~ "Xenium",
  TRUE ~ "Unknown"
)

all_moran_individual_NC$sample_factor <- factor(all_moran_individual_NC$sample, levels = bench_samples)

saveRDS(all_moran_individual_NC, file.path(results_table_dir, "all_moran_individual_NC.rds"))
