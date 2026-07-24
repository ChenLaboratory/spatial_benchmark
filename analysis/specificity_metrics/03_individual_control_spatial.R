# Cell-level negative control probe/codeword counts, hex-binned at 50um, for
# Xenium (all 6 samples) and MERSCOPE (the 5 samples with MERSCOPE) each
# separately. Feeds Fig3 panel j (MH0026 = "TNBC-02") and
# figS6_gene_specificity panel a (the other 5 samples). See README.md.
#
# Unlike the other specificity_metrics scripts (which bin to the 110um
# Visium-matched hexbins), this uses each platform's own finer 50um bins, as
# in the source analysis -- this panel is not cross-platform-aligned, it
# shows each platform's own raw spatial pattern per sample.

library(Seurat)
library(SpatialExperiment)
library(scider)
library(dplyr)
library(Matrix)
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

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)

ST_so2spe <- function(so) {
  counts <- GetAssayData(so, slot = "counts")
  spatial_coords <- so@reductions$spatial@cell.embeddings
  colData <- so@meta.data
  SpatialExperiment(assays = list(counts = counts), colData = colData, spatialCoords = spatial_coords)
}

BinSPE <- function(spe, grid.type = "hex", grid.length.x = 50, ct_group = "orig.ident") {
  ct1 <- spe@colData[[ct_group]]
  coi <- as.character(unique(ct1))
  spe_hex_cell <- scider::gridDensity(spe, coi = coi, grid.type = grid.type, grid.length.x = grid.length.x, id = ct_group)
  spe_hex <- scider::gridSPE(spe_hex_cell, cell.count = TRUE, id = ct_group)
  spe_hex[, spe_hex$cell_counts_overall > 0]
}

individual_control_spatial <- list()

for (sample in bench_samples) {
  message(sample)

  ## ---- Xenium: pull NegControl counts from raw transcripts, hex-bin ------

  iST_seu <- overlapped_data_list_AllSample[[sample]]$Xenium

  df <- read_parquet(xenium_raw_transcripts_files[sample])
  all_cells <- df %>% filter(cell_id != "UNASSIGNED") %>% pull(cell_id) %>% unique()
  neg_counts <- df %>%
    filter(grepl("NegControl", feature_name), cell_id != "UNASSIGNED") %>%
    group_by(feature_name, cell_id) %>%
    summarise(count = n(), .groups = "drop") %>%
    complete(feature_name, cell_id = all_cells, fill = list(count = 0))
  neg_matrix <- neg_counts %>%
    pivot_wider(names_from = cell_id, values_from = count, values_fill = 0) %>%
    as.data.frame()
  rownames(neg_matrix) <- neg_matrix$feature_name
  neg_matrix$feature_name <- NULL
  neg_matrix <- as.data.frame(t(neg_matrix[, colnames(iST_seu)]))
  neg_controls_xenium <- colnames(neg_matrix)

  mat <- GetAssayData(iST_seu, slot = "counts")
  for (neg in neg_controls_xenium) {
    neg_row <- Matrix(neg_matrix[[neg]], nrow = 1, sparse = TRUE)
    rownames(neg_row) <- neg
    mat <- rbind(mat, neg_row)
  }

  spe_xenium <- SpatialExperiment(
    assays = list(counts = mat),
    colData = iST_seu@meta.data,
    spatialCoords = iST_seu@reductions$spatial@cell.embeddings
  )
  spe_hex_xenium <- BinSPE(spe_xenium, grid.type = "hex", grid.length.x = 50, ct_group = "orig.ident")

  sample_result <- list(Xenium = list(spe_hex = spe_hex_xenium, neg_controls = neg_controls_xenium))

  ## ---- MERSCOPE: control rows already in the count matrix ----------------

  if (sample %in% bench_samples_with_MERSCOPE) {
    iST_seu <- overlapped_data_list_AllSample[[sample]]$MERSCOPE
    neg_controls_merscope <- rownames(iST_seu)[grepl("Sabrina", rownames(iST_seu)) |
                                                grepl("^FC[0-9]+$", rownames(iST_seu)) |
                                                grepl("Blank", rownames(iST_seu))]
    spe_merscope <- ST_so2spe(iST_seu)
    spe_hex_merscope <- BinSPE(spe_merscope, grid.type = "hex", grid.length.x = 50, ct_group = "orig.ident")

    sample_result$MERSCOPE <- list(spe_hex = spe_hex_merscope, neg_controls = neg_controls_merscope)
  }

  individual_control_spatial[[sample]] <- sample_result
}

saveRDS(individual_control_spatial, "individual_control_spatial.rds")
