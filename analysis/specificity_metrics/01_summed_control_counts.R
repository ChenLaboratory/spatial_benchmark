# Per-hexbin negative-control-probe and negative-control-codeword counts,
# summed across all probes/codewords per platform (i.e. total off-target-probe
# / decoding-error signal per bin, not broken down by individual probe).
# Feeds Fig3 panels b-h. See README.md.
#
# Xenium's per-cell `control_probe_counts`/`control_codeword_counts` are
# native 10x Xenium QC metadata columns, expected to already be present on
# the Xenium Seurat objects from preprocessing -- not computed here.
# MERSCOPE has no such native column; its control probes/codewords are
# identified by name pattern and summed per cell below.

library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)

overlapped_data_file <- "/path/to/analysis/process_aligned_data/overlapped_data_list_AllSample.rds"

bench_samples <- c("MH0007", "MH0026", "ER_0360", "ER_0114", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")
bench_samples_with_MERSCOPE_v1 <- c("MH0026", "MH0007")
bench_samples_with_MERSCOPE_v2 <- c("ER_0360", "ER_0114", "TN_0177")

# Xenium gene panel has 20 negative control probes and 41 negative control
# codewords (fixed panel design, used as the FDR denominator below)
XENIUM_N_PROBE <- 20
XENIUM_N_CODEWORD <- 41

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)

visium_list <- list()

for (sample in bench_samples) {
  message(sample)
  visium_seu <- overlapped_data_list_AllSample[[sample]]$Visium

  # MERSCOPE: sum control probe/codeword counts per cell, then per hexbin
  if (sample %in% bench_samples_with_MERSCOPE) {
    iST_seu <- overlapped_data_list_AllSample[[sample]]$MERSCOPE
    counts <- as.matrix(GetAssayData(iST_seu, layer = "counts"))
    genes <- rownames(counts)[!(grepl("Blank", rownames(counts)) |
                                 grepl("Sabrina", rownames(counts)) |
                                 grepl("^FC[0-9]+$", rownames(counts)))]
    gene_num <- length(genes)

    if (sample %in% bench_samples_with_MERSCOPE_v1) {
      control_codeword <- rownames(iST_seu)[grepl("Blank", rownames(iST_seu))]
      control_probe <- rownames(iST_seu)[grepl("Sabrina", rownames(iST_seu))]
    } else {
      control_codeword <- rownames(iST_seu)[grepl("Blank", rownames(iST_seu))]
      control_probe <- rownames(iST_seu)[grepl("^FC[0-9]+$", rownames(iST_seu))]
    }

    iST_seu$control_codeword_counts <- Matrix::colSums(counts[control_codeword, ])
    iST_seu$control_probe_counts <- Matrix::colSums(counts[control_probe, ])

    ncount_sum <- iST_seu@meta.data %>%
      group_by(Visium_spot_id) %>%
      summarise(total_control_codeword_counts_MERSCOPE = sum(control_codeword_counts, na.rm = TRUE),
                total_control_probe_counts_MERSCOPE = sum(control_probe_counts, na.rm = TRUE))
    rownames(ncount_sum) <- ncount_sum$Visium_spot_id

    visium_seu@meta.data <- visium_seu@meta.data %>%
      rownames_to_column(var = "Visium_spot_id") %>%
      left_join(ncount_sum, by = "Visium_spot_id") %>%
      column_to_rownames("Visium_spot_id")

    visium_seu$percent_control_codeword_counts_MERSCOPE <- visium_seu$total_control_codeword_counts_MERSCOPE / visium_seu$MERSCOPE_total_nCount_hexbin
    visium_seu$percent_control_probe_counts_MERSCOPE <- visium_seu$total_control_probe_counts_MERSCOPE / visium_seu$MERSCOPE_total_nCount_hexbin
    visium_seu$PerCell_control_codeword_counts_MERSCOPE <- visium_seu$total_control_codeword_counts_MERSCOPE / visium_seu$MERSCOPE_cell_count_filtered
    visium_seu$PerCell_control_probe_counts_MERSCOPE <- visium_seu$total_control_probe_counts_MERSCOPE / visium_seu$MERSCOPE_cell_count_filtered
    visium_seu$FDR_control_codeword_counts_MERSCOPE <- visium_seu$percent_control_codeword_counts_MERSCOPE * gene_num / length(control_codeword) * 100
    visium_seu$FDR_control_probe_counts_MERSCOPE <- visium_seu$percent_control_probe_counts_MERSCOPE * gene_num / length(control_probe) * 100
  }

  # Xenium: control_probe_counts/control_codeword_counts already per-cell
  iST_seu <- overlapped_data_list_AllSample[[sample]]$Xenium
  ncount_sum <- iST_seu@meta.data %>%
    group_by(Visium_spot_id) %>%
    summarise(total_control_codeword_counts_Xenium = sum(control_codeword_counts, na.rm = TRUE),
              total_control_probe_counts_Xenium = sum(control_probe_counts, na.rm = TRUE))
  rownames(ncount_sum) <- ncount_sum$Visium_spot_id

  visium_seu@meta.data <- visium_seu@meta.data %>%
    rownames_to_column(var = "Visium_spot_id") %>%
    left_join(ncount_sum, by = "Visium_spot_id") %>%
    column_to_rownames("Visium_spot_id")

  visium_seu$percent_control_codeword_counts_Xenium <- visium_seu$total_control_codeword_counts_Xenium / visium_seu$Xenium_total_nCount_hexbin
  visium_seu$percent_control_probe_counts_Xenium <- visium_seu$total_control_probe_counts_Xenium / visium_seu$Xenium_total_nCount_hexbin
  visium_seu$PerCell_control_codeword_counts_Xenium <- visium_seu$total_control_codeword_counts_Xenium / visium_seu$Xenium_cell_count_filtered
  visium_seu$PerCell_control_probe_counts_Xenium <- visium_seu$total_control_probe_counts_Xenium / visium_seu$Xenium_cell_count_filtered
  visium_seu$FDR_control_codeword_counts_Xenium <- visium_seu$percent_control_codeword_counts_Xenium * nrow(iST_seu) / XENIUM_N_CODEWORD * 100
  visium_seu$FDR_control_probe_counts_Xenium <- visium_seu$percent_control_probe_counts_Xenium * nrow(iST_seu) / XENIUM_N_PROBE * 100

  visium_seu@meta.data[is.na(visium_seu@meta.data)] <- 0

  visium_list[[sample]] <- visium_seu
}

# per-sample Visium objects (spatial image + summed control columns) --
# feeds the Fig3b spatial plots
saveRDS(visium_list, "visium_summed_control.rds")

# long-format table across samples/platforms, for each NC type (probe,
# codeword) and metric (total, PerCell, percent, FDR) -- feeds Fig3c-h
collect_control_counts <- function(NC, param) {
  All <- list()
  for (sample in bench_samples) {
    meta <- visium_list[[sample]]@meta.data
    meta$Sample <- sample

    cols <- intersect(
      c(paste0(param, "_control_", NC, "_counts_MERSCOPE"),
        paste0(param, "_control_", NC, "_counts_Xenium")),
      colnames(meta)
    )
    if (length(cols) == 0) next

    df_long <- reshape2::melt(meta[, c("Sample", cols), drop = FALSE],
                               id.vars = "Sample", variable.name = "PlatformMeasurement",
                               value.name = "Counts")
    df_long$Platform <- dplyr::case_when(
      grepl("MERSCOPE", df_long$PlatformMeasurement) & sample %in% bench_samples_with_MERSCOPE_v1 ~ "MERSCOPE_V1",
      grepl("MERSCOPE", df_long$PlatformMeasurement) & sample %in% bench_samples_with_MERSCOPE_v2 ~ "MERSCOPE_V2",
      grepl("Xenium", df_long$PlatformMeasurement) ~ "Xenium",
      TRUE ~ "Unknown"
    )
    df_long$NC <- NC
    df_long$metric <- param
    All[[sample]] <- df_long
  }
  dplyr::bind_rows(All)
}

summed_control_counts_long <- dplyr::bind_rows(lapply(
  c("probe", "codeword"), function(nc) {
    dplyr::bind_rows(lapply(c("total", "PerCell", "percent", "FDR"),
                             function(param) collect_control_counts(nc, param)))
  }
))

saveRDS(summed_control_counts_long, "summed_control_counts_long.rds")
