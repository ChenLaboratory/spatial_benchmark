# Purpose:  Cell type composition — proportion of cells per SingleR_labels202605 type —
#           per sample and platform (Xenium, MERSCOPE, and the matched sc/snRNA-seq
#           reference), after cross-platform alignment: restricted to the cells retained
#           in the overlapped hexbins used for cross-platform comparison, so the three
#           platforms are compared over the same tissue. See README.md.
# Inputs:   <processed_data_dir>/overlapped_data_list_AllSample.rds   from ../../integration
#           <processed_data_dir>/<platform>/<sample>/so.rds           annotated iST objects
#           <processed_data_dir>/scRNA/<sample>/seu_<sample>.rds      annotated reference
# Outputs:  <results_table_dir>/comp_df.rds
#                     feeds figures/Fig5_celltype_cci.R panel b
#
# The source analysis (cell_type_identification.Rmd) also computes the
# equivalent "before alignment" (all cells, no overlap restriction) version;
# not reproduced here since Fig5b only uses "after alignment". See README.md.

library(Seurat)
library(dplyr)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

bench_samples_with_MERSCOPE <- bench_samples_merscope

overlapped_data_file <- file.path(processed_data_dir, "overlapped_data_list_AllSample.rds")
xenium_dir   <- file.path(processed_data_dir, "Xenium")    # <sample>/so.rds
merscope_dir <- file.path(processed_data_dir, "MERSCOPE")  # <sample>/so.rds
sc_dir       <- file.path(processed_data_dir, "scRNA")     # <sample>/seu_<sample>.rds (snRNA for TNBC_01/TNBC_02)

overlapped_data_list_AllSample <- readRDS(overlapped_data_file)

get_prop_after_alignment <- function(sample) {
  out_list <- list()

  # scRNA / snRNA reference -- shown as-is, no before/after distinction
  ref_path <- file.path(sc_dir, sample, paste0("seu_", sample, ".rds"))
  ref_seu  <- readRDS(ref_path)
  ref_labels <- ref_seu$cell_type202605
  ref_labels <- ref_labels[ref_labels != "mixture"]
  out_list[["scRNA"]] <- data.frame(celltype = ref_labels, platform = "scRNA", sample = sample, stringsAsFactors = FALSE)

  # Xenium: cells retained in the cross-platform overlapped hexbins
  xen_path <- file.path(xenium_dir, sample, "so.rds")
  if (file.exists(xen_path)) {
    so_xen <- readRDS(xen_path)
    aligned_xen <- overlapped_data_list_AllSample[[sample]]$Xenium
    if (!is.null(aligned_xen) && "SingleR_labels202605" %in% colnames(so_xen@meta.data)) {
      shared_cells <- intersect(colnames(so_xen), colnames(aligned_xen))
      out_list[["Xenium"]] <- data.frame(
        celltype = so_xen$SingleR_labels202605[shared_cells], platform = "Xenium", sample = sample,
        stringsAsFactors = FALSE
      )
    }
  }

  # MERSCOPE: same, where available
  if (sample %in% bench_samples_with_MERSCOPE) {
    mer_path <- file.path(merscope_dir, sample, "so.rds")
    if (file.exists(mer_path)) {
      so_mer <- readRDS(mer_path)
      aligned_mer <- overlapped_data_list_AllSample[[sample]]$MERSCOPE
      if (!is.null(aligned_mer) && "SingleR_labels202605" %in% colnames(so_mer@meta.data)) {
        shared_cells <- intersect(colnames(so_mer), colnames(aligned_mer))
        out_list[["MERSCOPE"]] <- data.frame(
          celltype = so_mer$SingleR_labels202605[shared_cells], platform = "MERSCOPE", sample = sample,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- dplyr::bind_rows(out_list)
  if (nrow(df) == 0) return(NULL)

  df %>%
    dplyr::count(sample, platform, celltype, name = "n") %>%
    dplyr::group_by(sample, platform) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()
}

comp_df <- dplyr::bind_rows(lapply(bench_samples, get_prop_after_alignment))

platform_levels <- c("scRNA", "Xenium", "MERSCOPE")
comp_df$platform <- factor(comp_df$platform, levels = platform_levels)
comp_df$sample   <- factor(comp_df$sample,   levels = bench_samples)

saveRDS(comp_df, file.path(results_table_dir, "comp_df.rds"))
