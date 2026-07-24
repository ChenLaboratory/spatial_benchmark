# Cell type composition (proportion of cells per SingleR_labels202605 type)
# per sample and platform (Xenium, MERSCOPE, matched sc/snRNA-seq reference),
# after cross-platform alignment -- i.e. restricted to cells retained in the
# overlapped hexbins used for cross-platform comparison. Feeds Fig5 panel b.
# See README.md.
#
# The source analysis (cell_type_identification.Rmd) also computes the
# equivalent "before alignment" (all cells, no overlap restriction) version;
# not reproduced here since Fig5b only uses "after alignment". See README.md.

library(Seurat)
library(dplyr)

bench_samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")

overlapped_data_file <- "/path/to/analysis/process_aligned_data/overlapped_data_list_AllSample.rds"
xenium_dir   <- "/path/to/preprocessing/Xenium"    # <sample>/so.rds
merscope_dir <- "/path/to/preprocessing/MERSCOPE"  # <sample>/so.rds
sc_dir       <- "/path/to/preprocessing/scRNA"     # <sample>/seu_<sample>.rds (snRNA for MH0007/MH0026)

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

saveRDS(comp_df, "comp_df.rds")
