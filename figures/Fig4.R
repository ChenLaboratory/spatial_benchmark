# Fig4: comparing cell segmentation methods (Default, Baysor, Proseg v3).
#   b  dot plot: total cell count, median cell area, transcript assignment
#      rate, transcripts per cell -- across segmentation methods, Xenium
#      (top) / MERSCOPE (bottom), dots only (no box), linked + coloured by
#      sample
#   c  dot plot: median MECR (mutually exclusive co-expression rate) across
#      segmentation methods, faceted by platform (Xenium / MERSCOPE / scRNA),
#      linked + coloured by sample, after QC
#   e  UMAP of MH0007 across segmentation methods x platform, coloured by
#      cell type
#   f  stacked bar: cell type composition for the MH0007 datasets shown in e
#   g  dot plot: Average Silhouette Width (annotated cell type grouping,
#      PCA space), across segmentation methods, faceted by platform, linked
#      + coloured by sample
# See README.md.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(Seurat)

qc_df_file          <- "/path/to/analysis/cell_segmentation_metrics/qc_df.rds"
mecr_summary_file    <- "/path/to/analysis/cell_segmentation_metrics/mecr_summary.rds"
props_all_file       <- "/path/to/analysis/cell_segmentation_metrics/props_all.rds"
sil_ct_results_file  <- "/path/to/analysis/cell_segmentation_metrics/sil_ct_results.rds"
out_dir <- "/path/to/results/figures"

# MH0007 datasets for panels e/f: annotated (SingleR_labels202605) Seurat
# objects, one per segmentation method x platform
mh0007_paths <- c(
  "Default - Xenium"    = "/path/to/preprocessing/Xenium/MH0007/so.rds",
  "Baysor - Xenium"     = "/path/to/run_cell_segmentation/baysor/results/xenium/MH0007_so.rds",
  "Proseg v3 - Xenium"  = "/path/to/run_cell_segmentation/proseg_v3/results/xenium/MH0007_so.rds",
  "Default - MERSCOPE"  = "/path/to/preprocessing/MERSCOPE/MH0007/so.rds",
  "Baysor - MERSCOPE"   = "/path/to/run_cell_segmentation/baysor/results/merscope/MH0007_so.rds",
  "Proseg v3 - MERSCOPE"= "/path/to/run_cell_segmentation/proseg_v3/results/merscope/MH0007_so.rds"
)

fill_colors <- c("Default" = "#1b9e77", "Baysor" = "#d95f02", "Proseg" = "#7570b3")
sample_pal <- c(MH0026 = "#C0392B", MH0007 = "#1565C0", ER0114 = "#F9A825",
                 ER0360 = "#AD1457", TN0177 = "#2E7D32", TN0554 = "#4E342E")
platform_shapes <- c(Xenium = 16, MERSCOPE = 17, scRNA = 15)

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.5),
    panel.border     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    panel.spacing    = unit(1, "lines"),
    legend.position  = "right",
    legend.key.size  = unit(0.4, "cm"),
    axis.text.x      = element_text(angle = 30, hjust = 1)
  )

## ---- b: QC metrics dot plot, Xenium (top) / MERSCOPE (bottom) ------------

qc_df <- readRDS(qc_df_file)

seg_lvls <- c("Default", "Baysor", "Proseg")
plt_lvls <- c("Xenium", "MERSCOPE")

metric_labels_b <- c(
  total_cells                 = "Total cells (K)",
  median_cell_area_um2        = "Median cell area (µm²)",
  pct_assigned                = "% assigned",
  median_transcripts_per_cell = "Median tx/cell"
)

qc_long <- qc_df %>%
  mutate(
    pct_assigned = round(100 * assigned_transcripts / total_transcripts, 1),
    total_cells  = total_cells / 1e3,
    segmentation = factor(segmentation, levels = seg_lvls),
    platform     = factor(platform, levels = plt_lvls)
  ) %>%
  pivot_longer(cols = all_of(names(metric_labels_b)), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = names(metric_labels_b), labels = unname(metric_labels_b)))

y_ranges_b <- qc_long %>%
  group_by(metric) %>%
  summarise(ymin = min(value, na.rm = TRUE), ymax = max(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(ymin, ymax), names_to = NULL, values_to = "value")

f4b <- ggplot(qc_long, aes(x = segmentation, y = value)) +
  geom_blank(data = y_ranges_b, aes(y = value), inherit.aes = FALSE) +
  geom_line(aes(group = sample, color = sample), alpha = 0.55, linewidth = 0.7) +
  geom_point(aes(color = sample), size = 2.5, alpha = 0.9) +
  facet_grid(platform ~ metric, scales = "free") +
  scale_color_manual(values = sample_pal) +
  labs(x = NULL, y = NULL, color = "Sample",
       title = "QC metrics by segmentation method (dots linked per sample)") +
  theme_pub

## ---- c: median MECR dot plot, after QC ------------------------------------

mecr_summary <- readRDS(mecr_summary_file)
y_rng_med <- range(mecr_summary$median_mecr, na.rm = TRUE)

f4c <- mecr_summary %>%
  filter(qc_status == "After QC") %>%
  ggplot(aes(x = segmentation, y = median_mecr, colour = sample, group = sample)) +
  geom_line(linewidth = 0.6, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.9) +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = sample_pal) +
  coord_cartesian(ylim = y_rng_med) +
  labs(x = NULL, y = "Median MECR", colour = "Sample", title = "After QC") +
  theme_pub

## ---- e: UMAP, MH0007, by segmentation method x platform -------------------

cell_type_colors <- c(
  Tumor = "#1F77B4", Myeloid = "#AD8BC9", T = "#D62728", Endothelial = "#8C564B",
  Pericyte = "#C49C94", Fibroblast = "#FF7F0E", Adipocyte = "#E377C2", B = "#BCBD22",
  Plasma = "#17BECF", Epithelial = "#67BF5C", Lymph.vessel = "#A2A2A2",
  unsure = "black", mixture = "grey50"
)
ct_order <- c("Tumor", "Myeloid", "T", "Endothelial", "Pericyte", "Fibroblast",
              "Adipocyte", "B", "Plasma", "Epithelial", "Lymph.vessel", "unsure", "mixture")

theme_umap_ms <- theme_classic(base_size = 10) +
  theme(
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_blank(), axis.text = element_blank(),
    axis.title = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8), legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size = 9, face = "bold")
  )

make_umap <- function(path, title) {
  so <- readRDS(path)
  ct_present <- intersect(ct_order, unique(so$SingleR_labels202605))
  DimPlot(so, group.by = "SingleR_labels202605", label = FALSE,
          cols = cell_type_colors[ct_present], raster = TRUE, raster.dpi = c(2048, 2048), pt.size = 6) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2", color = "Cell type") +
    theme_umap_ms
}

mh0007_umaps <- lapply(names(mh0007_paths), function(nm) make_umap(mh0007_paths[nm], nm))
f4e <- wrap_plots(mh0007_umaps[1:3], nrow = 1) / wrap_plots(mh0007_umaps[4:6], nrow = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

## ---- f: MH0007 cell type composition, stacked bar -------------------------

props_all <- readRDS(props_all_file)

df_mh0007 <- props_all %>% filter(sample_id == "MH0007")
ct_present <- intersect(ct_order, unique(df_mh0007$cell_type))
extra      <- setdiff(unique(df_mh0007$cell_type), ct_order)
ct_levels  <- c(ct_present, extra)

df_mh0007$cell_type <- factor(df_mh0007$cell_type, levels = rev(ct_levels))
cols <- cell_type_colors[ct_levels]
cols[is.na(cols)] <- "grey70"
names(cols) <- ct_levels
cols <- rev(cols)

f4f <- ggplot(df_mh0007, aes(x = method, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = cols, drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  facet_wrap(~platform) +
  labs(title = "MH0007 - cell type composition by segmentation method",
       x = NULL, y = "Proportion", fill = "Cell type") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

## ---- g: Average Silhouette Width (cell type), dot plot -------------------

sil_ct_results <- readRDS(sil_ct_results_file)

sil_ct_plot <- sil_ct_results %>%
  mutate(
    segmentation = dplyr::recode(as.character(method), "default" = "Default", "baysor" = "Baysor", "proseg_v3" = "Proseg"),
    segmentation = factor(segmentation, levels = c("Default", "Baysor", "Proseg")),
    sample_key   = gsub("_", "", as.character(sample_id))
  )
y_rng_ct <- range(sil_ct_plot$value, na.rm = TRUE) * c(0.95, 1.05)

f4g <- ggplot(sil_ct_plot, aes(x = segmentation, y = value, colour = sample_key)) +
  geom_line(aes(group = sample_id), linewidth = 0.6, alpha = 0.7) +
  geom_point(aes(shape = platform), size = 3, alpha = 0.9) +
  scale_colour_manual(values = sample_pal) +
  scale_shape_manual(values = platform_shapes) +
  facet_grid(. ~ platform, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = y_rng_ct) +
  labs(x = NULL, y = "Average Silhouette Width (cell type)", colour = "Sample", shape = "Platform") +
  theme_pub

## ---- save ------------------------------------------------------------------

ggsave(file.path(out_dir, "f4b_qc_metrics_dot.pdf"), f4b, width = 14, height = 6)
ggsave(file.path(out_dir, "f4c_mecr_dot_after_qc.pdf"), f4c, width = 9, height = 4)
ggsave(file.path(out_dir, "f4e_umap_MH0007.pdf"), f4e, width = 12, height = 8)
ggsave(file.path(out_dir, "f4f_cell_type_composition_MH0007.pdf"), f4f, width = 8, height = 5)
ggsave(file.path(out_dir, "f4g_silhouette_width_celltype_dot.pdf"), f4g, width = 6, height = 4)
