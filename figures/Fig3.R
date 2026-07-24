# Fig3: evaluation of gene detection specificity (negative controls).
#   b   spatial density of summed NC probe/codeword counts, TNBC-02 (MH0026),
#       Xenium + MERSCOPE V1, 110um hexbins (log10 colour scale)
#   c,d violin: summed NC probe (c) / codeword (d) counts per 110um bin
#   e,f violin: summed NC probe (e) / codeword (f) counts per cell in bin
#   g,h violin: FDR of NC probe (g) / codeword (h) per 110um bin
#   i   boxplot: Global Moran's I for individual NC probes vs individual NC
#       codewords, log10(x+0.01), coloured by significance
#   j   spatial plots of individual NC probes/codewords, TNBC-02 (MH0026),
#       Xenium + MERSCOPE V1, 50um hexbins
#   k   scatter: Global Moran's I vs total signal count, TNBC-02 (MH0026),
#       gene-targeting probes (blue) vs NC probes (red) vs NC codewords
#       (orange)
# "TNBC-02" = MH0026 (the only sample with Xenium + MERSCOPE V1 in the
# aligned data). See README.md.

library(ggplot2)
library(dplyr)
library(cowplot)
library(scider)
library(Seurat)

visium_summed_control_file    <- "/path/to/analysis/specificity_metrics/visium_summed_control.rds"
summed_control_counts_file    <- "/path/to/analysis/specificity_metrics/summed_control_counts_long.rds"
moran_individual_nc_file      <- "/path/to/analysis/specificity_metrics/all_moran_individual_NC.rds"
individual_control_spatial_file <- "/path/to/analysis/specificity_metrics/individual_control_spatial.rds"
moran_genes_file              <- "/path/to/analysis/specificity_metrics/all_moran_genes.rds"
total_counts_file             <- "/path/to/analysis/specificity_metrics/total_counts_df.rds"
out_dir <- "/path/to/results/figures"

bench_samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")
tnbc02_sample <- "MH0026"  # "TNBC-02" in the manuscript

platform_colors <- c(
  "MERSCOPE_V1" = "#80B1D3",
  "MERSCOPE_V2" = "#0b6170",
  "Xenium"      = "#FDAE61",
  "Visium"      = "#BC80BD"
)

theme_common <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.3),
    panel.grid.minor = element_line(color = "grey92", size = 0.2),
    panel.border = element_rect(color = "black", size = 0.3),
    strip.background = element_rect(fill = "white", color = "black", size = 0.3),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(size = 13, hjust = 0.5),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key = element_rect(fill = "white", color = NA),
    panel.background = element_blank(),
    legend.background = element_blank()
  )

# white-background theme used for the scatter panels (k), matching
# moran_all_gene.Rmd's theme_light_scatter
theme_light_scatter <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor  = element_line(color = "grey92", linewidth = 0.2),
    panel.border      = element_rect(color = "black", linewidth = 0.3),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text        = element_text(size = 10, face = "bold"),
    axis.text         = element_text(size = 9),
    axis.title        = element_text(size = 11),
    plot.title        = element_text(hjust = 0.5, size = 13),
    legend.background = element_rect(fill = "white"),
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 10),
    legend.key        = element_rect(fill = "white", color = NA)
  )

format_counts <- function(x) {
  dplyr::case_when(x >= 1e6 ~ paste0(x / 1e6, "M"), x >= 1e3 ~ paste0(x / 1e3, "k"), TRUE ~ as.character(x))
}

## ---- b: spatial density, summed NC counts, TNBC-02 (log10) ---------------

visium_summed_control <- readRDS(visium_summed_control_file)

control_spatial_plot <- function(sample, visium_list, var, psf = 3) {
  visium_seu <- visium_list[[sample]]
  feat <- paste0("total_", var, "_counts")
  platforms <- c("Xenium", "MERSCOPE")[paste0(feat, "_", c("Xenium", "MERSCOPE")) %in% colnames(visium_seu@meta.data)]
  feat_cols <- paste0(feat, "_", platforms)
  global_limits <- range(unlist(visium_seu[[feat_cols]]), na.rm = TRUE)

  visium_seu[[feat_cols]] <- visium_seu[[feat_cols]] + 1
  plots <- lapply(platforms, function(plat) {
    SpatialFeaturePlot(visium_seu, features = paste0(feat, "_", plat), crop = TRUE,
                        pt.size.factor = psf, stroke = NA, image.alpha = 0) +
      scale_fill_viridis_c(option = "viridis", limits = global_limits + 1, trans = "log10",
                            breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      ggtitle(plat) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(fill = feat)
  })

  p <- cowplot::plot_grid(plotlist = plots, ncol = 2)
  title <- ggdraw() + draw_label(paste0(sample, " ", var), fontface = "bold", size = 18)
  cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
}

f3b_probe <- control_spatial_plot(tnbc02_sample, visium_summed_control, "control_probe", psf = 4)
f3b_codeword <- control_spatial_plot(tnbc02_sample, visium_summed_control, "control_codeword", psf = 4)

## ---- c-h: violin plots, summed NC counts per bin --------------------------

summed_control_counts <- readRDS(summed_control_counts_file)

nc_violin <- function(df, nc, metric, ylab, title) {
  df <- df %>%
    filter(NC == nc, metric == !!metric, Platform != "Unknown") %>%
    mutate(Sample = factor(Sample, levels = bench_samples))

  ggplot(df, aes(x = Platform, y = Counts, fill = Platform)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
    facet_grid(~ Sample, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = platform_colors, drop = TRUE) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    xlab("") +
    labs(title = title, y = ylab) +
    theme_common
}

f3c <- nc_violin(summed_control_counts, "probe", "total", "NC probe counts per bin (log10)",
                  "NC probe counts per 110um bin (Log10)")
f3d <- nc_violin(summed_control_counts, "codeword", "total", "NC codeword counts per bin (log10)",
                  "NC codeword counts per 110um bin (Log10)")
f3e <- nc_violin(summed_control_counts, "probe", "PerCell", "NC probe counts per cell in bin (log10)",
                  "NC probe counts per cell in bin (Log10)")
f3f <- nc_violin(summed_control_counts, "codeword", "PerCell", "NC codeword counts per cell in bin (log10)",
                  "NC codeword counts per cell in bin (Log10)")
f3g <- nc_violin(summed_control_counts, "probe", "FDR", "NC probe FDR per bin (log10)",
                  "NC probe FDR per 110um bin (Log10)")
f3h <- nc_violin(summed_control_counts, "codeword", "FDR", "NC codeword FDR per bin (log10)",
                  "NC codeword FDR per 110um bin (Log10)")

## ---- i: boxplot, Global Moran's I for individual NC probes/codewords -----
## log10(x+0.01), coloured by significance (controls only -- no gene panel)

all_moran_individual_NC <- readRDS(moran_individual_nc_file) %>%
  mutate(sig = pval < 0.05, sample_factor = factor(sample, levels = bench_samples))

all_moran_genes <- readRDS(moran_genes_file) %>%
  mutate(sample_factor = factor(sample, levels = bench_samples))

combined_log_df <- all_moran_individual_NC %>%
  mutate(data_type = case_when(type == "probe" ~ "Negative Control Probes",
                                type == "codeword" ~ "Negative Control Codewords")) %>%
  select(platform2, lisa, sig, sample_factor, data_type) %>%
  mutate(
    data_type = factor(data_type, levels = c("Negative Control Probes", "Negative Control Codewords")),
    sample_factor = factor(sample_factor, levels = bench_samples),
    platform2 = factor(platform2, levels = c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium"))
  )

# log10(x + 0.01): keeps near-zero/negative LISA values on the scale
log10p_scale <- scale_y_continuous(
  trans = scales::trans_new("log10p",
    transform = function(x) log10(x + 0.01), inverse = function(x) 10^x - 0.01),
  breaks = c(-0.007, -0.005, 0, 0.01, 0.05, 0.1, 0.2),
  labels = c("-0.007", "-0.005", "0", "0.01", "0.05", "0.1", "0.2")
)

f3i <- ggplot(combined_log_df, aes(x = platform2, y = lisa)) +
  geom_boxplot(aes(fill = platform2), width = 0.6, outlier.shape = NA, linewidth = 0.15, na.rm = TRUE) +
  geom_jitter(aes(color = sig), width = 0.15, size = 0.2, na.rm = TRUE) +
  facet_grid(~ data_type + sample_factor, scales = "free_x", space = "free_x") +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = platform_colors, drop = FALSE) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red"), name = "Significant\n(p < 0.05)") +
  log10p_scale +
  labs(x = "Platform", y = "Moran's I (LISA, log10(x+0.01))") +
  theme_common +
  theme(strip.text = element_text(size = 8, face = "plain"), strip.text.x.top = element_text(face = "bold"))

## ---- j: spatial plots, individual NC probes/codewords, TNBC-02 -----------
## (all individual NC features are plotted here; the manuscript panel shows a
## manually curated subset -- see README.md)

individual_control_spatial <- readRDS(individual_control_spatial_file)

nc_feature_grid <- function(spe_hex, neg_controls, title) {
  plots <- lapply(neg_controls, function(neg) {
    scider::plotGrid(spe_hex, feature = neg) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                    barheight = unit(2, "mm"), barwidth = unit(20, "mm"))) +
      theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size = 6),
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.line = element_blank())
  })
  combined <- cowplot::plot_grid(plotlist = plots, ncol = 7)
  title_grob <- ggdraw() + draw_label(title, fontface = "bold", size = 14)
  cowplot::plot_grid(title_grob, combined, ncol = 1, rel_heights = c(0.05, 1))
}

f3j_xenium <- nc_feature_grid(individual_control_spatial[[tnbc02_sample]]$Xenium$spe_hex,
                               individual_control_spatial[[tnbc02_sample]]$Xenium$neg_controls,
                               paste0(tnbc02_sample, " - Xenium Neg Controls"))
f3j_merscope <- nc_feature_grid(individual_control_spatial[[tnbc02_sample]]$MERSCOPE$spe_hex,
                                 individual_control_spatial[[tnbc02_sample]]$MERSCOPE$neg_controls,
                                 paste0(tnbc02_sample, " - MERSCOPE Neg Controls"))

## ---- k: scatter, Global Moran's I vs total signal, TNBC-02 ---------------

category_colors <- c(
  "Pre-designed Genes"        = "#8ebbde",
  "Negative Control Probes"   = "#ed0524",
  "Negative Control Codeword" = "#e87e27"
)

total_counts_df <- readRDS(total_counts_file)

scatter_genes <- all_moran_genes %>%
  mutate(feature_clean = make.names(gene), category = "Pre-designed Genes") %>%
  select(feature_clean, lisa, sample, platform2, category)

scatter_nc <- all_moran_individual_NC %>%
  mutate(
    feature_clean = make.names(sub("^(Xenium|MERSCOPE)_", "", feature)),
    category = case_when(type == "probe" ~ "Negative Control Probes",
                          type == "codeword" ~ "Negative Control Codeword",
                          TRUE ~ NA_character_)
  ) %>%
  select(feature_clean, lisa, sample, platform2, category)

scatter_all <- bind_rows(scatter_genes, scatter_nc) %>%
  left_join(total_counts_df, by = c("feature_clean", "sample", "platform2")) %>%
  filter(sample == tnbc02_sample) %>%
  mutate(
    category = factor(category, levels = names(category_colors)),
    platform2 = factor(platform2, levels = c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium"))
  ) %>%
  arrange(category)  # plot controls on top

f3k <- ggplot(scatter_all, aes(x = total_counts, y = lisa, color = category)) +
  geom_point(size = 0.8, alpha = 0.9) +
  scale_x_log10(breaks = c(10, 100, 1e3, 1e4, 1e5, 1e6), labels = format_counts) +
  scale_color_manual(values = category_colors, name = "Gene Category") +
  facet_wrap(~platform2, nrow = 1, drop = TRUE) +
  labs(title = paste0("Global Moran's I vs. Transcript Counts - ", tnbc02_sample),
       x = "Transcripts per gene", y = "Moran's I (LISA)") +
  theme_light_scatter +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

## ---- save ------------------------------------------------------------------

ggsave(file.path(out_dir, "f3b_spatial_density_NC_probe_TNBC02.pdf"), f3b_probe, width = 10, height = 6)
ggsave(file.path(out_dir, "f3b_spatial_density_NC_codeword_TNBC02.pdf"), f3b_codeword, width = 10, height = 6)
ggsave(file.path(out_dir, "f3c_NC_probe_per_bin.pdf"), f3c, width = 9, height = 3)
ggsave(file.path(out_dir, "f3d_NC_codeword_per_bin.pdf"), f3d, width = 9, height = 3)
ggsave(file.path(out_dir, "f3e_NC_probe_per_cell_in_bin.pdf"), f3e, width = 9, height = 3)
ggsave(file.path(out_dir, "f3f_NC_codeword_per_cell_in_bin.pdf"), f3f, width = 9, height = 3)
ggsave(file.path(out_dir, "f3g_NC_probe_FDR_per_bin.pdf"), f3g, width = 9, height = 3)
ggsave(file.path(out_dir, "f3h_NC_codeword_FDR_per_bin.pdf"), f3h, width = 9, height = 3)
ggsave(file.path(out_dir, "f3i_moran_boxplot_genes_vs_NC.pdf"), f3i, width = 24, height = 7)
ggsave(file.path(out_dir, "f3j_NC_spatial_Xenium_TNBC02.pdf"), f3j_xenium, width = 14, height = 4)
ggsave(file.path(out_dir, "f3j_NC_spatial_MERSCOPE_TNBC02.pdf"), f3j_merscope, width = 14, height = 4)
ggsave(file.path(out_dir, "f3k_moran_vs_counts_TNBC02.pdf"), f3k, width = 9, height = 4)
