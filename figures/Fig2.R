# Fig2: bin-level transcript/gene detection sensitivity across Xenium,
# MERSCOPE and Visium.
#   f2a  genes per bin,       all genes,    before alignment  (log10 y)
#   f2b  transcripts per bin, all genes,    before alignment  (log10 y)
#   f2c  genes per bin,       shared genes, before alignment  (log10 y)
#   f2d  transcripts per bin, shared genes, before alignment  (log10 y)
#   f2f  genes per bin,       shared genes, after alignment   (linear y)
#   f2g  transcripts per bin, shared genes, after alignment   (log10 y)
#   f2h  per-gene total count scatter (MH0026), Xenium/MERSCOPE vs Visium,
#        log scale, Visium area-scaled, 4 genes of interest highlighted


library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrepel)
library(yardstick)

# reads the tables produced by `../../analysis/sensitivity_metrics`
before_all_genes_file    <- "/path/to/analysis/sensitivity_metrics/df_bins_summary_110um_all_ST_all_genes_before_alignment.rds"
before_shared_genes_file <- "/path/to/analysis/sensitivity_metrics/df_bins_summary_110um_all_ST_213_shared_genes_before_alignment.rds"
bin_gene_counts_file     <- "/path/to/analysis/sensitivity_metrics/bin_gene_counts_shared_genes_after_alignment.rds"
bin_nCount_file          <- "/path/to/analysis/sensitivity_metrics/bin_nCount_shared_genes_unscaled_after_alignment.rds"
gene_counts_file         <- "/path/to/analysis/sensitivity_metrics/gene_counts_per_sample.rds"
out_dir <- "/path/to/results/figures"

bench_samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")
bench_samples_with_MERSCOPE_v1 <- c("MH0026", "MH0007")
bench_samples_with_MERSCOPE_v2 <- c("ER_0360", "ER_0114", "TN_0177")
platform_order <- c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium", "Visium")

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

platform2_of <- function(platform, sample) {
  dplyr::case_when(
    grepl("MERSCOPE", platform) & (sample %in% bench_samples_with_MERSCOPE_v1) ~ "MERSCOPE_V1",
    grepl("MERSCOPE", platform) & (sample %in% bench_samples_with_MERSCOPE_v2) ~ "MERSCOPE_V2",
    TRUE ~ platform
  )
}

## ---- f2a-d: before alignment ----------------------------------------------

bin_violin <- function(df, y, ylab, title) {
  df <- df %>%
    filter(sample %in% bench_samples) %>%
    mutate(
      sample = factor(sample, levels = bench_samples),
      platform = factor(platform, levels = platform_order)
    )

  ggplot(df, aes(x = platform, y = .data[[y]], fill = platform)) +
    geom_violin(scale = "width", trim = TRUE, color = "black", size = 0.3) +
    geom_boxplot(width = 0.1, fill = "white", color = "black", size = 0.3, outlier.size = 0.3) +
    facet_wrap(~sample, nrow = 1, scales = "fixed") +
    scale_y_log10() +
    scale_fill_manual(values = platform_colors) +
    labs(x = "", y = ylab, title = title, fill = "Platform") +
    theme_common
}

df_all_genes <- readRDS(before_all_genes_file)
df_shared_genes <- readRDS(before_shared_genes_file)

f2a <- bin_violin(df_all_genes, "n_genes", "Genes per bin (log10)",
                   "Detected genes per bin (all genes)")
f2b <- bin_violin(df_all_genes, "n_transcripts", "Transcripts per bin (log10)",
                   "Transcripts per bin (all genes)")
f2c <- bin_violin(df_shared_genes, "n_genes", "Genes per bin (log10)",
                   "Detected genes per bin (shared genes)")
f2d <- bin_violin(df_shared_genes, "n_transcripts", "Transcripts per bin (log10)",
                   "Transcripts per bin (shared genes)")

## ---- f2f: genes per bin, shared genes, after alignment (linear y) --------

bin_gene_counts <- readRDS(bin_gene_counts_file) %>%
  filter(sample %in% bench_samples) %>%
  mutate(
    platform2 = platform2_of(Platform, sample),
    platform2 = factor(platform2, levels = platform_order),
    sample = factor(sample, levels = bench_samples)
  )

f2f <- ggplot(bin_gene_counts, aes(x = platform2, y = Unique_Genes, fill = platform2)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_grid(~sample, scales = "fixed", drop = TRUE) +
  scale_fill_manual(values = platform_colors, drop = TRUE) +
  labs(x = "", y = "Gene counts per hexbin", title = "Detected genes per bin (shared genes)", fill = "Platform") +
  theme_common

## ---- f2g: transcripts per bin, shared genes, after alignment (log10 y) ---
## unscaled Visium values (no area rescaling)

bin_nCount <- readRDS(bin_nCount_file) %>%
  filter(sample %in% bench_samples) %>%
  mutate(
    platform2 = platform2_of(platform, sample),
    platform2 = factor(platform2, levels = platform_order),
    sample = factor(sample, levels = bench_samples)
  )

f2g <- ggplot(bin_nCount, aes(x = platform2, y = value, fill = platform2)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_grid(~sample, scales = "fixed", drop = TRUE) +
  scale_fill_manual(values = platform_colors, drop = TRUE) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "", y = "Total counts per hexbin (log10)",
       title = "Total transcript counts per bin (shared genes, unscaled Visium)", fill = "Platform") +
  theme_common

## ---- f2h: gene count scatter (MH0026), Xenium/MERSCOPE vs Visium ---------
## total counts per gene, log scale, Visium area-scaled to hexbin size; top 8
## outlier genes (farthest from y=x) labeled; PTPRC/NOTCH2/CD68/WNT3A
## highlighted as enlarged colored points (no extra text label for these)

highlight_genes <- c("PTPRC", "NOTCH2", "CD68", "WNT3A")
highlight_colors <- c(
  "PTPRC"  = "#E41A1C",
  "NOTCH2" = "#377EB8",
  "CD68"   = "#4DAF4A",
  "WNT3A"  = "#FF7F00"
)

gene_scatter_spearman <- function(df, meta1, meta2, point_highlight, point_highlight_colors) {
  var1 <- log(df[[meta1]] + 0.5)
  var2 <- log(df[[meta2]] + 0.5)
  keep <- !is.na(var1) & !is.na(var2)
  df_tmp <- data.frame(x = var1[keep], y = var2[keep], gene = df$genes[keep])
  ct <- cor.test(df_tmp$x, df_tmp$y, method = "spearman")
  rho_lab <- paste0("Spearman rho = ", round(ct$estimate, 3),
                     ", p = ", formatC(ct$p.value, format = "e", digits = 2))
  df_tmp$diff <- abs(df_tmp$y - df_tmp$x)
  outliers <- head(df_tmp[order(-df_tmp$diff), ], 8)

  hi_pts <- df_tmp[df_tmp$gene %in% point_highlight, ]
  hi_pts$pt_color <- point_highlight_colors[hi_pts$gene]

  ggplot(df_tmp, aes(x = x, y = y)) +
    geom_point(size = 0.6, color = "#2f4c7a") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text_repel(data = outliers, aes(x = x, y = y, label = gene), inherit.aes = FALSE,
                     size = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4, box.padding = 0.4) +
    geom_point(data = hi_pts, aes(x = x, y = y, color = pt_color), inherit.aes = FALSE, size = 1.5) +
    scale_color_identity(guide = "legend", name = NULL,
                          breaks = point_highlight_colors, labels = names(point_highlight_colors)) +
    labs(subtitle = rho_lab, x = paste0("log(", meta1, "+0.5)"), y = paste0("log(", meta2, "+0.5)")) +
    yardstick::coord_obs_pred() +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

gene_scatter_sample <- function(sample, gene_counts_df) {
  df <- gene_counts_df[gene_counts_df$sample == sample, ]
  df$Visium <- df$Visium_scaled
  platforms <- if (sample %in% bench_samples_with_MERSCOPE) c("Xenium", "MERSCOPE") else "Xenium"

  title <- ggdraw() + draw_label(paste0(sample, " (Visium scaled)"), fontface = "bold")
  plots <- list()
  plots[["XV"]] <- gene_scatter_spearman(df, "Xenium", "Visium", highlight_genes, highlight_colors)
  if ("MERSCOPE" %in% platforms) {
    plots[["XM"]] <- gene_scatter_spearman(df, "Xenium", "MERSCOPE", highlight_genes, highlight_colors)
    plots[["MV"]] <- gene_scatter_spearman(df, "MERSCOPE", "Visium", highlight_genes, highlight_colors)
  }
  p <- cowplot::plot_grid(plotlist = plots, ncol = 3)
  cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
}

gene_counts_per_sample <- readRDS(gene_counts_file)
f2h <- gene_scatter_sample("MH0026", gene_counts_per_sample)
