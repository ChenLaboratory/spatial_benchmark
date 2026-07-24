# Extended Data Fig 5: per-gene total count scatter (Xenium/MERSCOPE vs
# Visium), all 6 benchmarking samples, log scale, Visium area-scaled to
# hexbin size. Top 8 outlier genes (farthest from y=x) labeled; PTPRC,
# NOTCH2, CD68 and WNT3A highlighted as enlarged colored points. Same
# underlying data/scatter as figures/Fig2 f2h (MH0026 only). See README.md.

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrepel)
library(yardstick)

gene_counts_file <- "/path/to/analysis/sensitivity_metrics/gene_counts_per_sample.rds"
out_dir <- "/path/to/results/figures"

bench_samples <- c("MH0007", "MH0026", "ER_0114", "ER_0360", "TN_0177", "TN_0554")
bench_samples_with_MERSCOPE <- c("MH0026", "MH0007", "ER_0360", "ER_0114", "TN_0177")

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

plots <- lapply(bench_samples, gene_scatter_sample, gene_counts_df = gene_counts_per_sample)
figS5 <- cowplot::plot_grid(plotlist = plots, ncol = 1)

ggsave(file.path(out_dir, "figS5_gene_count_scatter_all_samples.pdf"), figS5, width = 15, height = 30, limitsize = FALSE)
