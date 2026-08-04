# Purpose:  Extended Data Fig 5 — per-gene total count scatter, Xenium/MERSCOPE vs
#           Visium, for all six benchmarking samples: log scale, Visium area-scaled to
#           hexbin size. The eight genes furthest from y=x are labelled, and PTPRC,
#           NOTCH2, CD68 and WNT3A are highlighted as enlarged coloured points. Same
#           underlying data and scatter as Fig2 panel h, which shows TNBC_02 alone.
#           See README.md.
# Inputs:   <results_table_dir>/gene_counts_per_sample.rds
#                     written by ../analysis/sensitivity_metrics/03_gene_count_scatter_data.R
# Outputs:  <figure_output_dir>/Ext_Fig5_gene_count_scatter_all_samples.pdf

library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrepel)
library(yardstick)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

bench_samples_with_MERSCOPE <- bench_samples_merscope

gene_counts_file <- file.path(results_table_dir, "gene_counts_per_sample.rds")
out_dir <- figure_output_dir

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
ext_fig5 <- cowplot::plot_grid(plotlist = plots, ncol = 1)

ggsave(file.path(out_dir, "Ext_Fig5_gene_count_scatter_all_samples.pdf"), ext_fig5, width = 15, height = 30, limitsize = FALSE)
