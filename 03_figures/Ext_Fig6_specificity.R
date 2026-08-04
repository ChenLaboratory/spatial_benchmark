# Purpose:  Extended Data Fig 6 — supplementary gene detection specificity, all samples.
#             a  spatial plots of individual NC probes/codewords, 50um bins, for the
#                five benchmarking samples not shown in Fig3 panel j (which shows TNBC_02)
#             b  scatter: Global Moran's I vs total signal count, shared gene panel only
#                (n=213) plus all negative controls merged in red, every sample x
#                platform, with the 75th-percentile NC-count threshold shaded and labelled
#           See README.md.
# Inputs:   <results_table_dir>/   written by ../analysis/specificity_metrics/
#             individual_control_spatial.rds, all_moran_individual_NC.rds
#             all_moran_genes.rds, total_counts_df.rds
#           <processed_data_dir>/shared_genes_withVisium.rds   from ../integration
# Outputs:  <figure_output_dir>/
#             f6a_NC_spatial_Xenium_<sample>.pdf     one per sample
#             f6a_NC_spatial_MERSCOPE_<sample>.pdf   one per sample
#             f6b_moran_vs_counts_shared_genes_all_samples.pdf

library(ggplot2)
library(dplyr)
library(cowplot)
library(scider)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))


individual_control_spatial_file <- file.path(results_table_dir, "individual_control_spatial.rds")
moran_individual_nc_file        <- file.path(results_table_dir, "all_moran_individual_NC.rds")
moran_genes_file                <- file.path(results_table_dir, "all_moran_genes.rds")
total_counts_file               <- file.path(results_table_dir, "total_counts_df.rds")
shared_genes_file               <- file.path(processed_data_dir, "shared_genes_withVisium.rds")
out_dir <- figure_output_dir

tnbc02_sample <- "TNBC_02"  # already shown in Fig3j -- excluded from panel a here

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

## ---- a: spatial plots of individual NC probes/codewords, other samples ---

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

f6a_samples <- setdiff(bench_samples, tnbc02_sample)

f6a_xenium <- lapply(f6a_samples, function(sample) {
  nc_feature_grid(individual_control_spatial[[sample]]$Xenium$spe_hex,
                   individual_control_spatial[[sample]]$Xenium$neg_controls,
                   paste0(sample, " - Xenium Neg Controls"))
})
names(f6a_xenium) <- f6a_samples

f6a_merscope_samples <- f6a_samples[sapply(f6a_samples, function(s) !is.null(individual_control_spatial[[s]]$MERSCOPE))]
f6a_merscope <- lapply(f6a_merscope_samples, function(sample) {
  nc_feature_grid(individual_control_spatial[[sample]]$MERSCOPE$spe_hex,
                   individual_control_spatial[[sample]]$MERSCOPE$neg_controls,
                   paste0(sample, " - MERSCOPE Neg Controls"))
})
names(f6a_merscope) <- f6a_merscope_samples

## ---- b: Moran's I vs total counts, shared genes only, all samples --------

all_moran_individual_NC <- readRDS(moran_individual_nc_file)
all_moran_genes <- readRDS(moran_genes_file)
total_counts_df <- readRDS(total_counts_file)
shared_genes <- readRDS(shared_genes_file)
shared_genes_clean <- make.names(shared_genes)

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
  mutate(
    category = factor(category, levels = c("Pre-designed Genes", "Negative Control Probes", "Negative Control Codeword")),
    sample_factor = factor(sample, levels = bench_samples),
    platform2 = factor(platform2, levels = c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium"))
  ) %>%
  arrange(category)

# 75th percentile of NC total counts per sample x platform (probes + codewords combined)
nc_q75 <- scatter_all %>%
  filter(category %in% c("Negative Control Probes", "Negative Control Codeword")) %>%
  group_by(sample_factor, platform2) %>%
  summarise(q75 = quantile(total_counts, 0.75, na.rm = TRUE), .groups = "drop")

# negative controls merged to one red category; restrict pre-designed genes to
# the 213-gene shared panel
nc_count_colors <- c("Pre-designed Genes" = "black", "Negative Control" = "#e41a1c")

scatter_nc_count <- scatter_all %>%
  mutate(
    color_cat = if_else(category %in% c("Negative Control Probes", "Negative Control Codeword"),
                         "Negative Control", "Pre-designed Genes"),
    color_cat = factor(color_cat, levels = c("Pre-designed Genes", "Negative Control"))
  ) %>%
  arrange(color_cat)

scatter_nc_shared <- scatter_nc_count %>%
  filter(color_cat == "Negative Control" |
         (color_cat == "Pre-designed Genes" & feature_clean %in% shared_genes_clean))

genes_in_nc_range_shared <- scatter_nc_shared %>%
  filter(color_cat == "Pre-designed Genes") %>%
  left_join(nc_q75, by = c("sample_factor", "platform2")) %>%
  filter(!is.na(total_counts), total_counts <= q75)

facet_label_shared <- genes_in_nc_range_shared %>%
  group_by(sample_factor, platform2) %>%
  summarise(n_within = n_distinct(feature_clean), .groups = "drop") %>%
  left_join(
    scatter_nc_shared %>%
      filter(color_cat == "Pre-designed Genes", !is.na(total_counts)) %>%
      group_by(sample_factor, platform2) %>%
      summarise(n_total = n_distinct(feature_clean), .groups = "drop"),
    by = c("sample_factor", "platform2")
  ) %>%
  mutate(pct = round(100 * n_within / n_total, 1), label = paste0(n_within, "\n(", pct, "%)"))

label_pts_shared <- scatter_nc_shared %>%
  distinct(sample_factor, platform2) %>%
  left_join(facet_label_shared %>% select(sample_factor, platform2, label), by = c("sample_factor", "platform2")) %>%
  filter(!is.na(label)) %>%
  mutate(x_pos = 0.9)

shade_df <- nc_q75 %>%
  mutate(xmin = 0.9, xmax = q75, ymin = -Inf, ymax = Inf)

f6b <- ggplot(scatter_nc_shared, aes(x = total_counts, y = lisa, color = color_cat)) +
  geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#ffcccc", alpha = 0.4) +
  geom_vline(data = nc_q75, aes(xintercept = q75), inherit.aes = FALSE,
             color = "#e41a1c", linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_text(data = label_pts_shared, aes(x = x_pos, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.2, size = 2.8, color = "black", fontface = "bold") +
  scale_x_log10(breaks = c(10, 100, 1e3, 1e4, 1e5, 1e6), labels = format_counts) +
  scale_color_manual(values = nc_count_colors, name = "Gene Category") +
  facet_grid(platform2 ~ sample_factor) +
  labs(title = "Moran's I vs. Transcript Counts - shared genes only (n = 213), labeled",
       subtitle = "Label: n shared genes within q75 threshold (percent of shared gene panel)",
       x = "Transcripts per gene", y = "Moran's I (LISA)") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_light_scatter

## ---- save ------------------------------------------------------------------

for (sample in names(f6a_xenium)) {
  ggsave(file.path(out_dir, paste0("f6a_NC_spatial_Xenium_", sample, ".pdf")), f6a_xenium[[sample]], width = 14, height = 4)
}
for (sample in names(f6a_merscope)) {
  ggsave(file.path(out_dir, paste0("f6a_NC_spatial_MERSCOPE_", sample, ".pdf")), f6a_merscope[[sample]], width = 14, height = 4)
}
ggsave(file.path(out_dir, "f6b_moran_vs_counts_shared_genes_all_samples.pdf"), f6b, width = 18, height = 8)
