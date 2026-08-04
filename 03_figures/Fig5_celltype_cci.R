# Purpose:  Fig 5 — cell type and cell-cell interaction analysis.
#             a  UMAP of the snRNA-seq reference plus spatial cell type maps (Xenium,
#                MERSCOPE) for TNBC_02, coloured by cell type
#             b  cell type proportions across all benchmarking samples, sc/snRNA-seq +
#                Xenium + MERSCOPE, after cross-platform alignment
#             d  heatmap: proportion of 110um spatial bins with a significant BLISA
#                ligand-receptor interaction, per pathologist-annotated tissue region x
#                LR pair x platform, TNBC_02
#             e  spatial plot of representative LR interactions (-log10 p), tissue
#                regions bordered and coloured by pathology type, empty bins dropped,
#                TNBC_02
#             f  cell-type-level LR interaction score heatmaps (BLISA-derived CCI),
#                Xenium vs MERSCOPE, TNBC_02, for the same LR pairs as panel e
#           Panel c (pathologist-annotated H&E image) is not code-generated.
#           See README.md.
# Inputs:   <results_table_dir>/comp_df.rds
#                     from ../analysis/cell_type_composition/
#           <results_table_dir>/   from ../analysis/cell_cell_communication/
#             seu_x_annotated.rds, seu_m_annotated.rds, seu_v_annotated.rds
#             annotation_df.rds, prop_all.rds
#             BLISA_x.rds, BLISA_m.rds, BLISA_v.rds, spot_sf.rds, d_theoretical.rds
#             CCI_x.rds, CCI_m.rds
#           <processed_data_dir>/scRNA/TNBC_02/seu_TNBC_02.rds
#                     snRNA reference, for the panel a UMAP
# Outputs:  <figure_output_dir>/
#             f5a_umap_spatial_celltype_TNBC02.pdf
#             f5b_celltype_composition_after_alignment.pdf
#             f5d_region_lr_proportion_heatmap_TNBC02.pdf
#             f5e_lr_spatial_<LR_pair>_TNBC02.pdf         one per representative LR pair
#             f5f_cci_heatmap_<LR_pair>_<platform>_TNBC02.pdf   one per pair x platform

library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(Seurat)
library(sf)
library(ComplexHeatmap)
library(ggnewscale)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))


ccc_dir <- results_table_dir
comp_df_file <- file.path(results_table_dir, "comp_df.rds")
sc_ref_file  <- file.path(processed_data_dir, "scRNA", "TNBC_02", "seu_TNBC_02.rds")  # snRNA reference for TNBC_02
out_dir <- figure_output_dir

sample <- "TNBC_02"  # "TNBC_02" in the manuscript
lr_pairs_plot <- c("CDH1_CDH1", "CCL5_CCR5", "TGFB1_TGFBR1_TGFBR2")  # representative LR pairs, panels e/f

## ---- a: snRNA UMAP + Xenium/MERSCOPE spatial cell type maps --------------

sc_so <- readRDS(sc_ref_file)
p_umap <- DimPlot(sc_so, group.by = "cell_type202605", raster = FALSE, cols = cell_type_colors) +
  ggtitle(paste0(sample, " - snRNA-seq"))

seu_x <- readRDS(file.path(ccc_dir, "seu_x_annotated.rds"))
seu_m <- readRDS(file.path(ccc_dir, "seu_m_annotated.rds"))

ct_spatial_map <- function(seu, pt_size = 0.8) {
  emb <- Embeddings(seu, "aligned_spatial")
  df  <- data.frame(x = emb[, 1], y = emb[, 2], ct = seu$SingleR_labels202605)
  bg  <- intersect(c("unsure", "mixture"), unique(df$ct))
  fg  <- setdiff(unique(df$ct), bg)
  df$ct <- factor(df$ct, levels = c(bg, fg))
  df <- df[order(df$ct), ]

  ggplot(df, aes(x, y, colour = ct)) +
    geom_point(size = pt_size, stroke = 0, shape = 16) +
    scale_colour_manual(values = cell_type_colors, na.value = "grey80",
                        labels = function(v) gsub("\\.", " ", v)) +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"), legend.spacing.y = unit(0.05, "cm"), plot.margin = margin(2, 2, 2, 2)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1))
}

f5a <- cowplot::plot_grid(p_umap, ct_spatial_map(seu_x), ct_spatial_map(seu_m), ncol = 3)

## ---- b: cell type composition, all samples, after alignment --------------

comp_df <- readRDS(comp_df_file)

ct_order <- setdiff(cell_type_order, "mixture")
ct_cols <- rev(cell_type_colors[ct_order])
ct_cols[is.na(ct_cols)] <- "grey70"

comp_df$celltype <- factor(comp_df$celltype, levels = rev(ct_order))

f5b <- ggplot(comp_df, aes(x = platform, y = prop, fill = celltype)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = ct_cols, drop = FALSE, guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~sample, nrow = 1) +
  labs(title = "Cell type composition - after alignment", x = NULL, y = "Proportion", fill = "Cell type") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        strip.background = element_rect(fill = "white"), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))

## ---- pathology colours (shared by d, e) -----------------------------------

annotation_df <- readRDS(file.path(ccc_dir, "annotation_df.rds"))

pathology_colors <- c(
  "Unannotated"                  = "lightgrey",
  "Tumour_more_cohesive"         = "#4A86C8",
  "Lymphocytes"                  = "#918034",
  "Tumour_with_reduced_adhesion" = "#C9579B",
  "Lymphoid_aggregate"           = "#66B832"
)

lighten_col <- function(col, amount = 0.5) { v <- col2rgb(col) / 255; rgb(v[1] + (1 - v[1]) * amount, v[2] + (1 - v[2]) * amount, v[3] + (1 - v[3]) * amount) }
darken_col  <- function(col, amount = 0.4) { v <- col2rgb(col) / 255; rgb(v[1] * (1 - amount), v[2] * (1 - amount), v[3] * (1 - amount)) }

region_colors <- c("Unannotated" = "lightgrey")
for (type in setdiff(names(pathology_colors), "Unannotated")) {
  regs <- sort(unique(annotation_df$pathology_region[annotation_df$pathology_type == type]))
  base <- pathology_colors[[type]]
  shades <- if (length(regs) == 1) base else colorRampPalette(c(lighten_col(base, 0.5), base, darken_col(base, 0.4)))(length(regs))
  region_colors <- c(region_colors, setNames(shades, regs))
}

lri_colors <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")

## ---- d: region x LR-pair proportion-of-significant-bins heatmap ----------

prop_all <- readRDS(file.path(ccc_dir, "prop_all.rds"))

region_type_df <- annotation_df %>%
  select(pathology_region, pathology_type) %>%
  distinct() %>%
  filter(pathology_region != "Unannotated") %>%
  arrange(pathology_type, desc(pathology_region))

region_order  <- region_type_df$pathology_region
type_rle      <- rle(region_type_df$pathology_type)
line_positions <- cumsum(type_rle$lengths)[-length(type_rle$lengths)] + 0.5
label_colors  <- unname(pathology_colors[region_type_df$pathology_type])
plain_labels  <- setNames(gsub("_", " ", region_order), region_order)

f5d <- ggplot(prop_all, aes(x = LR_pair, y = region, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = line_positions, color = "grey30", linewidth = 0.5) +
  facet_wrap(~platform, ncol = 1) +
  scale_y_discrete(limits = region_order, labels = plain_labels) +
  scale_fill_gradientn(colours = c("white", lri_colors), limits = c(0, 1), na.value = "grey80",
                        name = "Proportion\nsig. bins") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7, face = "bold"),
        axis.text.y = element_text(colour = label_colors, face = "bold"),
        strip.text = element_text(face = "bold"), panel.grid = element_blank()) +
  labs(x = "LR pair", y = "Pathology region",
       title = paste(sample, "- proportion of significant bins per region"))

## ---- e: spatial LR interaction significance, pathology-bordered ----------

seu_v   <- readRDS(file.path(ccc_dir, "seu_v_annotated.rds"))
BLISA_v <- readRDS(file.path(ccc_dir, "BLISA_v.rds"))
BLISA_x <- readRDS(file.path(ccc_dir, "BLISA_x.rds"))
BLISA_m <- readRDS(file.path(ccc_dir, "BLISA_m.rds"))
spot_sf <- readRDS(file.path(ccc_dir, "spot_sf.rds"))
d_theoretical <- readRDS(file.path(ccc_dir, "d_theoretical.rds"))

non_empty_x <- rownames(GetTissueCoordinates(seu_v))[spot_sf$Xenium_cell_count_filtered   > 0]
non_empty_m <- rownames(GetTissueCoordinates(seu_v))[spot_sf$MERSCOPE_cell_count_filtered > 0]

# keep only the exterior ring of each polygon part (drop interior holes)
drop_holes <- function(geom) {
  crds  <- as.data.frame(sf::st_coordinates(geom))
  lcols <- grep("^L", names(crds), value = TRUE)
  ext   <- crds[crds[[lcols[1]]] == 1, ]
  if (length(lcols) == 1) {
    sf::st_polygon(list(as.matrix(ext[, c("X", "Y")])))
  } else {
    polys <- lapply(split(ext, ext[[lcols[2]]]), function(s) list(as.matrix(s[, c("X", "Y")])))
    sf::st_multipolygon(unname(polys))
  }
}

# one smooth outline polygon per pathology region: expand each Visium spot,
# union within region, shrink back to smooth the boundary and drop holes
make_region_borders <- function(seu_v, d, region_col = "pathology_region", type_col = "pathology_type",
                                 expand = 0.62, shrink = 0.35, exclude = "Unannotated") {
  coords <- GetTissueCoordinates(seu_v)
  df <- data.frame(x = coords[, "imagecol"], y = coords[, "imagerow"],
                    region = seu_v[[region_col]][rownames(coords), 1], type = seu_v[[type_col]][rownames(coords), 1])
  df  <- df[!is.na(df$region) & !(df$region %in% exclude), ]
  pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = NA)

  parts <- lapply(unique(pts$region), function(rg) {
    sub  <- pts[pts$region == rg, ]
    poly <- sf::st_union(sf::st_buffer(sub, d * expand))
    poly <- sf::st_buffer(poly, -d * shrink)
    poly <- sf::st_sfc(drop_holes(poly))
    if (length(poly) == 0 || all(sf::st_is_empty(poly))) return(NULL)
    sf::st_sf(region = rg, type = unique(sub$type)[1], geometry = poly)
  })
  do.call(rbind, Filter(Negate(is.null), parts))
}

border_coords <- function(border_sf) {
  do.call(rbind, lapply(seq_len(nrow(border_sf)), function(i) {
    cc    <- as.data.frame(sf::st_coordinates(border_sf[i, ]))
    lcols <- grep("^L", names(cc), value = TRUE)
    grp   <- if (length(lcols)) do.call(paste, c(cc[lcols], sep = "-")) else "1"
    data.frame(x = cc$X, y = cc$Y, piece = paste0("f", i, "_", grp), region = border_sf$region[i], type = border_sf$type[i])
  }))
}

plot_pval_filled_borders_nonempty <- function(seu_v, BLISA_output, lr, d,
                                               border_by = "type", shape = "hex", sig_only = TRUE,
                                               pt_size = 1.8, r_factor = 0.5, angle_offset = -90,
                                               line_width = 1, fill_alpha = 0.15, show_border = TRUE,
                                               colors = NULL, non_empty_barcodes = NULL) {
  all_pval <- BLISA_output$LR_out[lr, "all_pval"][[1]]
  sig_idx  <- BLISA_output$LR_out[lr, "sig_index"][[1]]

  neglogp <- rep(NA_real_, length(all_pval))
  if (!sig_only) neglogp <- -log10(all_pval)
  neglogp[sig_idx] <- -log10(all_pval[sig_idx])
  neglogp <- setNames(neglogp, colnames(seu_v))

  coords   <- GetTissueCoordinates(seu_v)
  barcodes <- rownames(coords)
  bc  <- border_coords(make_region_borders(seu_v, d = d))
  pal <- if (!is.null(colors)) colors else if (border_by == "type") pathology_colors else region_colors

  if (shape == "circle") {
    plot_bc <- if (!is.null(non_empty_barcodes)) intersect(barcodes, non_empty_barcodes) else barcodes
    plot_df <- data.frame(x = coords[plot_bc, "imagecol"], y = coords[plot_bc, "imagerow"], neglogp = neglogp[plot_bc])
    p <- ggplot() +
      geom_polygon(data = bc, aes(x, y, group = piece, fill = .data[[border_by]]), alpha = fill_alpha, colour = NA) +
      scale_fill_manual(values = pal, guide = "none") +
      ggnewscale::new_scale_fill() +
      geom_point(data = plot_df, aes(x, y, fill = neglogp), shape = 21, stroke = 0, colour = "transparent", size = pt_size) +
      scale_fill_gradientn(colours = lri_colors, name = "-log10(p)", na.value = "grey90")
  } else {
    r      <- d * r_factor
    angles <- (seq(0, 5) * 60 + angle_offset) * pi / 180
    plot_bc  <- if (!is.null(non_empty_barcodes)) intersect(barcodes, non_empty_barcodes) else barcodes
    plot_pos <- coords[plot_bc, , drop = FALSE]
    hex_df <- data.frame(
      x = rep(plot_pos[, "imagecol"], each = 6) + r * cos(rep(angles, length(plot_bc))),
      y = rep(plot_pos[, "imagerow"], each = 6) + r * sin(rep(angles, length(plot_bc))),
      barcode = rep(plot_bc, each = 6), neglogp = rep(neglogp[plot_bc], each = 6)
    )
    p <- ggplot() +
      geom_polygon(data = hex_df, aes(x = x, y = y, group = barcode, fill = neglogp), color = NA) +
      scale_fill_gradientn(colours = lri_colors, name = "-log10(p)", na.value = "grey90") +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = bc, aes(x, y, group = piece, fill = .data[[border_by]]), alpha = fill_alpha, colour = NA) +
      scale_fill_manual(values = pal, guide = "none")
  }

  if (show_border) {
    p <- p +
      geom_path(data = bc, aes(x, y, group = piece, colour = .data[[border_by]]), linewidth = line_width) +
      scale_colour_manual(values = pal, name = "Pathology", labels = function(v) gsub("_", " ", v))
  }

  p + scale_y_reverse() + coord_fixed() + theme_void() + theme(legend.position = "right") +
    ggtitle(paste0(lr, " - -log10(p) [", shape, "] | ", border_by, " filled borders (non-empty only)"))
}

safe_pval_filled_nonempty <- function(seu_v, BLISA_out, lr, d, platform_name, shape, non_empty = NULL, ...) {
  if (!lr %in% BLISA_out$LR_out$interaction_name)
    return(ggplot() + theme_void() + ggtitle(paste0(platform_name, " - ", lr, " (not detected)")))
  plot_pval_filled_borders_nonempty(seu_v, BLISA_out, lr, d, shape = shape, non_empty_barcodes = non_empty, ...) +
    ggtitle(paste0(platform_name, " - ", lr))
}

f5e_plots <- lapply(lr_pairs_plot, function(lr) {
  wrap_plots(
    safe_pval_filled_nonempty(seu_v, BLISA_v, lr, d_theoretical, "Visium", shape = "circle", border_by = "type", line_width = 0.7, fill_alpha = 0.2),
    safe_pval_filled_nonempty(seu_v, BLISA_x, lr, d_theoretical, "Xenium", shape = "hex", r_factor = 0.42, border_by = "type", line_width = 0.7, fill_alpha = 0.2, non_empty = non_empty_x),
    safe_pval_filled_nonempty(seu_v, BLISA_m, lr, d_theoretical, "MERSCOPE", shape = "hex", r_factor = 0.42, border_by = "type", line_width = 0.7, fill_alpha = 0.2, non_empty = non_empty_m),
    ncol = 3
  )
})
names(f5e_plots) <- lr_pairs_plot

## ---- f: cell-type-level LR interaction score heatmaps (CCI) --------------

CCI_x <- readRDS(file.path(ccc_dir, "CCI_x.rds"))
CCI_m <- readRDS(file.path(ccc_dir, "CCI_m.rds"))

make_cci_heatmap_colored <- function(CCI_df, lr_pair, platform = "") {
  if (!lr_pair %in% colnames(CCI_df)) return(NULL)
  scores    <- CCI_df[[lr_pair]]
  pairs     <- strsplit(rownames(CCI_df), "->")
  receivers <- sapply(pairs, `[`, 2)
  senders   <- sapply(pairs, `[`, 1)
  all_r <- unique(receivers); all_s <- unique(senders)
  mat   <- matrix(NA, nrow = length(all_r), ncol = length(all_s), dimnames = list(all_r, all_s))
  for (i in seq_along(scores)) mat[receivers[i], senders[i]] <- scores[i]

  r_col <- cell_type_colors[rownames(mat)]; r_col[is.na(r_col)] <- "black"
  c_col <- cell_type_colors[colnames(mat)]; c_col[is.na(c_col)] <- "black"

  Heatmap(mat,
    name = "Interaction\nScore", col = viridisLite::viridis(10),
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_title = if (nchar(platform)) paste0(platform, " - Receiver") else "Receiver",
    column_title = if (nchar(platform)) paste0(platform, " - Sender") else "Sender",
    row_names_gp = gpar(fontsize = 12, fontface = "bold", col = r_col),
    column_names_gp = gpar(fontsize = 12, fontface = "bold", col = c_col),
    column_names_rot = 45, heatmap_legend_param = list(title_position = "topcenter")
  )
}

f5f_heatmaps <- lapply(lr_pairs_plot, function(lr) {
  list(Xenium = make_cci_heatmap_colored(CCI_x, lr, "Xenium"), MERSCOPE = make_cci_heatmap_colored(CCI_m, lr, "MERSCOPE"))
})
names(f5f_heatmaps) <- lr_pairs_plot

## ---- save ------------------------------------------------------------------

ggsave(file.path(out_dir, "f5a_umap_spatial_celltype_TNBC02.pdf"), f5a, width = 16, height = 6)
ggsave(file.path(out_dir, "f5b_celltype_composition_after_alignment.pdf"), f5b, width = 14, height = 5)
ggsave(file.path(out_dir, "f5d_region_lr_proportion_heatmap_TNBC02.pdf"), f5d, width = 12, height = 8)
for (lr in lr_pairs_plot) {
  ggsave(file.path(out_dir, paste0("f5e_lr_spatial_", lr, "_TNBC02.pdf")), f5e_plots[[lr]], width = 21, height = 7)
}
for (lr in lr_pairs_plot) {
  for (plat in c("Xenium", "MERSCOPE")) {
    ht <- f5f_heatmaps[[lr]][[plat]]
    if (is.null(ht)) next
    pdf(file.path(out_dir, paste0("f5f_cci_heatmap_", lr, "_", plat, "_TNBC02.pdf")), width = 6, height = 5)
    draw(ht, heatmap_legend_side = "right", padding = unit(c(5, 5, 10, 10), "mm"))
    dev.off()
  }
}
