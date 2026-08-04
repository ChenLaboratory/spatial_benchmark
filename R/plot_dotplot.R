# Purpose:  Marker gene reference and dot plot helpers used when annotating cell types.
#           Loads the marker table, appends markers it does not cover, defines the cell
#           type ranking, and provides plot_marker_dotplot() (annotation-time marker
#           inspection) and plot_manuscript_dotplot() (a fixed gene panel). Not run
#           directly: sourced by ../preprocessing/scRNA/03_annotate_cell_types.R.
# Inputs:   HsMarkers.txt   marker gene table, read at source time — its path is
#                           hardcoded below and must be edited before sourcing
# Outputs:  none (defines objects and functions in the calling environment)

# ============================================================================
# Load marker genes and define cell type ranking
# ============================================================================

# Load markers from HsMarkers.txt
markers_df <- read.table(
  "/vast/projects/Spatial/lei/Benchmarking/new_cell_type_annotation_202605/HsMarkers.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Append markers not in HsMarkers.txt
markers_df <- rbind(
  markers_df,
  data.frame(
    Cell_Type = c("LP",   "LP",    "ML",   "ML",    "Epithelial", "Epithelial"),
    Genes     = c("ELF5", "MMP7",  "PRLR", "AREG",  "EPCAM",      "ITGA6"),
    stringsAsFactors = FALSE
  )
)

# Create markers_list: cell type -> unique genes
markers_list <- markers_df |>
  dplyr::group_by(Cell_Type) |>
  dplyr::summarise(genes = list(unique(Genes)), .groups = "drop") |>
  tibble::deframe()

# Biological parent groups — controls cell type ordering in the dot plot
ct_parent_groups <- list(
  Tumor        = c("Tumour"),
  Myeloid      = c("Monocyte", "Macrophage", "TAM", "cDC", "pDC", "Mast", "Neutrophil"),
  T            = c("CD4_T", "CD8_T", "Treg", "Ex_T", "NK"),
  Endothelial  = c("Endothelial"),
  Pericyte     = c("Pericyte"),
  Fibroblast   = c("Fibroblast", "iCAF", "myCAF", "apCAF"),
  Adipocyte    = c("Adipocyte"),
  B            = c("Naive_B", "Memory_B"),
  Plasma       = c("Plasma"),
  Epithelial   = c("Basal", "LP", "ML", "Luminal", "Epithelial"),
  Lymph.vessel = c("Lymph.vessel")
)

# Flat rank: name -> position across all groups (preserves within-group order)
ct_rank <- setNames(seq_along(unlist(ct_parent_groups)), unlist(ct_parent_groups))

# ============================================================================
# Compute dot plot data without calling Seurat::DotPlot() (which uses the
# defunct formula interface of facet_grid() in recent ggplot2 versions).
.dot_data <- function(seu, features, group_by = NULL) {
  if (is.null(group_by)) {
    groups       <- as.character(Seurat::Idents(seu))
    group_levels <- levels(Seurat::Idents(seu))
  } else {
    col          <- seu@meta.data[[group_by]]
    groups       <- as.character(col)
    group_levels <- if (is.factor(col)) levels(col) else sort(unique(groups))
  }

  # Compatible with Seurat v4 (slot=) and v5 (layer=)
  data_use <- tryCatch(
    Seurat::GetAssayData(seu, layer = "data"),
    error = function(e) Seurat::GetAssayData(seu, slot = "data")
  )
  features  <- features[features %in% rownames(data_use)]
  data_use  <- as.matrix(data_use[features, , drop = FALSE])

  result <- do.call(rbind, lapply(group_levels, function(grp) {
    mat <- data_use[, groups == grp, drop = FALSE]
    data.frame(
      features.plot    = features,
      id               = grp,
      pct.exp          = rowMeans(mat > 0) * 100,
      avg.exp          = rowMeans(mat),
      stringsAsFactors = FALSE
    )
  }))

  result <- result |>
    dplyr::group_by(features.plot) |>
    dplyr::mutate(avg.exp.scaled = {
      s <- scale(avg.exp)[, 1]
      pmin(pmax(s, -2.5), 2.5)
    }) |>
    dplyr::ungroup() |>
    as.data.frame()

  result$id <- factor(result$id, levels = group_levels)
  result
}

# ============================================================================
#' Plot Marker Dotplot
#'
#' Creates a dotplot showing marker gene expression across cell types/clusters.
#' Uses embedded markers_list and ct_rank defined at the file level.
#'
#' @param seu A Seurat object with expression data and cluster/identity information
#' @param sample_id A character string representing the sample or group name for the plot title
#'
#' @return A ggplot object
#'
#' @details
#' The function:
#' - Filters markers to only include genes present in the Seurat object
#' - Sorts genes alphabetically within each cell type
#' - Orders cell types according to ct_rank (defined at file level)
#' - Ensures each gene appears only once (in its first cell type)
#' - Generates a dotplot with color (avg expression) and size (percent expressing)
#'
#' @examples
#' \dontrun{
#' p <- plot_marker_dotplot(seu, "Sample_001")
#' print(p)
#' }
#'
#' @importFrom Seurat DotPlot
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient scale_size
#'   facet_grid theme_bw theme element_text margin labs ggtitle
#'
#' @export
plot_marker_dotplot <- function(seu, sample_id, group_by = NULL) {

  # Validate inputs
  if (!inherits(seu, "Seurat")) {
    stop("seu must be a Seurat object")
  }
  if (!is.character(sample_id) || length(sample_id) != 1) {
    stop("sample_id must be a single character string")
  }
  if (!is.null(group_by) && !group_by %in% colnames(seu@meta.data)) {
    stop("group_by column '", group_by, "' not found in seu@meta.data")
  }
  
  # Process markers: keep only genes in seu, sort alphabetically within each cell type
  markers_plot <- lapply(markers_list, function(g) {
    sort(g[g %in% rownames(seu)])
  })
  
  # Remove cell types with no genes
  markers_plot <- markers_plot[lengths(markers_plot) > 0]
  
  # Sort by biological parent group order (ct_rank); 
  # unrecognized cell types go last
  rank_vals <- ct_rank[names(markers_plot)]
  rank_vals[is.na(rank_vals)] <- max(ct_rank) + seq_len(sum(is.na(rank_vals)))
  markers_plot <- markers_plot[order(rank_vals)]
  
  # Ensure unique genes across groups: keep each gene in its first cell type only
  seen <- character(0)
  markers_plot <- lapply(markers_plot, function(g) {
    out <- g[!g %in% seen]
    seen <<- c(seen, out)
    out
  })
  markers_plot <- markers_plot[lengths(markers_plot) > 0]
  
  # Flatten gene list and create mapping to cell types
  genes_flat <- unlist(markers_plot)
  gene_groups <- rep(names(markers_plot), lengths(markers_plot))
  
  # Compute pct.exp and avg.exp.scaled without calling Seurat::DotPlot()
  dp_data <- .dot_data(seu, genes_flat, group_by)
  
  # Add cell type information to dotplot data
  dp_data$cell_type <- factor(
    gene_groups[match(as.character(dp_data$features.plot), genes_flat)],
    levels = names(markers_plot)
  )
  
  # Reverse gene order so first cell type is at top
  dp_data$features.plot <- factor(dp_data$features.plot,
                                   levels = rev(genes_flat))
  
  # Create the dotplot
  p <- ggplot2::ggplot(dp_data, ggplot2::aes(x = id, y = features.plot,
                                              color = avg.exp.scaled, 
                                              size = pct.exp)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "lightgrey", high = "blue", 
                                   name = "Avg. Exp.\n(scaled)") +
    ggplot2::scale_size(range = c(0, 4), name = "Pct. Exp.") +
    ggplot2::facet_grid(rows = ggplot2::vars(cell_type), scales = "free_y", space = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0(sample_id, " — Marker Dot Plot")) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, 
                                           size = 12, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7),
      strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, size = 9),
      panel.spacing = ggplot2::unit(0.1, "lines"),
      plot.margin = ggplot2::margin(t = 10, r = 80, b = 10, l = 10, 
                                     unit = "pt")
    ) +
    ggplot2::labs(x = NULL, y = NULL)
  
  return(p)
}

# ============================================================================
#' Manuscript-style dotplot restricted to a shared gene panel
#'
#' @param seu         Seurat object with final cell type annotations
#' @param sample_id   Sample label for plot title
#' @param panel_genes Character vector of panel genes to restrict to
#' @param group_by    Metadata column with cell type labels (default "cell_type202605")
#'
#' @return A ggplot object
#' @export
plot_manuscript_dotplot <- function(seu, sample_id, panel_genes,
                                    group_by = "cell_type202605") {
  if (!inherits(seu, "Seurat"))
    stop("seu must be a Seurat object")

  # Collect genes per parent group (ct_parent_groups order preserved)
  markers_panel <- lapply(ct_parent_groups, function(sub_cts) {
    genes <- markers_df$Genes[markers_df$Cell_Type %in% sub_cts]
    sort(unique(genes[genes %in% panel_genes & genes %in% rownames(seu)]))
  })
  markers_panel <- markers_panel[lengths(markers_panel) > 0]

  # Deduplicate genes across parent groups (keep in first group)
  seen <- character(0)
  markers_panel <- lapply(markers_panel, function(g) {
    out <- g[!g %in% seen]; seen <<- c(seen, out); out
  })
  markers_panel <- markers_panel[lengths(markers_panel) > 0]

  if (length(markers_panel) == 0) {
    message("No panel genes found in ", sample_id, " — skipping dotplot")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  genes_flat  <- unlist(markers_panel)
  gene_groups <- rep(names(markers_panel), lengths(markers_panel))

  ct_present <- unique(as.character(seu@meta.data[[group_by]]))
  ct_order   <- c(
    names(ct_parent_groups)[names(ct_parent_groups) %in% ct_present],
    setdiff(ct_present, names(ct_parent_groups))
  )

  seu_tmp <- seu
  seu_tmp@meta.data[[group_by]] <- factor(seu_tmp@meta.data[[group_by]], levels = ct_order)

  dp_data <- .dot_data(seu_tmp, genes_flat, group_by)
  dp_data$cell_type <- factor(
    gene_groups[match(as.character(dp_data$features.plot), genes_flat)],
    levels = names(markers_panel)
  )
  dp_data$features.plot <- factor(dp_data$features.plot, levels = rev(genes_flat))

  ggplot2::ggplot(dp_data, ggplot2::aes(x = id, y = features.plot,
                                         color = avg.exp.scaled, size = pct.exp)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "lightgrey", high = "blue",
                                   name = "Avg. Exp.\n(scaled)") +
    ggplot2::scale_size(range = c(0, 4), name = "Pct. Exp.") +
    ggplot2::facet_grid(rows = ggplot2::vars(cell_type),
                         scales = "free_y", space = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y      = ggplot2::element_text(size = 7),
      strip.text.y     = ggplot2::element_text(angle = 0, hjust = 0, size = 9),
      panel.spacing    = ggplot2::unit(0.1, "lines"),
      plot.margin      = ggplot2::margin(t = 10, r = 80, b = 10, l = 10, unit = "pt")
    ) +
    ggplot2::labs(x = NULL, y = NULL)
}