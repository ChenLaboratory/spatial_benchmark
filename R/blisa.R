# Purpose:  BLISA/CCI method library — cell-cell communication analysis for spatial
#           transcriptomics. Covers ligand-receptor pair filtering, spatial binning to
#           hexagons or Visium spots, the BLISA spatial LR test, cell-type-level CCI
#           scoring, and the associated heatmap/spatial plots. Not run directly:
#           sourced by ../analysis/cell_cell_communication/02_run_blisa.R and
#           03_run_cci.R.
# Inputs:   none (functions only — source into an analysis script)
# Outputs:  none (defines functions in the calling environment)
#
# Sections:
#   1. LR pair filtering
#   2. Expression utilities
#   3. Spatial binning (hex / Visium spots)
#   4. BLISA core
#   5. CCI scoring
#   6. Visualization

library(sf)
library(spdep)
library(Matrix)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(viridisLite)
library(grid)
library(Seurat)


# =============================================================================
# 1. LR pair filtering
# =============================================================================

# Remove an interaction if any subunit of the L or R is absent from gene_list
filter_genes <- function(gene_column, gene_list) {
  sapply(strsplit(as.character(gene_column), ", "), function(genes) {
    filtered <- genes[genes %in% gene_list]
    if (length(filtered) != length(genes)) NA
    else paste(filtered, collapse = ", ")
  })
}

# Return subset of LR_df whose L and R genes are all in gene_panel
getLRpairs <- function(gene_panel, LR_df = CellChatDB.human$interaction) {
  LR_df$ligand.symbol   <- filter_genes(LR_df$ligand.symbol,   gene_panel)
  LR_df$receptor.symbol <- filter_genes(LR_df$receptor.symbol, gene_panel)
  LR_df[!(is.na(LR_df$ligand.symbol) | is.na(LR_df$receptor.symbol)), ]
}

# Keep LR pairs where at least one spot/bin has >= min_ligand / min_receptor counts
filterLRpairs <- function(counts,
                          min_ligand   = 10,
                          min_receptor = 10,
                          LR_df        = CellChatDB.human$interaction) {
  lig_panel <- rownames(counts)[rowSums(counts >= min_ligand)   > 0]
  rec_panel <- rownames(counts)[rowSums(counts >= min_receptor) > 0]

  LR_lig <- getLRpairs(lig_panel, LR_df)
  LR_rec <- getLRpairs(rec_panel, LR_df)

  shared <- intersect(rownames(LR_lig), rownames(LR_rec))
  LR_lig[shared, ]
}


# =============================================================================
# 2. Expression utilities
# =============================================================================

# Split "A, B" or "A_B" into a character vector of gene symbols
parse_units <- function(s) {
  s <- gsub("\\s+", "", as.character(s))
  unlist(strsplit(s, "[,_]"))
}

# Per-bin minimum expression across all subunits of a (possibly multi-unit) gene
# hex_gene_counts: genes x bins matrix
get_min_expr <- function(gene_str, hex_gene_counts) {
  genes <- parse_units(gene_str)
  genes <- genes[genes %in% rownames(hex_gene_counts)]
  n     <- ncol(hex_gene_counts)
  if (length(genes) == 0) return(rep(0, n))  # gene(s) absent -> 0
  if (length(genes) == 1) return(as.numeric(hex_gene_counts[genes, ]))
  # per-bin minimum across subunits
  mat <- hex_gene_counts[genes, , drop = FALSE]
  apply(mat, 2, min)
}

# Assign ccc_mode ("diffuse" / "nearby") based on annotation column
LR_df_add_mode <- function(LR_df,
                           col              = "annotation",
                           default_mode     = "diffuse",
                           diffuse_category = c("Secreted Signaling",
                                                "Non-protein Signaling")) {
  # Case 1: annotation column missing -> all diffuse
  if (!col %in% colnames(LR_df)) {
    LR_df$ccc_mode <- default_mode
    message(col, " column missing - setting ccc_mode='", default_mode, "' for all.")
    return(LR_df)
  }

  # Case 2: annotation exists
  LR_df$ccc_mode <-
    ifelse(
      LR_df[[col]] %in% diffuse_category,
      "diffuse",
      "nearby"
    )

  message("ccc_mode: 'diffuse' for [",
          paste(diffuse_category, collapse = ", "),
          "]; 'nearby' for others.")
  LR_df
}


# =============================================================================
# 3. Spatial binning and Visium aggregation
# =============================================================================

# Bin single cells into hexagonal tiles and aggregate counts
hex_binning_cells <- function(coords_df     = NULL,
                              counts_matrix = NULL,
                              spe           = NULL,
                              hex_size      = 50) {
  ## ---------------------------
  ## 0. Input handling
  ## ---------------------------
  if (!is.null(spe)) {
    coords_df     <- as.data.frame(spe@int_colData$spatialCoords)
    counts_matrix <- spe@assays@data$counts
  } else {
    stopifnot(!is.null(coords_df), !is.null(counts_matrix))
    stopifnot(inherits(counts_matrix, c("matrix", "dgCMatrix")))
  }

  coords_df$cell_id <- rownames(coords_df)
  stopifnot(all(c("cell_id", "x_centroid", "y_centroid") %in% colnames(coords_df)))
  stopifnot(all(colnames(counts_matrix) %in% coords_df$cell_id))

  ## ---------------------------
  ## 1. sf cells
  ## ---------------------------
  cell_sf <- sf::st_as_sf(coords_df,
                          coords = c("x_centroid", "y_centroid"),
                          crs    = NA)
  ## ---------------------------
  ## 2. Hex binning
  ## ---------------------------
  # create bins
  hex_geom <- sf::st_make_grid( # only polygon info
    cell_sf,
    cellsize = hex_size,
    what = "polygons",
    square = FALSE
  )

  hex_sf <- sf::st_sf( # sf object
    hex_id = seq_along(hex_geom),
    geometry = hex_geom
  )
  dim(hex_sf)

  # map cells to bins
  cell_hex_sf <- sf::st_join(cell_sf, hex_sf, join = sf::st_intersects)
  cell_hex_df <- sf::st_drop_geometry(cell_hex_sf)
  dim(cell_hex_df)
  length(unique(cell_hex_df$hex_id))

  ## ---------------------------
  ## 3. Aggregate counts
  ## ---------------------------

  cell_to_hex <- setNames(cell_hex_df$hex_id, cell_hex_df$cell_id)
  cell_to_hex <- cell_to_hex[colnames(counts_matrix)]

  n_hex      <- nrow(hex_sf) # hex_sf contains ALL hexes including empty bins
  hex_factor <- factor(cell_to_hex, levels = seq_len(n_hex))
  H          <- sparse.model.matrix(~ hex_factor - 1)

  hex_gene_counts       <- counts_matrix %*% H
  dim(hex_gene_counts)

  hex_sf$n_cells <- as.numeric(Matrix::colSums(H))
  dim(hex_sf)

  list(
    hex_sf          = hex_sf,
    hex_gene_counts = hex_gene_counts,
    cell_to_hex     = cell_to_hex   # named vector: cell_id -> hex_id
  )

}

# Aggregate sub-spot cells (Xenium / MERSCOPE) to Visium spot resolution
# Returns a genes x spots sparse matrix aligned to spot_ids_keep
aggregate_cells_to_spots <- function(seu_cells,
                                     gene_panel,
                                     spot_ids_keep = NULL,
                                     assay         = "RNA",
                                     slot          = "counts",
                                     spot_col      = "Visium_spot_id") {
  # gene x cell
  counts <- GetAssayData(seu_cells, assay = assay, slot = slot)
  counts <- counts[intersect(rownames(counts), gene_panel), , drop = FALSE]

  # cell -> spot mapping
  spot_id <- seu_cells[[spot_col]][, 1]
  names(spot_id) <- colnames(counts)

  # build sparse cell->spot aggregation matrix
  f <- factor(spot_id, levels = spot_ids_keep %||% unique(spot_id))
  S <- sparse.model.matrix(~ f - 1)  # cells x spots

  # gene x spots
  counts_spot <- counts %*% S
  colnames(counts_spot) <- sub("^f", "", colnames(counts_spot))  # clean colnames

  counts_spot
}

# Theoretical Visium centre-to-centre distance in image pixels
# (100 um physical / 55 um spot diameter) * spot_diameter_fullres * scale_factor
visium_spot_distance <- function(seu, scale = "lowres") {
  sf <- seu@images[[1]]@scale.factors
  spot <- sf$spot
  scale_factor <- sf[[scale]]

  spot_distance <- spot * (100/55) * scale_factor
  return(spot_distance)
}


# =============================================================================
# 4. BLISA core
# =============================================================================

# --- Helper: assign a random neighbour to isolated bins so spdep never gets
#     an empty neighbourhood ---
.fix_isolates <- function(nb, seed = 123) {
  isolates <- which(spdep::card(nb) == 0)
  if (length(isolates)) {
    set.seed(seed)
    for (i in isolates) nb[[i]] <- sample(setdiff(seq_along(nb), i), 1)
    message(length(isolates), " isolated bin(s) given random neighbour.")
  }
  nb
}

# Run BLISA on pre-computed hex bins (spe / hex_sf + hex_gene_counts workflow)
runBLISA <- function(hex_sf,
                     hex_gene_counts,
                     LR_df,
                     spe              = NULL,
                     hex_size         = 100 / 3.63,
                     dmax             = 250 / 3.63,
                     nsim             = 999,
                     p_cutoff         = 0.05,
                     min_ligand       = 10,
                     min_receptor     = 10,
                     col              = "annotation",
                     default_mode     = "diffuse",
                     diffuse_category = c("Secreted Signaling",
                                          "Non-protein Signaling")) {
  centroids <- sf::st_centroid(hex_sf)
  coords    <- sf::st_coordinates(centroids)

  dist_nb  <- .fix_isolates(spdep::dnearneigh(coords, 0, dmax))
  queen_nb <- .fix_isolates(spdep::dnearneigh(centroids, 0, 1.2 * hex_size))

  weight_at_dmax <- 0.01
  dist_wt  <- spdep::nb2listwdist(dist_nb, hex_sf, type = "exp", style = "W",
                                   alpha = -log(weight_at_dmax) / dmax)
  queen_wt <- spdep::nb2listwdist(queen_nb, centroids, type = "idw",
                                   style = "W", zero.policy = TRUE)

  LR_df_filtered <- filterLRpairs(hex_gene_counts, min_ligand, min_receptor, LR_df)
  LR_out         <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_out$sig_numbers <- integer(nrow(LR_out))
  LR_out$sig_index   <- vector("list", nrow(LR_out))
  LR_out$sig_pval    <- vector("list", nrow(LR_out))
  LR_out$all_pval    <- vector("list", nrow(LR_out))
  LR_out$all_lisa    <- vector("list", nrow(LR_out))

  for (i in seq_len(nrow(LR_out))) {
    message(rownames(LR_out)[i])
    mode <- LR_out$ccc_mode[i]
    wt   <- if (mode == "nearby") queen_wt else dist_wt
    message("  mode: ", mode)

    x      <- get_min_expr(LR_out$receptor.symbol[i], hex_gene_counts)
    y      <- get_min_expr(LR_out$ligand.symbol[i],   hex_gene_counts)
    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    hs  <- hotspot(res_bv, Prname = "Pr(folded) Sim", cutoff = p_cutoff,
                   quadrant.type = "pysal", p.adjust = "none")
    idx <- !is.na(hs) & hs == "High-High"

    LR_out$sig_numbers[i] <- sum(idx)
    LR_out$sig_index[[i]] <- which(idx)
    LR_out$sig_pval[[i]]  <- res_bv[idx, "Pr(folded) Sim"]
    LR_out$all_pval[[i]]  <- res_bv[, "Pr(folded) Sim"]
    LR_out$all_lisa[[i]]  <- res_bv[, "Ibvi"]
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]
  front  <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")
  LR_out <- LR_out[, c(front, setdiff(colnames(LR_out), front))]

  list(LR_out = LR_out, hex_sf = hex_sf, hex_gene_counts = hex_gene_counts)
}

# Run BLISA directly on a pre-formed bin/spot count matrix (Visium spots or
# any pre-aggregated matrix). This is the entry point used in CCC_sample_analysis.Rmd.
runBLISA.default <- function(counts_matrix,
                             bin_sf,
                             LR_df,
                             hex_size         = 50,
                             dmax             = 250,
                             nsim             = 999,
                             p_cutoff         = 0.05,
                             min_ligand       = 10,
                             min_receptor     = 10,
                             col              = "annotation",
                             default_mode     = "diffuse",
                             diffuse_category = c("Secreted Signaling",
                                                  "Non-protein Signaling")) {
  centroids <- sf::st_centroid(bin_sf)
  coords    <- sf::st_coordinates(centroids)

  queen_nb <- .fix_isolates(spdep::dnearneigh(coords, 0, 1.2 * hex_size))
  dist_nb  <- .fix_isolates(spdep::dnearneigh(coords, 0, dmax))

  queen_wt <- spdep::nb2listwdist(queen_nb, centroids, type = "idw",
                                   style = "W", zero.policy = TRUE)

  weight_at_dmax <- 0.01
  dist_wt <- spdep::nb2listwdist(dist_nb, bin_sf, type = "exp", style = "W",
                                  alpha = -log(weight_at_dmax) / dmax)

  # Isolate indices (p-val forced to 1, LISA to 0)
  isolate_queen <- which(spdep::card(queen_nb) == 0)
  isolate_dist  <- which(spdep::card(dist_nb)  == 0)

  LR_df_filtered <- filterLRpairs(counts_matrix, min_ligand, min_receptor, LR_df)
  LR_out         <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_out$sig_numbers <- integer(nrow(LR_out))
  LR_out$sig_index   <- vector("list", nrow(LR_out))
  LR_out$sig_pval    <- vector("list", nrow(LR_out))
  LR_out$all_pval    <- vector("list", nrow(LR_out))
  LR_out$all_lisa    <- vector("list", nrow(LR_out))

  for (i in seq_len(nrow(LR_out))) {
    message(rownames(LR_out)[i])
    mode         <- LR_out$ccc_mode[i]
    wt           <- if (mode == "nearby") queen_wt else dist_wt
    isolate_idx  <- if (mode == "nearby") isolate_queen else isolate_dist
    message("  mode: ", mode)

    x      <- get_min_expr(LR_out$receptor.symbol[i], counts_matrix)
    y      <- get_min_expr(LR_out$ligand.symbol[i],   counts_matrix)
    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    # Silence isolated-bin artefacts
    res_bv[isolate_idx, "Pr(folded) Sim"] <- 1
    res_bv[isolate_idx, "Ibvi"]           <- 0

    hs  <- hotspot(res_bv, Prname = "Pr(folded) Sim", cutoff = p_cutoff,
                   quadrant.type = "pysal", p.adjust = "none")
    idx <- !is.na(hs) & hs == "High-High"

    LR_out$sig_numbers[i] <- sum(idx)
    LR_out$sig_index[[i]] <- which(idx)
    LR_out$sig_pval[[i]]  <- res_bv[idx, "Pr(folded) Sim"]
    LR_out$all_pval[[i]]  <- res_bv[, "Pr(folded) Sim"]
    LR_out$all_lisa[[i]]  <- res_bv[, "Ibvi"]
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]
  front  <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")
  LR_out <- LR_out[, c(front, setdiff(colnames(LR_out), front))]

  list(LR_out = LR_out, bin_sf = bin_sf)
}

# Variant of runBLISA.default that truly excludes isolates and low-cell bins
# from the Moran computation (rather than assigning them random neighbours).
# Excluded bins are returned with p=1 and LISA=0 in the full-length output vectors.
# Set n_cells_col = NA (default) to skip cell-count filtering.
runBLISA.default.isolates.removed <- function(
    counts_matrix,
    bin_sf,
    LR_df,
    hex_size          = 50,
    dmax              = 250,
    nsim              = 999,
    p_cutoff          = 0.05,
    min_ligand        = 10,
    min_receptor      = 10,
    min_cells_per_bin = 1,
    n_cells_col       = NA,   # column name in bin_sf for cell counts; NA = skip filtering
    col               = "annotation",
    default_mode      = "diffuse",
    diffuse_category  = c("Secreted Signaling", "Non-protein Signaling")
) {

  centroids <- sf::st_centroid(bin_sf)
  coords    <- sf::st_coordinates(centroids)
  n_bins    <- nrow(bin_sf)

  ## ---------------------------
  ## Filter low-cell bins
  ## ---------------------------
  if (!is.na(n_cells_col)) {
    if (!n_cells_col %in% colnames(bin_sf))
      stop("Column '", n_cells_col, "' not found in bin_sf.")
    low_cell_idx <- which(bin_sf[[n_cells_col]] < min_cells_per_bin)
    message(length(low_cell_idx), " bins removed: < ", min_cells_per_bin,
            " cells (column: '", n_cells_col, "').")
  } else {
    low_cell_idx <- integer(0)
    message("n_cells_col = NA - cell-count filtering skipped.")
  }

  ## ---------------------------
  ## Helper: second-pass isolation check
  ## After subset.nb, bins whose neighbours were all excluded become new isolates.
  ## Detect them and further shrink keep_idx before building weights.
  ## ---------------------------
  resolve_new_isolates <- function(nb, nb_full, keep_idx) {
    new_iso_sub  <- which(spdep::card(nb) == 0)   # positions in subset space
    if (length(new_iso_sub) == 0) return(list(nb = nb, keep_idx = keep_idx))
    new_iso_full <- keep_idx[new_iso_sub]          # back to full-space indices
    message(length(new_iso_full), " bins became isolated after subset - excluded.")
    keep_idx     <- setdiff(keep_idx, new_iso_full)
    nb           <- spdep::subset.nb(nb_full, subset = seq_len(n_bins) %in% keep_idx)
    list(nb = nb, keep_idx = keep_idx)
  }

  ## ---------------------------
  ## Queen spatial weights  (for "nearby" mode)
  ## ---------------------------
  queen_nb_full     <- spdep::dnearneigh(coords, 0, 1.2 * hex_size)
  isolate_idx_queen <- which(spdep::card(queen_nb_full) == 0)
  message(length(isolate_idx_queen), " isolated bins with no nearby neighbours: ",
          paste(isolate_idx_queen, collapse = ","))

  keep_idx_queen <- setdiff(seq_len(n_bins), union(isolate_idx_queen, low_cell_idx))
  queen_nb       <- spdep::subset.nb(queen_nb_full,
                                     subset = seq_len(n_bins) %in% keep_idx_queen)

  r              <- resolve_new_isolates(queen_nb, queen_nb_full, keep_idx_queen)
  queen_nb       <- r$nb;  keep_idx_queen <- r$keep_idx

  queen_wt <- spdep::nb2listwdist(queen_nb, centroids[keep_idx_queen, ],
                                  type = "idw", style = "W", zero.policy = TRUE)

  ## ---------------------------
  ## Distance spatial weights  (for "diffuse" mode)
  ## ---------------------------
  dist_nb_full     <- spdep::dnearneigh(coords, 0, dmax)
  isolate_idx_dist <- which(spdep::card(dist_nb_full) == 0)
  message(length(isolate_idx_dist), " isolated bins with no neighbours within ",
          dmax, " um: ", paste(isolate_idx_dist, collapse = ","))

  keep_idx_dist <- setdiff(seq_len(n_bins), union(isolate_idx_dist, low_cell_idx))
  dist_nb       <- spdep::subset.nb(dist_nb_full,
                                    subset = seq_len(n_bins) %in% keep_idx_dist)

  r             <- resolve_new_isolates(dist_nb, dist_nb_full, keep_idx_dist)
  dist_nb       <- r$nb;  keep_idx_dist <- r$keep_idx

  weight_at_dmax <- 0.01
  dist_wt <- spdep::nb2listwdist(dist_nb, bin_sf[keep_idx_dist, ],
                                 type = "exp", style = "W", zero.policy = TRUE,
                                 alpha = -log(weight_at_dmax) / dmax)

  ## ---------------------------
  ## Filter LR pairs
  ## ---------------------------
  LR_df_filtered <- filterLRpairs(counts_matrix, min_ligand, min_receptor, LR_df)
  LR_out         <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_out$sig_numbers <- integer(nrow(LR_out))
  LR_out$sig_index   <- vector("list", nrow(LR_out))
  LR_out$sig_pval    <- vector("list", nrow(LR_out))
  LR_out$all_pval    <- vector("list", nrow(LR_out))
  LR_out$all_lisa    <- vector("list", nrow(LR_out))

  ## ---------------------------
  ## Local bivariate Moran
  ## ---------------------------
  for (i in seq_len(nrow(LR_out))) {
    message(rownames(LR_out)[i])
    mode <- LR_out$ccc_mode[i]

    if (mode == "nearby") {
      wt       <- queen_wt
      keep_idx <- keep_idx_queen
    } else {
      wt       <- dist_wt
      keep_idx <- keep_idx_dist
    }
    message("  ccc mode: ", mode)

    x <- get_min_expr(LR_out$receptor.symbol[i], counts_matrix)[keep_idx]
    y <- get_min_expr(LR_out$ligand.symbol[i],   counts_matrix)[keep_idx]

    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    # Restore to full-length vectors; excluded/isolated bins default to p=1, LISA=0
    full_pval <- rep(1, n_bins)
    full_lisa <- rep(0, n_bins)
    full_pval[keep_idx] <- res_bv[, "Pr(folded) Sim"]
    full_lisa[keep_idx] <- res_bv[, "Ibvi"]

    hs     <- spdep::hotspot(res_bv, Prname = "Pr(folded) Sim", cutoff = p_cutoff,
                             quadrant.type = "pysal", p.adjust = "none")
    idx_hh <- !is.na(hs) & hs == "High-High"
    HH_idx <- keep_idx[which(idx_hh)]   # convert subset-space to full-space indices

    LR_out$sig_numbers[i] <- length(HH_idx)
    LR_out$sig_index[[i]] <- HH_idx
    LR_out$sig_pval[[i]]  <- full_pval[HH_idx]
    LR_out$all_pval[[i]]  <- full_pval
    LR_out$all_lisa[[i]]  <- full_lisa
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]
  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")
  LR_out <- LR_out[, c(front_cols, setdiff(colnames(LR_out), front_cols))]

  list(
    LR_out            = LR_out,
    bin_sf            = bin_sf,
    queen_nb_full     = queen_nb_full,   # full nb (all bins) - for CCI niche lookup
    dist_nb_full      = dist_nb_full,
    keep_idx_queen    = keep_idx_queen,  # bins used in Moran (queen mode)
    keep_idx_dist     = keep_idx_dist,
    isolate_idx_queen = isolate_idx_queen,
    isolate_idx_dist  = isolate_idx_dist
  )
}


# =============================================================================
# 5. CCI scoring
# =============================================================================

# Map BLISA hotspot spots to single cells, score sender/receiver by cell type.
# Returns a wide data frame: rows = CellPair (Sender->Receiver),
#                            cols = LR_pair, values = geometric-mean score.
runCCI_SpotAligned_MatchingWeights <- function(
    seu,
    BLISA_res,
    spot_ids,
    spot_col  = "Visium_spot_id",
    ct_group  = "SingleR_labels"
) {
  LRI_sum <- BLISA_res$LR_out

  # Neighbour lists and keep indices carried over from BLISA -- no recomputation
  queen_nb_full  <- BLISA_res$queen_nb_full
  dist_nb_full   <- BLISA_res$dist_nb_full
  keep_idx_queen <- BLISA_res$keep_idx_queen  # NULL when using runBLISA.default
  keep_idx_dist  <- BLISA_res$keep_idx_dist

  cell_to_spot <- seu[[spot_col]][, 1]
  names(cell_to_spot) <- colnames(seu)
  all_cts <- unique(as.character(seu[[ct_group]][, 1]))

  get_ct_sums <- function(target_cells, gene_str, seu_obj, ct_col) {
    if (length(target_cells) == 0) return(numeric(0))
    genes      <- unlist(strsplit(gsub("\\s+", "", gene_str), "[,_]"))
    raw        <- GetAssayData(seu_obj, slot = "counts")[
      intersect(genes, rownames(seu_obj)), target_cells, drop = FALSE]
    counts_vec <- if (nrow(raw) > 1) apply(raw, 2, min) else as.numeric(raw)
    tapply(counts_vec, seu_obj[[ct_col]][target_cells, 1], sum)
  }

  interaction_list <- lapply(seq_len(nrow(LRI_sum)), function(idx) {
    if (LRI_sum$sig_numbers[idx] == 0) return(NULL)

    mode      <- LRI_sum$ccc_mode[idx]
    sigHH_idx <- LRI_sum$sig_index[[idx]]
    sigHH_id  <- spot_ids[sigHH_idx]

    nb_full  <- if (mode == "nearby") queen_nb_full else dist_nb_full
    keep_idx <- if (mode == "nearby") keep_idx_queen else keep_idx_dist
    # Niche: hotspot bins + neighbours; intersect with keep_idx to drop excluded bins
    ng_all <- unique(c(sigHH_idx, unlist(nb_full[sigHH_idx])))
    ng_idx <- if (!is.null(keep_idx)) intersect(ng_all, keep_idx) else ng_all
    ng_id  <- spot_ids[ng_idx]

    cells_rec  <- names(cell_to_spot)[cell_to_spot %in% sigHH_id]
    cells_send <- names(cell_to_spot)[cell_to_spot %in% ng_id]

    rec_sums  <- get_ct_sums(cells_rec,  LRI_sum$receptor.symbol[idx], seu, ct_group)
    send_sums <- get_ct_sums(cells_send, LRI_sum$ligand.symbol[idx],   seu, ct_group)

    res_df <- data.frame(
      receiver_sum = as.numeric(rec_sums[all_cts]),
      sender_sum   = as.numeric(send_sums[all_cts]),
      row.names    = all_cts
    )
    res_df[is.na(res_df)] <- 0

    score_mat           <- sqrt(outer(log10(res_df$receiver_sum + 1),
                                      log10(res_df$sender_sum   + 1),
                                      FUN = "*"))
    dimnames(score_mat) <- list(all_cts, all_cts)

    df           <- as.data.frame(as.table(score_mat))
    colnames(df) <- c("Receiver", "Sender", "Score")
    df$LR_pair   <- rownames(LRI_sum)[idx]
    df$CellPair  <- paste(df$Sender, df$Receiver, sep = "->")
    df
  })

  all_df <- do.call(rbind, interaction_list)
  wide   <- reshape2::dcast(all_df, CellPair ~ LR_pair, value.var = "Score", fill = 0)
  rownames(wide) <- wide$CellPair
  wide[, -1, drop = FALSE]
}


# =============================================================================
# 6. Visualization
# =============================================================================

# Dot plot: LR pairs ranked by number of significant hotspot spots
plotLRIsum <- function(LR_out, top = NULL, pt.size = 4) {
  if (!is.null(top)) LR_out <- LR_out[seq_len(min(top, nrow(LR_out))), ]

  LR_out$LR_pair <- rownames(LR_out)
  LR_out <- LR_out %>%
    dplyr::arrange(dplyr::desc(sig_numbers)) %>%
    dplyr::mutate(LR_pair = factor(LR_pair, levels = rev(LR_pair))) # reverse so high at top

  ggplot(LR_out, aes(x = sig_numbers, y = LR_pair, color = annotation)) +
    geom_point(size = pt.size) +
    scale_color_manual(values = c(
      "Secreted Signaling"     = "#90a955",
      "ECM-Receptor"           = "#219ebc",
      "Cell-Cell Contact"      = "#f7b801",
      "Non-protein Signaling"  = "#9f86c0"
    )) +
    scale_x_continuous(expand = expansion(add = 100)) +
    labs(x = "Sig Spot Numbers", y = "Ligand-Receptor Pair", color = "Annotation") +
    theme_minimal() +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 10),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text.y      = element_text(face = "bold", size = 12),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    ) +
    coord_cartesian(clip = "off")
}

# Overlay BLISA -log10(p-val) on Visium tissue image for one LR pair
plot_LRI <- function(seu_v, result_df, LR) {
  res      <- result_df[LR, ]
  sig_idx  <- res$sig_index[[1]]
  sig_pval <- res$sig_pval[[1]]

  pvals <- rep(NA, nrow(seu_v))
  pvals[sig_idx] <- sig_pval
  seu_v$BLISA_pval <- -log10(pvals)

  lri_colors <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")
  SpatialFeaturePlot(seu_v, features = "BLISA_pval",
                     pt.size.factor = 3, crop = TRUE) +
    scale_fill_gradientn(colours = lri_colors) +
    ggtitle(LR)
}

# Heatmap of CCI scores: rows = cell-type pairs, cols = LR pairs
# cell_type_colors comes from config.R, sourced by the calling script
CCIheatmap <- function(CCI_df, include_celltypes = NULL) {
  if (!is.null(include_celltypes)) {
    pairs     <- strsplit(rownames(CCI_df), "->")
    senders   <- sapply(pairs, `[`, 1)
    receivers <- sapply(pairs, `[`, 2)
    keep      <- senders %in% include_celltypes | receivers %in% include_celltypes
    CCI_df    <- CCI_df[keep, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No matching cell-type pairs for include_celltypes.")
  }

  pairs     <- strsplit(rownames(CCI_df), "->")
  senders   <- sapply(pairs, `[`, 1)
  receivers <- sapply(pairs, `[`, 2)

  row_ha <- rowAnnotation(
    Sender   = senders,
    Receiver = receivers,
    col      = list(Sender = cell_type_colors, Receiver = cell_type_colors),
    annotation_legend_param = list(Sender   = list(show = FALSE),
                                   Receiver = list(show = FALSE)),
    annotation_name_gp   = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "top"
  )

  ht <- Heatmap(
    as.matrix(CCI_df),
    name                  = "Interaction\nScore",
    col                   = viridisLite::viridis(10),
    cluster_rows          = TRUE,
    cluster_columns       = TRUE,
    row_names_gp          = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp       = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot      = 45,
    heatmap_legend_param  = list(title_position = "topcenter"),
    left_annotation       = row_ha
  )

  draw(ht, heatmap_legend_side = "right",
       padding = unit(c(5, 15, 5, 10), "mm"))
  invisible(ht)
}

# Qualitative colour palette (30 colours, distinct-first ordering)
cols <- c(
  "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7",
  "#117a77", "#F0E442", "#0b2b5e", "#33A02C", "#E31A1C",
  "#6A3D9A", "#B15928", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#CAB2D6", "#FB9A99", "#a18e6a", "#FF7F00",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FCCDE5", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#A6CEE3", "#B2DF8A", "#FFB3BA"
)

# Fallback: map cells in a SPE to hex bin IDs via spatial join.
# Prefer passing hex_binning_cells()$cell_to_hex directly to CCIspatial.
get_cell_hex_mapping <- function(spe, hex_sf) {
  coords_df <- as.data.frame(spe@int_colData$spatialCoords)
  cell_sf   <- sf::st_as_sf(coords_df,
                             coords = c("x_centroid", "y_centroid"), crs = NA)
  joined    <- sf::st_drop_geometry(
    sf::st_join(cell_sf, hex_sf["hex_id"], join = sf::st_intersects)
  )
  setNames(joined$hex_id, rownames(joined))
}

# Spatial map of the dominant sender->receiver cell-type pair at each hotspot hex bin.
# BLISA_output must come from runBLISA (hex workflow), which returns hex_sf,
# queen_nb_full, and dist_nb_full.
# Pass cell_to_hex = hex_binning_cells()$cell_to_hex to avoid recomputation.
CCIspatial <- function(
    spe,
    BLISA_output,
    index,
    cell_to_hex = NULL,  # from hex_binning_cells()$cell_to_hex; recomputed if NULL
    ct_group    = "cell_type",
    top         = 30,
    hex_size    = 50,    # fallback only, if nb not stored in BLISA_output
    dmax        = 250
) {
  LRI_sum <- BLISA_output$LR_out
  hex_sf  <- BLISA_output$hex_sf

  interaction <- unname(unlist(LRI_sum[index, c("ligand.symbol", "receptor.symbol")]))
  sigHH       <- LRI_sum$sig_index[[index]]
  mode        <- LRI_sum$ccc_mode[index]

  # Use pre-computed nb from BLISA if available, else recompute
  if (mode == "nearby" && !is.null(BLISA_output$queen_nb_full)) {
    nb_list <- BLISA_output$queen_nb_full
  } else if (mode != "nearby" && !is.null(BLISA_output$dist_nb_full)) {
    nb_list <- BLISA_output$dist_nb_full
  } else {
    centroids <- sf::st_centroid(hex_sf)
    coords    <- sf::st_coordinates(centroids)
    nb_list   <- if (mode == "nearby") spdep::dnearneigh(centroids, 0, 1.2 * hex_size)
                 else                  spdep::dnearneigh(coords, 0, dmax)
  }

  # Use pre-computed cell-to-hex mapping if provided, else recompute via spatial join
  if (is.null(cell_to_hex)) cell_to_hex <- get_cell_hex_mapping(spe, hex_sf)

  cell_data <- data.table::data.table(
    hex_id        = as.integer(cell_to_hex),
    ct            = as.character(SummarizedExperiment::colData(spe)[[ct_group]]),
    ligand_expr   = as.numeric(counts(spe)[interaction[1], ]),
    receptor_expr = as.numeric(counts(spe)[interaction[2], ])
  )

  # Receiver cells inside hotspot bins; sender cells in hotspot neighbourhood
  rcpt_summary <- cell_data[hex_id %in% sigHH,
                             .(r_sum = sum(receptor_expr)), by = .(hex_id, ct_r = ct)]

  hh_nb_map <- data.table::data.table(
    hh_hex = rep(sigHH, sapply(nb_list[sigHH], length)),
    nb_hex  = unlist(nb_list[sigHH])
  )

  lig_summary <- merge(hh_nb_map, cell_data,
                       by.x = "nb_hex", by.y = "hex_id", allow.cartesian = TRUE)
  lig_summary <- lig_summary[, .(l_sum = sum(ligand_expr)), by = .(hh_hex, ct_l = ct)]

  merged_scores <- merge(rcpt_summary, lig_summary,
                         by.x = "hex_id", by.y = "hh_hex", allow.cartesian = TRUE)
  merged_scores[, product   := (log10(r_sum + 1) + log10(l_sum + 1)) / 2]
  merged_scores[, cell_pair := paste(ct_l, ct_r, sep = " -> ")]

  # Dominant pair per hotspot bin
  top_pairs <- merged_scores[merged_scores[, .I[which.max(product)], by = hex_id]$V1]

  # Legend: top N pairs; remainder labelled "rare pairs"
  tbl            <- sort(table(top_pairs$cell_pair), decreasing = TRUE)
  filtered_pairs <- names(tbl[seq_len(min(top, length(tbl)))])
  legend_title   <- if (length(tbl) <= top) "All pairs" else paste0("Top ", top, " pairs")

  top_pairs[, cell_pair_plot := ifelse(cell_pair %in% filtered_pairs,
                                       cell_pair, "rare pairs")]

  # Assign category to every hex bin
  hex_sf$cell_pair_plot <- ifelse(hex_sf$n_cells > 0, "Non-Significant", "Empty")
  hex_sf$cell_pair_plot[top_pairs$hex_id] <- top_pairs$cell_pair_plot

  legend_levels         <- c(filtered_pairs, "rare pairs", "Non-Significant", "Empty")
  hex_sf$cell_pair_plot <- factor(hex_sf$cell_pair_plot, levels = legend_levels)

  ggplot(hex_sf) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(values = c(
      "Empty"           = "#F0F0F0",
      "Non-Significant" = "#D3D3D3",
      "rare pairs"      = "#818589",
      setNames(cols[seq_along(filtered_pairs)], filtered_pairs)
    )) +
    theme_void() +
    guides(fill = guide_legend(title = legend_title)) +
    labs(
      title    = paste0(rownames(LRI_sum)[index], ": Interacting Hotspots"),
      subtitle = "Grey: light = empty bins, medium = non-significant bins"
    )
}

# Heatmap of CCI scores for a single LR pair (receiver x sender matrix)
CCIheatmapOneLR <- function(CCI_df, lr_pair) {
  if (!lr_pair %in% colnames(CCI_df))
    stop("LR pair '", lr_pair, "' not found in CCI data frame.")

  scores    <- CCI_df[[lr_pair]]
  pairs     <- strsplit(rownames(CCI_df), "->")
  receivers <- sapply(pairs, `[`, 2)
  senders   <- sapply(pairs, `[`, 1)

  all_r <- unique(receivers)
  all_s <- unique(senders)
  mat   <- matrix(NA, nrow = length(all_r), ncol = length(all_s),
                  dimnames = list(all_r, all_s))
  for (i in seq_along(scores)) mat[receivers[i], senders[i]] <- scores[i]

  Heatmap(
    mat,
    name                 = "Interaction\nScore",
    col                  = viridisLite::viridis(10),
    cluster_rows         = TRUE,
    cluster_columns      = TRUE,
    row_title            = paste0(lr_pair, " - Receiver"),
    column_title         = paste0(lr_pair, " - Sender"),
    row_names_gp         = gpar(fontsize = 12, fontface = "bold"),
    column_names_gp      = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot     = 45,
    heatmap_legend_param = list(title_position = "topcenter")
  )
}
