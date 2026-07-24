# Per-(sample, platform, segmentation-method) QC metrics: total/assigned
# transcripts, total cells, median cell area, median transcripts per cell --
# the metrics used by Fig4b (total cell count, median cell area, transcript
# assignment rate, transcripts per cell). See README.md.

library(Seurat)
library(tidyverse)
library(data.table)
library(arrow)  # read_parquet
library(sf)     # st_as_sfc, st_area

control_pattern <- "(?i)^(Blank|FC\\d+$|Sabrina|NegControl)"

samples_xenium   <- c("MH0007", "MH0026", "ER0114", "ER0360", "TN0177", "TN0554")
samples_merscope <- c("MH0007", "MH0026", "ER0114", "ER0360", "TN0177")

# Xenium default so_raw.rds, per sample (preprocessing output)
xenium_default_paths <- setNames(
  file.path("/path/to/preprocessing/Xenium", samples_xenium, "so_raw.rds"), samples_xenium)
# MERSCOPE default so_raw.rds (unfiltered, volume not yet replaced by 2D area), per sample
merscope_default_paths <- setNames(
  file.path("/path/to/preprocessing/MERSCOPE", samples_merscope, "so_raw.rds"), samples_merscope)
# MERSCOPE cell_boundaries.parquet, used to compute 2D XY cell area (MERSCOPE's
# default output stores 3D volume instead)
merscope_cell_boundaries_paths <- setNames(
  file.path("/path/to/raw/MERSCOPE", samples_merscope, "cell_boundaries.parquet"), samples_merscope)

# raw transcript files, for total (pre-segmentation) transcript counts --
# shared between Default and Proseg (same raw file, QC/control-filtered)
xenium_tx_paths <- setNames(
  file.path("/path/to/raw/Xenium", samples_xenium, "transcripts.csv.gz"), samples_xenium)
merscope_tx_paths <- setNames(
  file.path("/path/to/raw/MERSCOPE", samples_merscope, "detected_transcripts.csv"), samples_merscope)

baysor_xenium_dir   <- "/path/to/run_cell_segmentation/baysor/results/xenium"    # <sample>_so_raw.rds
baysor_merscope_dir <- "/path/to/run_cell_segmentation/baysor/results/merscope"  # <sample>_so_raw.rds
proseg_xenium_dir    <- "/path/to/run_cell_segmentation/proseg_v3/results/xenium"    # <sample>_so_raw.rds, <sample>_proseg_xy_area.csv
proseg_merscope_dir  <- "/path/to/run_cell_segmentation/proseg_v3/results/merscope"  # <sample>_so_raw.rds, <sample>_proseg_xy_area.csv

# ---- cell-area helpers -----------------------------------------------------

# MERSCOPE default stores 3D volume; replace with the 2D XY projected area at
# the median z-slice, computed from the cell boundary polygons
load_merscope_seu_2d_area <- function(sample) {
  bounds <- read_parquet(merscope_cell_boundaries_paths[[sample]])
  bounds_sf <- st_sf(
    EntityID = as.character(bounds$EntityID),
    ZIndex   = bounds$ZIndex,
    geometry = st_as_sfc(bounds$Geometry, crs = NA)
  )
  mid_z <- median(unique(bounds_sf$ZIndex))
  xy_per_cell <- bounds_sf %>%
    filter(ZIndex == mid_z) %>%
    mutate(xy_area_um2 = as.numeric(st_area(geometry))) %>%
    st_drop_geometry() %>%
    select(EntityID, xy_area_um2)
  so <- readRDS(merscope_default_paths[[sample]])
  so$cell_area <- xy_per_cell$xy_area_um2[match(colnames(so), xy_per_cell$EntityID)]
  so
}

# Proseg stores 3D volume/surface_area by default; attach the 2D XY area
# (computed from Proseg's boundary polygons, see run_cell_segmentation/proseg_v3/03_run_cell_area_2d.sh)
load_proseg_seu_2d_area <- function(sample, platform = c("xenium", "merscope")) {
  platform <- match.arg(platform)
  dir <- if (platform == "xenium") proseg_xenium_dir else proseg_merscope_dir
  so <- readRDS(file.path(dir, paste0(sample, "_so_raw.rds")))
  xy_per_cell <- read.csv(file.path(dir, paste0(sample, "_proseg_xy_area.csv")))
  xy_per_cell$cell <- as.character(xy_per_cell$cell)
  so$cell_area <- xy_per_cell$xy_area_um2[match(colnames(so), xy_per_cell$cell)]
  so
}

# ---- core metrics row builder ----------------------------------------------

compute_qc_row <- function(so, sample, platform, segmentation, total_transcripts,
                            area_col = "cell_area", counts_layer = "counts") {
  counts  <- GetAssayData(so, assay = "RNA", layer = counts_layer)
  is_ctrl <- grepl(control_pattern, rownames(counts), perl = TRUE)
  counts  <- counts[!is_ctrl, , drop = FALSE]
  nCount  <- colSums(counts)

  assigned_tx <- sum(nCount)
  total_cells <- ncol(so)

  cell_area <- if (area_col %in% colnames(so@meta.data)) so@meta.data[[area_col]] else rep(NA_real_, ncol(so))

  data.frame(
    sample = sample, platform = platform, segmentation = segmentation,
    total_transcripts           = total_transcripts,
    assigned_transcripts        = assigned_tx,
    total_cells                 = total_cells,
    median_cell_area_um2        = round(median(cell_area, na.rm = TRUE), 2),
    median_transcripts_per_cell = round(median(nCount),   1),
    stringsAsFactors = FALSE
  )
}

# ---- Step 1: total (pre-segmentation) transcripts per sample/platform -----
# shared by Default and Proseg (both count from the same raw file); Baysor
# reuses the same totals too (its own segmentation.csv row count is not used)

xenium_total <- sapply(samples_xenium, function(s)
  fread(xenium_tx_paths[s], select = c("qv", "feature_name"))[
    qv >= 20 & !grepl(control_pattern, feature_name, perl = TRUE), .N])

merscope_total <- sapply(samples_merscope, function(s)
  fread(merscope_tx_paths[s], select = "gene")[
    !grepl(control_pattern, gene, perl = TRUE), .N])

# ---- Step 2: per-dataset QC metrics ----------------------------------------

res_xenium_default <- lapply(samples_xenium, function(s)
  compute_qc_row(readRDS(xenium_default_paths[s]), s, "Xenium", "Default", xenium_total[s], area_col = "cell_area"))

res_merscope_default <- lapply(samples_merscope, function(s)
  compute_qc_row(load_merscope_seu_2d_area(s), s, "MERSCOPE", "Default", merscope_total[s], area_col = "cell_area"))

res_baysor_xenium <- lapply(samples_xenium, function(s)
  compute_qc_row(readRDS(file.path(baysor_xenium_dir, paste0(s, "_so_raw.rds"))), s, "Xenium", "Baysor",
                  xenium_total[s], area_col = "area", counts_layer = "counts.Gene Expression"))

res_baysor_merscope <- lapply(samples_merscope, function(s)
  compute_qc_row(readRDS(file.path(baysor_merscope_dir, paste0(s, "_so_raw.rds"))), s, "MERSCOPE", "Baysor",
                  merscope_total[s], area_col = "area", counts_layer = "counts.Gene Expression"))

res_proseg_xenium <- lapply(samples_xenium, function(s)
  compute_qc_row(load_proseg_seu_2d_area(s, "xenium"), s, "Xenium", "Proseg", xenium_total[s], area_col = "cell_area"))

res_proseg_merscope <- lapply(samples_merscope, function(s)
  compute_qc_row(load_proseg_seu_2d_area(s, "merscope"), s, "MERSCOPE", "Proseg", merscope_total[s], area_col = "cell_area"))

# ---- Step 3: assemble and save ---------------------------------------------

qc_df <- do.call(rbind, c(
  res_xenium_default, res_merscope_default,
  res_baysor_xenium, res_baysor_merscope,
  res_proseg_xenium, res_proseg_merscope
))
rownames(qc_df) <- NULL
qc_df <- qc_df %>% arrange(platform, segmentation, sample)

saveRDS(qc_df, "qc_df.rds")
