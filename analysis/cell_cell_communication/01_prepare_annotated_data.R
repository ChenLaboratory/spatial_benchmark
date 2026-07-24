# Cross-platform, pathology-annotated Seurat objects for one sample (MH0026 /
# "TNBC-02"), covering the full Visium tissue region (not restricted to the
# Xenium-MERSCOPE overlap used elsewhere in this repo). Feeds
# 02_run_blisa.R, 03_run_cci.R, 04_region_lr_proportion.R and figures/Fig5.R
# (panels a, d, e, f). See README.md.

library(Seurat)
library(dplyr)
library(sf)

sample <- "MH0026"  # "TNBC-02" in the manuscript

raw_data_file <- "/path/to/analysis/process_aligned_data/raw_data_list_AllSample.rds"
visium_qc_file  <- "/path/to/preprocessing/Visium/MH0026/seu_ST_Filtered.rds"
xenium_qc_file   <- "/path/to/preprocessing/Xenium/MH0026/so.rds"
merscope_qc_file <- "/path/to/preprocessing/MERSCOPE/MH0026/so.rds"
# QuPath cell-level exports, one CSV per pathologist-annotated region;
# each file starts with a "#Selection name: <region>" header line
annotation_dir <- "/path/to/raw/pathology_annotation/MH0026"

## ---- Step 1: raw data with alignment results (no QC) ----------------------
## (Visium_spot_id, aligned_x, aligned_y_reversed, aligned_spatial reduction
## already added by process_aligned_data)

raw_data_list_AllSample <- readRDS(raw_data_file)
raw_x <- raw_data_list_AllSample[[sample]]$Xenium
raw_m <- raw_data_list_AllSample[[sample]]$MERSCOPE
raw_v <- raw_data_list_AllSample[[sample]]$Visium

## ---- Step 2: cells passing QC ----------------------------------------------

visium_proc <- readRDS(visium_qc_file)
xen_proc    <- readRDS(xenium_qc_file)
mer_proc    <- readRDS(merscope_qc_file)

## ---- Step 3: subset Xenium/MERSCOPE to cells covered by Visium spots ------
## (all iST cells in the Visium tissue region, not just the cross-platform
## overlapped hexbins used for benchmarking)

valid_spots <- colnames(visium_proc)

seu_x <- raw_x[, intersect(colnames(raw_x), colnames(xen_proc))]
seu_x$Visium_spot_id <- ifelse(seu_x$Visium_spot_id %in% valid_spots, seu_x$Visium_spot_id, "")
seu_x <- seu_x[, seu_x$Visium_spot_id != ""]

seu_m <- raw_m[, intersect(colnames(raw_m), colnames(mer_proc))]
seu_m$Visium_spot_id <- ifelse(seu_m$Visium_spot_id %in% valid_spots, seu_m$Visium_spot_id, "")
seu_m <- seu_m[, seu_m$Visium_spot_id != ""]

seu_v <- visium_proc

seu_x$SingleR_labels202605 <- xen_proc$SingleR_labels202605[colnames(seu_x)]
seu_m$SingleR_labels202605 <- mer_proc$SingleR_labels202605[colnames(seu_m)]

## ---- Pathology annotation ---------------------------------------------------
## region = specific annotated region (e.g. "Tumour_more_cohesive_R1");
## type   = region with the trailing "_R<n>" replicate suffix stripped

annotation_files <- list.files(annotation_dir, pattern = "_cells_stats\\.csv$", full.names = TRUE)

read_annotation_file <- function(file) {
  lines <- readLines(file)
  region_name <- trimws(sub("^#Selection name:\\s*", "", lines[1]))
  region_type <- sub("_[Rr]\\d+$", "", region_name)
  region_type <- gsub("\\s+", "_", region_type)
  region_type <- paste0(toupper(substr(region_type, 1, 1)), substr(region_type, 2, nchar(region_type)))
  region_name <- gsub("\\s+", "_", region_name)
  region_name <- paste0(toupper(substr(region_name, 1, 1)), substr(region_name, 2, nchar(region_name)))
  df <- read.csv(textConnection(paste(lines[3:length(lines)], collapse = "\n")), check.names = FALSE)
  df$pathology_region <- region_name
  df$pathology_type   <- region_type
  df
}

annotation_df <- do.call(rbind, lapply(annotation_files, read_annotation_file))
colnames(annotation_df)[colnames(annotation_df) == "Cell ID"] <- "cell_id"

# 1. Xenium: matched directly by cell ID; unmatched cells are "Unannotated"
xen_anno_df <- data.frame(
  pathology_region = annotation_df$pathology_region,
  pathology_type   = annotation_df$pathology_type,
  row.names        = annotation_df$cell_id
)
seu_x <- AddMetaData(seu_x, metadata = xen_anno_df)
seu_x$pathology_region[is.na(seu_x$pathology_region)] <- "Unannotated"
seu_x$pathology_type[is.na(seu_x$pathology_type)]     <- "Unannotated"

# 2. Visium: majority vote of annotated Xenium cells per spot
xen_annotated <- seu_x@meta.data %>%
  filter(pathology_region != "Unannotated") %>%
  select(Visium_spot_id, pathology_region, pathology_type)

spot_to_region <- xen_annotated %>%
  group_by(Visium_spot_id) %>%
  summarise(
    pathology_region = names(which.max(table(pathology_region))),
    pathology_type   = names(which.max(table(pathology_type))),
    n_xenium_cells   = n(),
    .groups = "drop"
  )

vis_anno_df <- data.frame(
  pathology_region = spot_to_region$pathology_region,
  pathology_type   = spot_to_region$pathology_type,
  row.names        = spot_to_region$Visium_spot_id
)
seu_v <- AddMetaData(seu_v, metadata = vis_anno_df)
seu_v$pathology_region[is.na(seu_v$pathology_region)] <- "Unannotated"
seu_v$pathology_type[is.na(seu_v$pathology_type)]     <- "Unannotated"

# 3. MERSCOPE: inherits its assigned Visium spot's majority-vote annotation
spot_region_map <- setNames(spot_to_region$pathology_region, spot_to_region$Visium_spot_id)
spot_type_map   <- setNames(spot_to_region$pathology_type,   spot_to_region$Visium_spot_id)

mer_regions <- spot_region_map[seu_m$Visium_spot_id]
mer_types   <- spot_type_map[seu_m$Visium_spot_id]
names(mer_regions) <- names(mer_types) <- colnames(seu_m)

seu_m <- AddMetaData(seu_m, metadata = mer_regions, col.name = "pathology_region")
seu_m <- AddMetaData(seu_m, metadata = mer_types,   col.name = "pathology_type")
seu_m$pathology_region[is.na(seu_m$pathology_region)] <- "Unannotated"
seu_m$pathology_type[is.na(seu_m$pathology_type)]     <- "Unannotated"

## ---- save ------------------------------------------------------------------

saveRDS(annotation_df, "annotation_df.rds")  # cell-level QuPath table, used by Fig5 for region/type ordering
saveRDS(seu_x, "seu_x_annotated.rds")
saveRDS(seu_m, "seu_m_annotated.rds")
saveRDS(seu_v, "seu_v_annotated.rds")
