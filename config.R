## 1. Directories -------------------------------------------------------------

processed_data_dir <- "processed_data_dir"  # Seurat objects written by preprocessing/
results_table_dir  <- "results_table_dir"   # result tables written by analysis/
figure_output_dir  <- "figure_output_dir"   # figure files written by figures/


## 2. Samples -----------------------------------------------------------------

# Canonical display order. Every factor(levels = ) in the repo derives from this,
# so changing the order here changes panel and legend order in every figure.
sample_order <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "TNBC_05",
                  "ER_01", "ER_02")

bench_samples <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "ER_01", "ER_02")

xenium_samples   <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "ER_01", "ER_02")
visium_samples   <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "ER_01", "ER_02")
merscope_samples <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "TNBC_05",
                      "ER_01", "ER_02")

merscope_v1_samples <- c("TNBC_01", "TNBC_02", "TNBC_03", "TNBC_04", "TNBC_05",
                         "ER_01", "ER_02")
merscope_v2_samples <- c("TNBC_03", "TNBC_05", "ER_01", "ER_02")

bench_samples_xenium   <- intersect(bench_samples, xenium_samples)
bench_samples_merscope <- c("TNBC_01", "TNBC_02", "TNBC_03", "ER_01", "ER_02")

bench_samples_merscope_v1 <- c("TNBC_01", "TNBC_02")
bench_samples_merscope_v2 <- c("TNBC_03", "ER_01", "ER_02")


## 3. Colours and factor levels -----------------------------------------------

sample_colors <- c(TNBC_01 = "#1565C0", TNBC_02 = "#C0392B", TNBC_03 = "#2E7D32",
                   TNBC_04 = "#4E342E", TNBC_05 = "#00838F",
                   ER_01   = "#AD1457", ER_02   = "#F9A825")

platform_colors <- c(
  MERSCOPE_V1 = "#80B1D3",
  MERSCOPE_V2 = "#0b6170",
  Xenium      = "#FDAE61",
  Visium      = "#BC80BD"
)
platform_order <- c("MERSCOPE_V1", "MERSCOPE_V2", "Xenium", "Visium")

platform_shapes <- c(Xenium = 16, MERSCOPE = 17, scRNA = 15)

cell_type_colors <- c(
  Tumor = "#1F77B4", Myeloid = "#AD8BC9", T = "#D62728", Endothelial = "#8C564B",
  Pericyte = "#C49C94", Fibroblast = "#FF7F0E", Adipocyte = "#E377C2", B = "#BCBD22",
  Plasma = "#17BECF", Epithelial = "#67BF5C", Lymph.vessel = "#A2A2A2",
  unsure = "black", mixture = "grey50"
)
cell_type_order <- names(cell_type_colors)


## 4. Raw data directories ----------------------------------------------------
#
# One directory per (sample, dataset). Point each at wherever that dataset lives
# -- nothing above the sample directory is assumed, so they need not share a
# parent. Inside each directory the vendor's standard output layout is assumed,
# and the scripts read the usual filenames from it:
#
#   Xenium     cell_feature_matrix.h5, cells.csv.gz,
#              transcripts.parquet, transcripts.csv.gz
#   MERSCOPE   cell_by_gene.csv, cell_metadata.csv,
#              detected_transcripts.csv, cell_boundaries.parquet
#   Visium     the Space Ranger sample directory, containing outs/
#   sc/snRNA   the Cell Ranger matrix directory (barcodes/features/matrix)

raw_xenium_dir <- c(
  TNBC_01 = "TNBC_01_xenium_raw_data_dir",
  TNBC_02 = "TNBC_02_xenium_raw_data_dir",
  TNBC_03 = "TNBC_03_xenium_raw_data_dir",
  TNBC_04 = "TNBC_04_xenium_raw_data_dir",
  ER_01   = "ER_01_xenium_raw_data_dir",
  ER_02   = "ER_02_xenium_raw_data_dir"
)

# MERSCOPE V1 acquisitions. For TNBC_01/TNBC_02/TNBC_04 this is the sample's only
# MERSCOPE run; for the rest it is the secondary run shown only in Fig2 a-d.
raw_merscope_v1_dir <- c(
  TNBC_01 = "TNBC_01_merscope_v1_raw_data_dir",
  TNBC_02 = "TNBC_02_merscope_v1_raw_data_dir",
  TNBC_03 = "TNBC_03_merscope_v1_raw_data_dir",
  TNBC_04 = "TNBC_04_merscope_v1_raw_data_dir",
  TNBC_05 = "TNBC_05_merscope_v1_raw_data_dir",
  ER_01   = "ER_01_merscope_v1_raw_data_dir",
  ER_02   = "ER_02_merscope_v1_raw_data_dir"
)

# MERSCOPE V2 acquisitions.
raw_merscope_v2_dir <- c(
  TNBC_03 = "TNBC_03_merscope_v2_raw_data_dir",
  TNBC_05 = "TNBC_05_merscope_v2_raw_data_dir",
  ER_01   = "ER_01_merscope_v2_raw_data_dir",
  ER_02   = "ER_02_merscope_v2_raw_data_dir"
)

raw_visium_dir <- c(
  TNBC_01 = "TNBC_01_visium_raw_data_dir",
  TNBC_02 = "TNBC_02_visium_raw_data_dir",
  TNBC_03 = "TNBC_03_visium_raw_data_dir",
  TNBC_04 = "TNBC_04_visium_raw_data_dir",
  ER_01   = "ER_01_visium_raw_data_dir",
  ER_02   = "ER_02_visium_raw_data_dir"
)

raw_scrna_dir <- c(
  TNBC_01 = "TNBC_01_snrna_raw_data_dir",
  TNBC_02 = "TNBC_02_snrna_raw_data_dir",
  TNBC_03 = "TNBC_03_scrna_raw_data_dir",
  TNBC_04 = "TNBC_04_scrna_raw_data_dir",
  ER_01   = "ER_01_scrna_raw_data_dir",
  ER_02   = "ER_02_scrna_raw_data_dir"
)

# Pathology annotation polygons (Visium spot regions); TNBC_02 only.
raw_pathology_annotation_dir <- c(TNBC_02 = "TNBC_02_pathology_annotation_dir")

# Gene panel definition (single file), shared by both iST platforms.
raw_gene_panel_file <- "Xenium_380_gene_panel.xlsx"


## 5. Convenience -------------------------------------------------------------

# Each sample's MERSCOPE directory: its V2 run where one exists, otherwise V1.
raw_merscope_dir <- stats::setNames(
  ifelse(merscope_samples %in% merscope_v2_samples,
         raw_merscope_v2_dir[merscope_samples],
         raw_merscope_v1_dir[merscope_samples]),
  merscope_samples)

# Raw sample directory for a given iST platform, used by scripts that handle
# both (e.g. preprocessing/iST/01_preprocess_iST.R).
raw_ist_dir <- function(platform, sample_id) {
  stopifnot(platform %in% c("Xenium", "MERSCOPE"))
  lookup <- if (platform == "Xenium") raw_xenium_dir else raw_merscope_dir
  if (!sample_id %in% names(lookup))
    stop(sample_id, " has no ", platform, " data; see config.R section 2.")
  unname(lookup[[sample_id]])
}
