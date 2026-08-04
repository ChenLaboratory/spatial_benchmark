# Purpose:  Shared QC, clustering and SingleR annotation step for the Baysor and Proseg v3
#           Seurat objects, both platforms — step 3 in ../README.md. Removes control
#           features, applies the cell-level QC filter, runs normalise/scale/PCA/UMAP/
#           cluster, transfers cell type labels from the matched sc/snRNA-seq reference,
#           then saves. Loops over the method x platform x sample manifest defined below;
#           samples whose input is missing are skipped with a message.
# Inputs:   <processed_data_dir>/segmentation/<method>/<platform>/<sample>_so_raw.rds
#                     raw segmented objects from baysor/ and proseg_v3/
#           <processed_data_dir>/scRNA/<ref_id>/seu_<ref_id>.rds
#                     annotated reference (snRNA for TNBC_01/TNBC_02, scRNA otherwise),
#                     supplying the cell_type202605 labels
# Outputs:  <processed_data_dir>/segmentation/<method>/<platform>/<sample>_so.rds
#                     read by analysis/segmentation_metrics/ and
#                     figures/Fig4_segmentation.R (panels e/f/g)
#
# Trimmed to the steps that feed those figures: control-feature removal, QC
# filter, normalise/scale/PCA/UMAP/cluster, SingleR annotation, save. The
# source analysis also produced a battery of diagnostic QC/cluster/marker
# plots at this step (histograms, spatial QC maps, per-cluster count
# distributions, a marker dot plot, cell-morphology-by-type boxplots) -- none
# of those feed a manuscript figure, so they're not reproduced here.

library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(dplyr)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

nCount_RNA_cutoff   <- 10
nFeature_RNA_cutoff <- 5
control_pattern     <- "(?i)^(Blank|FC\\d+$|Sabrina|NegControl)"

baysor_xenium_dir    <- file.path(processed_data_dir, "segmentation", "baysor", "Xenium")      # <sample>_so_raw.rds -> <sample>_so.rds
baysor_merscope_dir  <- file.path(processed_data_dir, "segmentation", "baysor", "MERSCOPE")    # <sample>_so_raw.rds -> <sample>_so.rds
proseg_xenium_dir     <- file.path(processed_data_dir, "segmentation", "proseg_v3", "Xenium")    # <sample>_so_raw.rds -> <sample>_so.rds
proseg_merscope_dir   <- file.path(processed_data_dir, "segmentation", "proseg_v3", "MERSCOPE")  # <sample>_so_raw.rds -> <sample>_so.rds
sc_dir <- file.path(processed_data_dir, "scRNA")  # <ref_id>/seu_<ref_id>.rds (snRNA for TNBC_01/TNBC_02, scRNA otherwise); cell_type202605 labels

# method x platform x sample manifest; ref_id is the reference sample ID as it
# appears in the reference file path (e.g. ER_02)
samples_meta <- tibble::tribble(
  ~method,      ~platform,   ~sample,   ~ref_id,
  "baysor",     "merscope",  "ER_02",  "ER_02",
  "baysor",     "merscope",  "ER_01",  "ER_01",
  "baysor",     "merscope",  "TNBC_01",  "TNBC_01",
  "baysor",     "merscope",  "TNBC_02",  "TNBC_02",
  "baysor",     "merscope",  "TNBC_03",  "TNBC_03",
  "baysor",     "xenium",    "ER_02",  "ER_02",
  "baysor",     "xenium",    "ER_01",  "ER_01",
  "baysor",     "xenium",    "TNBC_01",  "TNBC_01",
  "baysor",     "xenium",    "TNBC_02",  "TNBC_02",
  "baysor",     "xenium",    "TNBC_03",  "TNBC_03",
  "baysor",     "xenium",    "TNBC_04",  "TNBC_04",
  "proseg_v3",  "merscope",  "ER_02",  "ER_02",
  "proseg_v3",  "merscope",  "ER_01",  "ER_01",
  "proseg_v3",  "merscope",  "TNBC_01",  "TNBC_01",
  "proseg_v3",  "merscope",  "TNBC_02",  "TNBC_02",
  "proseg_v3",  "merscope",  "TNBC_03",  "TNBC_03",
  "proseg_v3",  "xenium",    "ER_02",  "ER_02",
  "proseg_v3",  "xenium",    "ER_01",  "ER_01",
  "proseg_v3",  "xenium",    "TNBC_01",  "TNBC_01",
  "proseg_v3",  "xenium",    "TNBC_02",  "TNBC_02",
  "proseg_v3",  "xenium",    "TNBC_03",  "TNBC_03",
  "proseg_v3",  "xenium",    "TNBC_04",  "TNBC_04"
)

so_path <- function(method, platform, sample) {
  dir <- switch(paste(method, platform),
    "baysor merscope"    = baysor_merscope_dir,
    "baysor xenium"      = baysor_xenium_dir,
    "proseg_v3 merscope" = proseg_merscope_dir,
    "proseg_v3 xenium"   = proseg_xenium_dir
  )
  file.path(dir, paste0(sample, "_so_raw.rds"))
}
out_path <- function(method, platform, sample) sub("_so_raw\\.rds$", "_so.rds", so_path(method, platform, sample))

for (i in seq_len(nrow(samples_meta))) {
  method   <- samples_meta$method[i]
  platform <- samples_meta$platform[i]
  sample   <- samples_meta$sample[i]
  ref_id   <- samples_meta$ref_id[i]

  path <- so_path(method, platform, sample)
  if (!file.exists(path)) { message("so_raw not found, skipping: ", path); next }
  message(method, " | ", platform, " | ", sample)

  # ---- control-feature removal ---------------------------------------------
  so <- readRDS(path)
  control_genes <- grep(control_pattern, rownames(so), value = TRUE, perl = TRUE)
  so <- so[!rownames(so) %in% control_genes, ]

  # ---- QC filter -------------------------------------------------------------
  keep_cell <- so@meta.data$nCount_RNA >= nCount_RNA_cutoff & so@meta.data$nFeature_RNA >= nFeature_RNA_cutoff
  so <- so[, keep_cell]

  # ---- Seurat pipeline: normalise -> scale -> PCA -> UMAP -> cluster --------
  # scale.factor = 100 (transcripts per 100), following the Xenium pipeline convention;
  # all genes used for PCA (targeted panel, no HVG selection needed)
  options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB -- large objects passed to FindNeighbors workers
  so <- NormalizeData(so, scale.factor = 100)
  so <- ScaleData(so, features = rownames(so))
  dimUsed <- 30
  so <- RunPCA(so, features = rownames(so), dims = 1:dimUsed)
  so <- RunUMAP(so, dims = 1:dimUsed, seed.use = 2021, verbose = FALSE)
  so <- FindNeighbors(so, dims = 1:dimUsed)
  so <- FindClusters(so, resolution = 0.3)

  # ---- SingleR annotation ----------------------------------------------------
  # matched sc/snRNA-seq reference, cell_type202605 labels; "mixture" (ambiguous
  # doublet) cells excluded from the reference
  ref_seu <- readRDS(file.path(sc_dir, ref_id, paste0("seu_", ref_id, ".rds")))
  ref_sce <- as.SingleCellExperiment(ref_seu)
  ref_sce <- ref_sce[, ref_sce$cell_type202605 != "mixture"]

  sce  <- as.SingleCellExperiment(so)
  pred <- SingleR::SingleR(test = sce, ref = ref_sce, labels = ref_sce$cell_type202605, de.method = "wilcox")
  so$SingleR_labels202605 <- pred$labels

  # ---- save -------------------------------------------------------------------
  saveRDS(so, out_path(method, platform, sample))
}
