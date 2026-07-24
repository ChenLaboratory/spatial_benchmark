# Shared QC + clustering + SingleR annotation step for Baysor and Proseg v3
# Seurat objects (both platforms) -- the step referenced as step 3 in
# ../README.md. Produces the "<sample>_so.rds" objects that
# analysis/cell_segmentation_metrics and figures/Fig4.R (panels e/f/g) read.
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

nCount_RNA_cutoff   <- 10
nFeature_RNA_cutoff <- 5
control_pattern     <- "(?i)^(Blank|FC\\d+$|Sabrina|NegControl)"

baysor_xenium_dir    <- "/path/to/run_cell_segmentation/baysor/results/xenium"      # <sample>_so_raw.rds -> <sample>_so.rds
baysor_merscope_dir  <- "/path/to/run_cell_segmentation/baysor/results/merscope"    # <sample>_so_raw.rds -> <sample>_so.rds
proseg_xenium_dir     <- "/path/to/run_cell_segmentation/proseg_v3/results/xenium"    # <sample>_so_raw.rds -> <sample>_so.rds
proseg_merscope_dir   <- "/path/to/run_cell_segmentation/proseg_v3/results/merscope"  # <sample>_so_raw.rds -> <sample>_so.rds
sc_dir <- "/path/to/preprocessing/scRNA"  # <ref_id>/seu_<ref_id>.rds (snRNA for MH0007/MH0026, scRNA otherwise); cell_type202605 labels

# method x platform x sample manifest; ref_id is the reference sample ID as it
# appears in the reference file path (underscore style, e.g. ER_0114)
samples_meta <- tibble::tribble(
  ~method,      ~platform,   ~sample,   ~ref_id,
  "baysor",     "merscope",  "ER0114",  "ER_0114",
  "baysor",     "merscope",  "ER0360",  "ER_0360",
  "baysor",     "merscope",  "MH0007",  "MH0007",
  "baysor",     "merscope",  "MH0026",  "MH0026",
  "baysor",     "merscope",  "TN0177",  "TN_0177",
  "baysor",     "xenium",    "ER0114",  "ER_0114",
  "baysor",     "xenium",    "ER0360",  "ER_0360",
  "baysor",     "xenium",    "MH0007",  "MH0007",
  "baysor",     "xenium",    "MH0026",  "MH0026",
  "baysor",     "xenium",    "TN0177",  "TN_0177",
  "baysor",     "xenium",    "TN0554",  "TN_0554",
  "proseg_v3",  "merscope",  "ER0114",  "ER_0114",
  "proseg_v3",  "merscope",  "ER0360",  "ER_0360",
  "proseg_v3",  "merscope",  "MH0007",  "MH0007",
  "proseg_v3",  "merscope",  "MH0026",  "MH0026",
  "proseg_v3",  "merscope",  "TN0177",  "TN_0177",
  "proseg_v3",  "xenium",    "ER0114",  "ER_0114",
  "proseg_v3",  "xenium",    "ER0360",  "ER_0360",
  "proseg_v3",  "xenium",    "MH0007",  "MH0007",
  "proseg_v3",  "xenium",    "MH0026",  "MH0026",
  "proseg_v3",  "xenium",    "TN0177",  "TN_0177",
  "proseg_v3",  "xenium",    "TN0554",  "TN_0554"
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
