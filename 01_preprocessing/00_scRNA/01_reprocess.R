# Purpose:  Build a Seurat object from one raw sc/snRNA-seq sample: load the Cell Ranger
#           matrix, apply cell-level QC, then normalise, embed and cluster, and print a
#           top-10-marker heatmap for manual cluster inspection. Stage 1 of the reference
#           arm — its clusters are annotated in 03_annotate_cell_types.R.
#           Run once per sample; `sample_id` and the three QC thresholds are set by hand
#           at the top — see the table in README.md.
# Inputs:   raw_scrna_dir[<sample_id>]   Cell Ranger matrix directory, set in config.R
#                                              (barcodes/features/matrix)
# Outputs:  seu.rds   QC-filtered, normalised, embedded and clustered object,
#                     written to the working directory
#           A marker heatmap is drawn to the active graphics device, not saved.

library(Seurat)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

sample_id <- "TNBC_02"   # any ID in `sample_order`; see config.R

# QC thresholds are sample-specific; see the table in README.md.
min_nFeature <- NA
max_nCount <- NA
max_percent_mt <- NA

## load data
data = Read10X(raw_scrna_dir[[sample_id]])
seu = CreateSeuratObject(data)

## Quality Control
seu$percent.mt = PercentageFeatureSet(seu, pattern = "^MT-")

# We remove cells with low detected genes, high total counts and high mitochondrial gene expression.
seu = seu[, seu$nFeature_RNA > min_nFeature & seu$nCount_RNA < max_nCount & seu$percent.mt < max_percent_mt]

## Data exploration
seu = NormalizeData(seu)

seu = FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

seu = ScaleData(seu)

seu = RunPCA(seu, verbose = F, features = VariableFeatures(seu))
seu = RunUMAP(seu, dims = 1:25)

seu = FindNeighbors(seu, dims = 1:25)

seu = FindClusters(seu, resolution = 0.3)

# The heatmap shows the top 10 markers of each cluster.
Idents(seu) <- seu$seurat_clusters
markers <- FindAllMarkers(seu, only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)

markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(seu, features = top10$gene) + NoLegend()

saveRDS(seu, "seu.rds")
