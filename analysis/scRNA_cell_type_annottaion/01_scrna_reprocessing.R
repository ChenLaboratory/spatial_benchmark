# Human breast cancer sc/snRNA-seq preprocessing

library(Seurat)

# QC thresholds are sample-specific; see the table in README.md.
min_nFeature <- NA
max_nCount <- NA
max_percent_mt <- NA

## load data
data = Read10X("/path/to/sample")
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
