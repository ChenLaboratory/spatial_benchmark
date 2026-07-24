# Human breast cancer Visium preprocessing

sample_id <- "sample_id"
visium_data_path <- "/path/to/spaceranger/outs"

library(Seurat)

# QC thresholds and clustering resolution are sample-specific; see the table in README.md.
min_nCount_Spatial <- NA
min_nFeature_Spatial <- NA
max_percent_mt <- NA
resolution <- NA

## Preliminary analysis

### Read in the data

seu <- Load10X_Spatial(visium_data_path, filter.matrix = FALSE)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

saveRDS(seu, "so_raw.rds")

### Quality control

#### At spot level
# We filter spots by the number of reads, the number of expressed genes and the percentage of
# mitochondrial gene expression.

ST <- subset(seu, subset = nCount_Spatial > min_nCount_Spatial & nFeature_Spatial > min_nFeature_Spatial & percent.mt < max_percent_mt)

#### At gene level
# We also removed lowly expressed genes as follows.

ST_Filtered <- ST

counts <- ST@assays$Spatial@layers$counts
keep1 <- rowSums(counts > 0) > 10 # expressed in at least 10 spots
keep2 <- rowSums(counts) > 30 # total detected expression > 30
keep <- keep1 & keep2
genes.use <- rownames(ST)[keep]
ST_Filtered <- ST[genes.use, ]

### Seurat pre-processing

ST_Filtered <- NormalizeData(ST_Filtered)
ST_Filtered <- FindVariableFeatures(ST_Filtered, selection.method = "vst", nfeatures = 2000)
ST_Filtered <- ScaleData(ST_Filtered, features = rownames(ST_Filtered))

### Spot clustering

ST_Filtered <- RunPCA(ST_Filtered)
ST_Filtered <- RunUMAP(ST_Filtered, reduction = "pca", dims = 1:30)
ST_Filtered <- FindNeighbors(ST_Filtered, reduction = "pca", dims = 1:30)

ST_Filtered <- FindClusters(ST_Filtered, verbose = FALSE, resolution = resolution)

saveRDS(ST_Filtered, "seu_ST_Filtered.rds")
