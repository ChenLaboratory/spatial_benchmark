# Human breast cancer Xenium / MERSCOPE preprocessing

library(Seurat)

platform <- "Xenium" # or "MERSCOPE"
sample_id <- "sample_id"
data_path <- "/path/to/platform/output"

# dimensionality (RunPCA/RunUMAP/FindNeighbors) and clustering resolution are
# platform/sample-specific; see the table in README.md
dimUsed <- NA
resolution <- NA

## Read in the data

if (platform == "Xenium") {
  raw_data <- Read10X_h5(paste0(data_path, "/cell_feature_matrix.h5"))
  counts <- raw_data[[1]]

  cell_info <- read.csv(paste0(data_path, "/cells.csv.gz"))
  rownames(cell_info) <- cell_info$cell_id
  cell_info <- cell_info[, -1]

  x_col <- "x_centroid"
  y_col <- "y_centroid"
  size_col <- "cell_area"
} else if (platform == "MERSCOPE") {
  raw_data <- read.csv(paste0(data_path, "/cell_by_gene.csv"))
  raw_data_chr <- read.csv(paste0(data_path, "/cell_by_gene.csv"), colClasses = "character")
  rownames(raw_data) <- raw_data_chr$cell
  counts <- t(raw_data[, -1])

  cell_info <- read.csv(paste0(data_path, "/cell_metadata.csv"))
  cell_info_chr <- read.csv(paste0(data_path, "/cell_metadata.csv"), colClasses = "character")
  rownames(cell_info) <- cell_info_chr$EntityID
  cell_info <- cell_info[, -1]

  x_col <- "center_x"
  y_col <- "center_y"
  size_col <- "volume"
}

so <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = cell_info,
                         min.cells = 0, min.features = 0)

# add spatial map
sp_embeddings <- so@meta.data[, c(x_col, y_col)]
colnames(sp_embeddings) <- c("Sp_1", "Sp_2")
spatial <- CreateDimReducObject(embeddings = as.matrix(sp_embeddings), assay = "RNA", key = "Sp_")
so@reductions[["spatial"]] <- spatial

# add feature
so$density <- so$nCount_RNA / so@meta.data[[size_col]]

saveRDS(so, "so_raw.rds")

## Quality control
# Cell-level QC thresholds are fixed across all samples and both platforms (see README.md).

nCount_RNA_cutoff <- 10
nFeature_RNA_cutoff <- 5

keep.cell <- so@meta.data$nCount_RNA >= nCount_RNA_cutoff & so@meta.data$nFeature_RNA >= nFeature_RNA_cutoff
so <- so[, keep.cell]

## Data exploration

so <- NormalizeData(so, scale.factor = 100)
so <- ScaleData(so, features = rownames(so))
so <- RunPCA(so, features = rownames(so), dims = 1:dimUsed)
so <- RunUMAP(so, dims = 1:dimUsed, seed.use = 2021, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:dimUsed)
so <- FindClusters(so, resolution = resolution)

saveRDS(so, "so.rds")
