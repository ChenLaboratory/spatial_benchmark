# Purpose:  Assign manual cell type labels to the sc/snRNA-seq reference. One block per
#           sample (TNBC_01…04, ER_01, ER_02), each mapping that sample's Seurat clusters
#           to cell types from marker dot plots, sub-clustering ambiguous clusters and
#           using the 02_inferCNV.R instability score to call tumour clusters. The result
#           is the `cell_type202605` label that every other platform is annotated against.
#           Interactive/exploratory: run a sample's block at a time and inspect the plots,
#           rather than sourcing the file end to end.
# Inputs:   <sample>/seu_<sample>.rds   clustered object from 01_reprocess.R, per sample
#           ../../R/plot_dotplot.R      marker dot plot helper
#           inferCNV instability scores from 02_inferCNV.R (read by eye, not loaded here)
#           Sample paths are hardcoded in each block — edit them before running.
# Outputs:  <sample>/seu_<sample>.rds   same object, updated in place with the
#                     cell_type202605 label (plus intermediate cell_type<res> and
#                     subcluster columns)
#           DimPlots and dot plots are drawn to the active graphics device, not saved.

library(Seurat)
library(RColorBrewer)
library(magrittr)
library(dplyr)

cell_type_colors <- c(
  Tumor         = "#1F77B4",
  Myeloid       = "#AD8BC9",
  T             = "#D62728",
  Endothelial   = "#8C564B",
  Pericyte      = "#C49C94",
  Fibroblast    = "#FF7F0E",
  Adipocyte     = "#E377C2",
  B             = "#BCBD22",
  Plasma        = "#17BECF",
  Epithelial    = "#67BF5C",
  Lymph.vessel  = "#A2A2A2",
  unsure        = "black",
  B.Plasma      = "green"
)

ct_parent_groups <- list(
  Tumor        = c("Tumour"),
  Myeloid      = c("Monocyte", "Macrophage", "TAM", "cDC", "pDC", "Mast", "Neutrophil"),
  T            = c("CD4_T", "CD8_T", "Treg", "Ex_T", "NK"),
  Endothelial  = c("Endothelial"),
  Pericyte     = c("Pericyte"),
  Fibroblast   = c("Fibroblast", "iCAF", "myCAF", "apCAF"),
  Adipocyte    = c("Adipocyte"),
  B            = c("Naive_B", "Memory_B"),
  Plasma       = c("Plasma"),
  Epithelial   = c("Basal", "Luminal"),
  Lymph.vessel = c("Lymph.vessel")
)

source("../../R/plot_dotplot.R")

# ============================================================
# TNBC_01
# ============================================================

seu_TNBC_01 <- readRDS("/vast/projects/Spatial/lei/Benchmarking/snRNA/TNBC_01/seu_TNBC_01.rds")

## Initial annotation - res 0.6
# Uncertain clusters:
# - c2: mixture of Plasma, Basal, Fibroblast, Tumour -> assigned `unsure`
# - c12: Pericyte / Basal / Fibroblast markers -> assigned `Epithelial`

seu_TNBC_01@meta.data %<>% mutate(cell_type0.6 = case_match(RNA_snn_res.0.6,
  c("0", "1", "3", "5", "11") ~ "Tumor",   # supported by inferCNV
  "4"                          ~ "B.Plasma",
  "6"                          ~ "Myeloid",
  "8"                          ~ "T",
  "10"                         ~ "Endothelial",
  "12"                         ~ "Epithelial",
  "13"                         ~ "Adipocyte",
  c("7", "9")                  ~ "Fibroblast",
  "2"                          ~ "unsure"
))

DimPlot(seu_TNBC_01, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Higher resolution - res 1.2
# Reason: distinguish Pericyte from Endothelial.

seu_TNBC_01 <- FindClusters(seu_TNBC_01, resolution = 1.2, verbose = FALSE)
DimPlot(seu_TNBC_01, group.by = "RNA_snn_res.1.2",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu_TNBC_01, "TNBC_01", group_by = "RNA_snn_res.1.2")

seu_TNBC_01@meta.data %<>% mutate(cell_type1.2 = case_when(
  RNA_snn_res.1.2 == "15" ~ "Pericyte",
  TRUE                     ~ cell_type0.6
))

DimPlot(seu_TNBC_01, group.by = "cell_type1.2", cols = cell_type_colors, label = TRUE)

# Full re-annotation at res 1.2:

seu_TNBC_01@meta.data %<>% mutate(cell_type1.2 = case_match(RNA_snn_res.1.2,
  c("1", "2", "3", "5", "6", "7") ~ "Tumor",
  "10"                             ~ "T",
  "8"                              ~ "Myeloid",
  "4"                              ~ "Plasma",
  "12"                             ~ "Endothelial",
  "15"                             ~ "Pericyte",
  "14"                             ~ "Adipocyte",
  "13"                             ~ "Epithelial",
  c("11", "9")                     ~ "Fibroblast",
  "0"                              ~ "mixture"
))

DimPlot(seu_TNBC_01, group.by = "cell_type1.2", cols = cell_type_colors, label = TRUE)

## Subcluster B.Plasma (cluster 4 at res 1.2)
# Attempting to separate B cells from Plasma. Result: subclustering failed to resolve a clean B population.

Idents(seu_TNBC_01) <- seu_TNBC_01$RNA_snn_res.1.2

seu_TNBC_01 <- FindSubCluster(
  seu_TNBC_01,
  cluster        = 4,
  graph.name     = "RNA_snn",
  resolution     = 0.8,
  subcluster.name = "subcluster"
)

DimPlot(seu_TNBC_01, group.by = "subcluster",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu_TNBC_01, "TNBC_01", group_by = "subcluster")

## Save

seu_TNBC_01$cell_type202605 <- seu_TNBC_01$cell_type1.2
saveRDS(seu_TNBC_01, "/vast/projects/Spatial/lei/Benchmarking/snRNA/TNBC_01/seu_TNBC_01.rds")

# ============================================================
# TNBC_02
# ============================================================

seu <- readRDS("/vast/projects/Spatial/lei/Benchmarking/snRNA/TNBC_02/seu_TNBC_02.rds")

## Initial annotation - res 0.6

seu@meta.data %<>% mutate(cell_type0.6 = case_match(RNA_snn_res.0.6,
  c("0", "1", "2", "5", "10", "11") ~ "Tumor",   # supported by inferCNV
  "4"                                ~ "Fibroblast",
  "6"                                ~ "Myeloid",
  "7"                                ~ "T",
  "8"                                ~ "Endothelial",
  "9"                                ~ "Plasma",
  "13"                               ~ "Pericyte",
  "14"                               ~ "Adipocyte",
  c("3", "12")                       ~ "Epithelial" # used for Xenium validation
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Higher resolution - res 1.7
# Reason: attempt to distinguish Lymph.vessel from Endothelial, and Plasma from B cells.
# Result: Endothelial splits into 2 clusters but both co-express Endothelial and Lymph.vessel markers; B/Plasma not resolved.

seu <- FindClusters(seu, resolution = 1.7, verbose = FALSE)
DimPlot(seu, group.by = "RNA_snn_res.1.7",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

## Subcluster Myeloid (cluster 6 at res 0.6)

Idents(seu) <- seu$RNA_snn_res.0.6

seu <- FindSubCluster(
  seu,
  cluster         = 6,
  graph.name      = "RNA_snn",
  resolution      = 0.3,
  subcluster.name = "subcluster"
)

DimPlot(seu, group.by = "subcluster",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "TNBC_02", group_by = "subcluster")

# Subcluster 6_4 shows B cell markers -> reassigned.

seu@meta.data %<>% mutate(cell_type0.6 = case_when(
  subcluster == "6_4" ~ "B",
  TRUE                ~ cell_type0.6
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Save

seu$cell_type202605 <- seu$cell_type0.6
saveRDS(seu, "/vast/projects/Spatial/lei/Benchmarking/snRNA/TNBC_02/seu_TNBC_02.rds")

# ============================================================
# ER_02
# ============================================================

seu <- readRDS("/vast/projects/Spatial/lei/Benchmarking/scRNA/ER_02/seu_ER_02.rds")

## Initial annotation - res 0.6
# Cluster notes:
# - c3, c12: co-express Basal, Fibroblast, Pericyte markers -> assigned `Epithelial`
# - c0, c1, c2, c4, c6, c10: Tumour / Epithelial markers

seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
DimPlot(seu, group.by = "RNA_snn_res.0.6",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

seu@meta.data %<>% mutate(cell_type0.6 = case_match(RNA_snn_res.0.6,
  # confident
  "9"                              ~ "Fibroblast",
  c("8", "14")                     ~ "Myeloid",
  "13"                             ~ "Plasma",
  "11"                             ~ "Endothelial",
  "5"                              ~ "T",
  # less confident
  "7"                              ~ "Pericyte",
  c("0", "1", "2", "4", "6", "10") ~ "Tumor",
  c("3", "12")                     ~ "Epithelial"
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Subcluster Plasma (cluster 13 at res 0.6)

Idents(seu) <- seu$RNA_snn_res.0.6

seu <- FindSubCluster(
  seu,
  cluster         = 13,
  graph.name      = "RNA_snn",
  resolution      = 0.6,
  subcluster.name = "subcluster"
)

DimPlot(seu, group.by = "subcluster",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "ER_02", group_by = "subcluster")

# Subcluster 13_0 shows B cell markers -> reassigned.

seu@meta.data %<>% mutate(cell_type0.6 = case_when(
  subcluster == "13_0" ~ "B",
  TRUE                  ~ cell_type0.6
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Save

seu$cell_type202605 <- seu$cell_type0.6
saveRDS(seu, "/vast/projects/Spatial/lei/Benchmarking/scRNA/ER_02/seu_ER_02.rds")

# ============================================================
# ER_01
# ============================================================

seu <- readRDS("/vast/projects/Spatial/lei/Benchmarking/scRNA/ER_01/seu_ER_01.rds")

## Initial annotation - res 0.6
# Uncertain clusters:
# - c14: mixed Plasma / Myeloid markers -> assigned `Myeloid`

seu@meta.data %<>% mutate(cell_type0.6 = case_match(RNA_snn_res.0.6,
  # confident
  c("2", "3", "8")   ~ "Tumor",
  c("0", "1", "5")   ~ "T",
  "4"                ~ "B",
  "7"                ~ "Myeloid",
  "9"                ~ "Endothelial",
  "10"               ~ "Pericyte",
  "12"               ~ "Fibroblast",
  # less confident
  "13"               ~ "Plasma",
  # not confident
  c("11", "6")       ~ "Epithelial",
  "14"               ~ "Myeloid"
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

plot_marker_dotplot(seu, "ER_01", group_by = "RNA_snn_res.0.6")

## Save

seu$cell_type202605 <- seu$cell_type0.6
saveRDS(seu, "/vast/projects/Spatial/lei/Benchmarking/scRNA/ER_01/seu_ER_01.rds")

# ============================================================
# TNBC_03
# ============================================================

seu <- readRDS("/vast/projects/Spatial/lei/Benchmarking/scRNA/TNBC_03/seu_TNBC_03.rds")

## Initial annotation - res 0.3
# Uncertain clusters:
# - c9: co-express Pericyte / Fibroblast markers
# - c13: co-express Basal / Fibroblast / Pericyte markers
# - c14, c15: co-express Myeloid / Fibroblast / T markers

seu@meta.data %<>% mutate(cell_type0.3 = case_match(RNA_snn_res.0.3,
  # confident
  "6"              ~ "Tumor",
  c("1", "10")     ~ "Myeloid",
  c("0", "2", "3", "7") ~ "T",
  c("4", "8")      ~ "Plasma",
  "5"              ~ "B",
  "11"             ~ "Fibroblast",
  # less confident
  "12"             ~ "Endothelial",  # also shows Adipocyte markers
  "9"              ~ "Fibroblast",
  "13"             ~ "Epithelial",
  c("14", "15")    ~ "Myeloid"
))

DimPlot(seu, group.by = "cell_type0.3", cols = cell_type_colors, label = TRUE)

## Higher resolution - res 0.6
# Reason: distinguish Pericyte from Fibroblast.

seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)
DimPlot(seu, group.by = "RNA_snn_res.0.6",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "TNBC_03", group_by = "RNA_snn_res.0.6")

# Cluster 18 at res 0.6 shows Pericyte markers:

seu@meta.data %<>% mutate(cell_type0.6 = case_when(
  RNA_snn_res.0.6 == "18" ~ "Pericyte",
  TRUE                     ~ cell_type0.3
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

# Full re-annotation at res 0.6:

seu@meta.data %<>% mutate(cell_type0.6 = case_match(RNA_snn_res.0.6,
  # confident
  "7"                          ~ "Tumor",
  c("1", "10", "12", "19")     ~ "Myeloid",
  c("0", "2", "3", "4", "8")  ~ "T",
  c("5", "15", "9")            ~ "Plasma",
  "6"                          ~ "B",
  c("11", "14", "17")          ~ "Fibroblast",
  "13"                         ~ "Endothelial",  # also shows Adipocyte markers
  "16"                         ~ "Epithelial",
  "18"                         ~ "Pericyte"
))

DimPlot(seu, group.by = "cell_type0.6", cols = cell_type_colors, label = TRUE)

## Higher resolution - res 2.5
# Reason: look for Adipocytes. Result: cluster 12 (Endothelial/Adipocyte) did not split cleanly.

seu <- FindClusters(seu, resolution = 2.5, verbose = FALSE)
DimPlot(seu, group.by = "RNA_snn_res.2.5",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "TNBC_03", group_by = "RNA_snn_res.2.5")

## Subcluster Endothelial (cluster 13 at res 0.6)
# Reason: attempt to separate Adipocytes from Endothelial. Result: Adipocyte markers not resolved from Endothelial.

Idents(seu) <- seu$RNA_snn_res.0.6

seu <- FindSubCluster(
  seu,
  cluster         = 13,
  graph.name      = "RNA_snn",
  resolution      = 0.6,
  subcluster.name = "subcluster"
)

DimPlot(seu, group.by = "subcluster",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "TNBC_03", group_by = "subcluster")

## Save

seu$cell_type202605 <- seu$cell_type0.6
saveRDS(seu, "/vast/projects/Spatial/lei/Benchmarking/scRNA/TNBC_03/seu_TNBC_03.rds")

# ============================================================
# TNBC_04
# ============================================================

seu <- readRDS("/vast/projects/Spatial/lei/Benchmarking/scRNA/TNBC_04/seu_TNBC_04.rds")

## Initial annotation - res 0.3
# Uncertain clusters:
# - c13: Endothelial with Adipocyte co-expression

seu@meta.data %<>% mutate(cell_type0.3 = case_match(RNA_snn_res.0.3,
  # confident
  c("1", "4")        ~ "Tumor",
  c("0", "10")       ~ "T",
  c("2", "5", "11")  ~ "Myeloid",
  "3"                ~ "B",
  c("6", "7")        ~ "Fibroblast",
  "8"                ~ "Plasma",
  "13"               ~ "Endothelial",  # also shows Adipocyte markers
  # less confident
  c("9", "12", "14") ~ "Epithelial"
))

DimPlot(seu, group.by = "cell_type0.3", cols = cell_type_colors, label = TRUE)

## Higher resolution - res 0.7
# Reason: distinguish Pericyte from other stromal populations.

DimPlot(seu, group.by = "RNA_snn_res.0.7",
        reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5)

plot_marker_dotplot(seu, "TNBC_04", group_by = "RNA_snn_res.0.7")

# Cluster 17 at res 0.7 shows Pericyte markers:

seu@meta.data %<>% mutate(cell_type0.7 = case_when(
  RNA_snn_res.0.7 == "17" ~ "Pericyte",
  TRUE                     ~ cell_type0.3
))

DimPlot(seu, group.by = "cell_type0.7", cols = cell_type_colors, label = TRUE)

# Full re-annotation at res 0.7:

seu@meta.data %<>% mutate(cell_type0.7 = case_match(RNA_snn_res.0.7,
  # confident
  c("1", "3", "12")          ~ "Tumor",
  c("2", "4", "9", "14", "20") ~ "T",
  c("5", "6", "7", "15")     ~ "Myeloid",
  "0"                         ~ "B",
  c("8", "11")                ~ "Fibroblast",
  "10"                        ~ "Plasma",
  "18"                        ~ "Endothelial",  # also shows Adipocyte markers
  "17"                        ~ "Pericyte",
  # less confident
  c("13", "16", "19")         ~ "Epithelial"
))

DimPlot(seu, group.by = "cell_type0.7", cols = cell_type_colors, label = TRUE)

## Save

seu$cell_type202605 <- seu$cell_type0.7
saveRDS(seu, "/vast/projects/Spatial/lei/Benchmarking/scRNA/TNBC_04/seu_TNBC_04.rds")
