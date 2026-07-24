# h5ad -> Seurat object, for one sample.
# Usage: Rscript 03_generate_seu.R
#
# Note: running this as a batch Rscript job previously hit a readH5AD error of
# unclear cause. The identical logic below is confirmed working when run
# interactively (e.g. as an .Rmd chunk) -- if the batch error recurs, fall back
# to running this interactively instead of via sbatch.

library(Seurat)
library(zellkonverter)
library(dplyr)

sce <- readH5AD("proseg-anndata.h5ad")
seu <- as.Seurat(sce, counts = "X", data = NULL)
seu <- RenameAssays(seu, originalexp = "RNA")

saveRDS(seu, "so_raw.rds")
message("Saved: so_raw.rds")
