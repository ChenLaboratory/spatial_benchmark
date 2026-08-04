# Purpose:  Convert one sample's Proseg v3 output from h5ad to a raw Seurat object.
#           Run from the sample's Proseg output directory — both paths below are
#           relative to the working directory.
# Inputs:   proseg-anndata.h5ad   Proseg output, converted by zarr_to_h5ad.py
# Outputs:  so_raw.rds            raw object, taken up by ../03_postprocess_seu.R
# Usage:    Rscript generate_seu.R
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
