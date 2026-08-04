# Purpose:  Build the fixed reference model matrix that guides FICTURE's factor
#           decomposition of the MERSCOPE data, in place of FICTURE's default
#           unsupervised LDA. Pseudobulks each annotated cell type
#           (SingleR_labels202605) with edgeR, converts to log CPM, and writes the
#           matrix and its colour scheme in the format FICTURE expects. Cell types
#           absent from a sample are dropped, along with their colours.
# Inputs:   <MERSCOPE processed dir>/<sample>/so.rds   annotated MERSCOPE objects
#           celltype_all.rgb.tsv                       cell type colour scheme
#           Both paths are hardcoded below — edit them before running.
# Outputs:  <out_dir>/23<sample>/model_matrix.tsv.gz   reference model for FICTURE
#           <out_dir>/23<sample>/celltype.rgb.tsv      matching colour scheme
#           The "23" filename prefix is a legacy sample-naming artefact.
#
# Modified version of Andy's script to create model matrices from specified samples into tsv.gz file

library(Seurat)
library(edgeR)
SampleIDs <- c("TNBC_01",
               "TNBC_02")
levels <- c("Adipocyte","B","Endothelial","Epithelial","Fibroblast",
            "Myeloid","Pericyte","Plasma","T","Tumor")
for (SampleID in SampleIDs) {
  so <- readRDS(paste0("/vast/projects/Spatial/lei/Benchmarking/MERSCOPE/",SampleID,"/so.rds"))
  Idents(so) <-factor(Idents(so),
                               levels=levels)
  
  y  <- Seurat2PB(so, sample="orig.ident",cluster = "SingleR_labels202605")
  colnames(y) <- gsub("^.*cluster", "", colnames(y))
  CPM <- edgeR::cpm(y, log=FALSE)
  mat <- log1p(CPM/1e2)

  # Some sample may not have all cell types. Removing those cell types & colors
  m <- match(levels, colnames(mat))
  mat <- mat[, m]
  keep <- !is.na(colnames(mat))
  mat <- mat[,keep]
  # Color
  rgb <- read.csv("/vast/ai_projects/IBD/nguyen.q/Pilot_Human_Breast_2_withref/celltype_all.rgb.tsv",
                  sep="\t")
  rgb <- rgb[keep,]
  rgb$Name <- 0:(nrow(rgb)-1)
  rgb$Color_index <- rgb$Name
  
  # Write model_matrix. Append 23 to sampleID so I dont have to change previous code
  out_path <- paste0("/vast/ai_projects/IBD/nguyen.q/Pilot_Human_Breast_2_withref/Ficture_Merscope/predefine_model_m/model_matrix/23",SampleID)
  dir.create(out_path,recursive = TRUE)
  mat <- data.frame(gene = rownames(mat), mat)
  write.table(mat, file = gzfile(paste0(out_path,"/model_matrix.tsv.gz")),
              sep = "\t",
              row.names = FALSE)
  # Write color scheme
  write.table(rgb, file = paste0(out_path,"/celltype.rgb.tsv"),
              sep = "\t",
              row.names = FALSE)
}