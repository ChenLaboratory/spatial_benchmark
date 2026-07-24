# Modified version of Andy's script to create model matrices from specified samples into tsv.gz file

library(Seurat)
library(edgeR)
SampleIDs <- c("MH0007",
               "MH0026")
levels <- c("Adipocyte","B","Endothelial","Epithelial","Fibroblast",
            "Myeloid","Pericyte","Plasma","T","Tumor")
for (SampleID in SampleIDs) {
  so <- readRDS(paste0("/vast/projects/Spatial/lei/Benchmarking/Xenium/",SampleID,"/so.rds"))
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
  out_path <- paste0("/vast/ai_projects/IBD/nguyen.q/Pilot_Human_Breast_2_withref/Ficture_Xenium/predefine_model_m/model_matrix/23",SampleID)
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