# Purpose:  Convert the raw Xenium transcripts.parquet files to the CSV form FICTURE's
#           Snakemake workflow expects. Drops rows missing transcript ID, feature name
#           or coordinates, and removes control probes (BLANK, NegControlCodeword).
# Inputs:   <xenium_path>/transcripts.parquet   one per sample; paths are hardcoded
#                                               below and must be edited before running
# Outputs:  <sample_name>/<sample_name>.csv     one directory per sample, written under
#                                               the working directory set by setwd() below
#           The "23" sample-name prefix is a legacy naming artefact.

library(arrow)

# Output directory
setwd("/vast/ai_projects/IBD/nguyen.q/Pilot_Human_Breast_2/Ficture_Xenium/transcripts_qc")

xenium_path <- c(
  "/vast/projects/lab_chen/Dataset/Xenium_Pilot_Human_Breast_2/output-XETG00068__0003647__Region_2__20230727__025552/",
  "/vast/projects/lab_chen/Dataset/Xenium_Pilot_Human_Breast_2/output-XETG00068__0003647__Region_3__20230727__025552/"
)
sample_name <- c(
  "23TNBC_01",
  "23TNBC_02"
)
transcript_path = paste0(xenium_path,"transcripts.parquet")

for (i in seq_along(transcript_path)) {
  transcript <- as.data.frame(arrow::read_parquet(transcript_path[i]))
  # Grab out sample name. regex made by chatgpt
  s_n <- sample_name[i]
  if (!file.exists(s_n)) {
    dir.create(s_n)
  }
  
  # Remove rows with NA in certain columns
  transcript = transcript[complete.cases(transcript[,c("transcript_id","feature_name","x_location","y_location")]),]
  
  # Remove control probes
  kp <- !grepl("(BLANK)|(NegControlCodeword)",
              transcript$feature_name)
  transcript = transcript[kp,]
  
  write.csv(transcript, file = paste0(s_n,"/",s_n,".csv"),
              row.names = FALSE)
}