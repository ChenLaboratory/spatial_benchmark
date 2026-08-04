# Purpose:  Infer large-scale copy-number alterations in one sc/snRNA-seq sample with
#           inferCNV, using a normal human breast reference as the diploid baseline.
#           Produces the CNV heatmap and the per-cluster genomic-instability score used
#           to decide which clusters are tumour in 03_annotate_cell_types.R.
#           Run once per sample.
# Inputs:   args[1]   path to the sample's seu.rds (clustered object from 01_reprocess.R)
#           args[2]   sample name, used in the heatmap title
#           seu_normal_human_breast.rds   normal breast scRNA-seq reference
#           hg38_gencode_v27.txt          gene order file (GENCODE v27 / hg38)
#           The last two are not covered by the GEO accession, and their paths are
#           hardcoded below — edit them before running.
# Outputs:  inferCNV_res/ in the working directory:
#             run.final.infercnv_obj                     inferCNV's own output object
#             ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png
#             heatmap_list_combined.rds, heatmap_list_raw.rds
#             instability_score_from_infercnv.rds        per-cluster instability score
#             Boxplot_infercnv_instability_score.png
# Usage:    Rscript 02_inferCNV.R <path/to/seu.rds> <sample_name>

library(infercnv)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)

mycolor <-  c(brewer.pal(8, "Set2"), brewer.pal(8, "Set1"), brewer.pal(8, "Set3"))

args = commandArgs(trailingOnly=TRUE)
print(args)

obs_seu_dir = args[1]
obs_sample = args[2]

ref_seu_dir = "/vast/projects/lab_chen/lei/Data/homo_sapiens/scRNA_normal_human_breast/seu_normal_human_breast.rds"
ref_sample = "N0372" 
gene_order_file_dir="/vast/projects/lab_chen/lei/Data/homo_sapiens/hg38_gencode_v27.txt"
infercnv_cutoff = 0.1
out_dir = "inferCNV_res"

print("Preparing input for inferCNV") #----------------------------------------

# read in normal cell reference
ref_seu <- readRDS(ref_seu_dir)
ref_seu$group <- "Normal"

# read in seurat object
obs_seu <- readRDS(obs_seu_dir)
obs_seu$group <- obs_seu$seurat_clusters

# combine ref and obs
seu <- merge(ref_seu, y = obs_seu, add.cell.ids = c("ref", "obs"))
seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
seu$sample <- sapply(strsplit(colnames(seu),"_"), `[`, 1)

seu$group <- factor(seu$group, levels = c(levels(obs_seu$group), unique(ref_seu$group)))

# specify result path
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# prepare input for inferCNV
counts_matrix = as.matrix(seu[["RNA"]]$counts[, colnames(seu)]) 

anno <- seu@meta.data[,"group", drop=FALSE]
anno <- as.matrix(anno)
colnames(anno) <- NULL

print("Running inferCNV") #----------------------------------------

# create inferCNV object
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=counts_matrix,
  annotations_file=anno,
  delim="\t",
  gene_order_file = gene_order_file_dir,
  ref_group_names = "Normal")

# run inferCNV
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=infercnv_cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out_dir,
  cluster_by_groups=TRUE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  output_format = "pdf"
)

res_title <- paste0("inferCNV results for sample ", obs_sample, " using normal reference ", ref_sample)

print("generating manual heatmap") #----------------------------------------

filename_infercnv_obj = paste0(out_dir,"/","run.final.infercnv_obj")
data_for_plot_cnv = readRDS(filename_infercnv_obj)

ordered_groups = do.call(c, list(data_for_plot_cnv@observation_grouped_cell_indices,
                                 data_for_plot_cnv@reference_grouped_cell_indices) )

ordered_groups = cbind.data.frame( group = rep(names(ordered_groups),
                                               c(lapply(ordered_groups, length))),
                                   position = do.call(rbind,
                                                      lapply(ordered_groups, data.frame) )[,1] )

ordered_groups$group = factor(ordered_groups$group, levels = levels(seu$group) ) # use ordered levels in the raw data

data_for_heatmap = data_for_plot_cnv@expr.data[, ordered_groups$position]

ordered_groups$sample = seu@meta.data[colnames(data_for_heatmap),"sample"]# add sample information

# colors for group used in infercnv
my_colors_col = mycolor
color_groups_col = my_colors_col[1:nlevels(ordered_groups$group)]
names(color_groups_col) = levels(ordered_groups$group)

# colors for subsetting data
ordered_groups$group = droplevels(ordered_groups$group)
color_groups_col = color_groups_col[levels(ordered_groups$group)]

# colors for sample information from dge data
my_colors_col = mycolor
color_samples_col = my_colors_col[1:length(unique(ordered_groups$sample))]
names(color_samples_col) = unique(ordered_groups$sample)

list_data_for_heatmap = list()
list_ordered_groups = list()
list_column_ha = list()
list_ht_raw = list()

for (i in levels(ordered_groups$group) ) {
  keep_columns <- (ordered_groups$group==i)
  list_data_for_heatmap[[i]] = data_for_heatmap[,keep_columns]
  list_ordered_groups[[i]] = ordered_groups[ordered_groups$group==i,]
  
  list_column_ha[[i]] = HeatmapAnnotation( group  = list_ordered_groups[[i]]$group,
                                           sample = list_ordered_groups[[i]]$sample,
                                           col = list(group  = color_groups_col,
                                                      sample = color_samples_col),
                                           show_legend = TRUE, # show legend
                                           show_annotation_name = FALSE )
  
  ht =  Heatmap(matrix = list_data_for_heatmap[[i]],
                
                width = 7,
                
                # column anno #
                top_annotation =  list_column_ha[[i]],
                cluster_columns = TRUE,
                show_column_dend = FALSE, # show denrogram or not
                show_parent_dend_line = FALSE, # remove dashed line of parent dendrogram
                show_heatmap_legend  = FALSE,
                column_split = list_ordered_groups[[i]]$sample,
                column_title = i,    # show column title or not (level names in `column_split`)
                column_gap = unit(0, "mm"), # column split gap
                column_title_rot = 1,
                
                # row anno #
                cluster_rows = FALSE,
                row_split = data_for_plot_cnv@gene_order$chr,
                # row_title = NULL,      # show row title or not (level names in `row_split`)
                row_title_rot = 0,
                
                use_raster = TRUE,
                
                show_column_names = FALSE,
                show_row_names = FALSE)
  
  
  list_ht_raw[[i]] = ht
  
  if(i==levels(ordered_groups$group)[1]){
    list_ht = ht
  }else{
    list_ht = list_ht + ht
  }
  
}

output_filename = paste0(out_dir,"/","ComplexHeatmap_infercnv_with_sampleInfo_equal_width.png")

date()

png(filename = output_filename,
    width = 8,height = 10,units = "in",res = 300)


ComplexHeatmap::draw(list_ht,ht_gap = unit(1, "mm"), column_title = res_title)

dev.off()

date()

filename_rds = paste0(out_dir,"/","heatmap_list_combined.rds")
saveRDS(list_ht,file = filename_rds)

filename_rds = paste0(out_dir,"/","heatmap_list_raw.rds")
saveRDS(list_ht_raw,file = filename_rds)

print("generating instability score") #----------------------------------------

# calculate instability score
instability_score = apply(( data_for_heatmap-1)**2, 2, mean)

saveRDS(instability_score, file = paste0(out_dir,"/","instability_score_from_infercnv.rds") )

# boxplot 
data_for_plot = cbind.data.frame(instability_score, ordered_groups)

my_colors =  c(mycolor, 
               rev(viridis::plasma(10)) ) 
p1 <- ggplot(data = data_for_plot,
             aes(x=group, 
                 y=instability_score, color=group ) ) + 
  geom_boxplot(show.legend=FALSE) + 
  scale_color_manual(values=my_colors ) + 
  labs(x="", y="inferCNV instability score",title = res_title)+ 
  Seurat::RotatedAxis() +
  # coord_cartesian(ylim=c(0,0.02) ) +
  theme_classic()+ theme(plot.title = element_text(size=8))


n_scale = nlevels(ordered_groups$group)/10

my_output_filename = paste0(out_dir,"/", "Boxplot_infercnv_instability_score.png")

png(my_output_filename,
    height = 4*1,
    width = 4*1*n_scale,
    units = "in",res = 300)
p1
dev.off()



