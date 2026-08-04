# Purpose:  Combine Visium with Xenium/MERSCOPE, using the STalign + hexbin alignment
#           from 01 and 02, into per-sample cross-platform Seurat object lists — one
#           list of all cells/spots, one restricted to the region the platforms overlap.
#           Attaches alignment metrics to each Visium spot (binned iST counts, per hexbin)
#           and to each iST cell (its Visium spot ID and aligned coordinates), and derives
#           the shared gene panels. Everything in ../analysis that compares platforms
#           reads these objects.
# Inputs:   <processed_data_dir>/Visium/<sample>/so_raw.rds, seu_ST_Filtered.rds
#           <processed_data_dir>/<platform>/<sample>/so_raw.rds, so.rds
#           <processed_data_dir>/integration/<sample>_<platform>_hexbin.csv.gz
#                                                       hexbin assignment from 02_hexbin.py
#           raw_gene_panel_file                         panel definition, both platforms
# Outputs:  <processed_data_dir>/shared_genes_iST.rds
#                     genes shared by Xenium and MERSCOPE
#           <processed_data_dir>/shared_genes_withVisium.rds
#                     the 213 genes shared by all three platforms
#           <processed_data_dir>/raw_data_list_AllSample.rds
#                     per-sample, per-platform raw objects with alignment metadata
#           <processed_data_dir>/overlapped_data_list_AllSample.rds
#                     the same, restricted to the overlapped region
#           See README.md for the structure of the two object lists.

library(Seurat)
library(tidyverse)
library(dplyr)
library(magrittr)
library(data.table)

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))

# note: `samples` in config.R is the long-form (sample, platform) manifest; the
# loop below iterates over sample IDs, so use bench_samples directly
samples_with_MERSCOPE <- bench_samples_merscope

determine_platforms <- function(x) {
  if (x %in% samples_with_MERSCOPE) {
    return(c("Xenium", "MERSCOPE"))
  } else {
    return("Xenium")
  }
}

visium_dir      <- file.path(processed_data_dir, "Visium")        # <sample>/so_raw.rds, <sample>/seu_ST_Filtered.rds
ist_dir         <- processed_data_dir            # <platform>/<sample>/so_raw.rds, so.rds
hexbin_dir      <- file.path(processed_data_dir, "integration")  # <sample>_<platform>_hexbin.csv.gz
gene_panel_file <- raw_gene_panel_file

## Gene panels
# 213 shared genes across Xenium, MERSCOPE and Visium

iST_gene_panel <- readxl::read_excel(gene_panel_file)
unique(iST_gene_panel$...5[-1]) # MERSCOPE
unique(iST_gene_panel$...12[-1]) # Xenium

# shared gene panel between MERSCOPE and Xenium
shared_genes <- intersect(iST_gene_panel$MERSCOPE[-1], iST_gene_panel$Xenium[-1]) # Xenium & MERSCOPE
saveRDS(shared_genes, file.path(processed_data_dir, "shared_genes_iST.rds"))
iST_gene_panel$...5[iST_gene_panel$MERSCOPE %in% shared_genes] <- iST_gene_panel$...12[iST_gene_panel$MERSCOPE %in% shared_genes] # refer to xenium annotation for the shared genes

# shared gene panel between MERSCOPE and Xenium and Visium
visium_raw_seu <- readRDS(file.path(visium_dir, "TNBC_02", "so_raw.rds"))
shared_genes <- intersect(shared_genes, rownames(visium_raw_seu)) # iST & Visium
saveRDS(shared_genes, file.path(processed_data_dir, "shared_genes_withVisium.rds"))

## Store alignment results to seu
# Visium: iST_total_nCount_hexbin, iST_avg_nCount_hexbin (per cell),
#         iST_total_nCount_hexbin_shared_genes, iST_avg_nCount_hexbin_shared_genes
# iST:    Visium_spot_id, aligned spatial map

calculate_CellNumber_hexbin <- function(hexbin, visium.seu, sample, platform, name = "cell_count") {
  col_name = paste0(platform, "_", name)

  print(paste0("There are ", sum(hexbin$Visium_spot_id == ""), " out of ", nrow(hexbin), " cells not assigned to hexbins. (", sum(hexbin$Visium_spot_id == "")/nrow(hexbin), ")"))

  counts <- as.data.frame(table(hexbin$Visium_spot_id))
  colnames(counts) <- c("Visium_spot_id", col_name)
  counts %<>% filter(Visium_spot_id != "")

  visium.seu@meta.data <- visium.seu@meta.data %>%
      rownames_to_column(var = "Visium_spot_id")%>%
      left_join(counts, by = "Visium_spot_id") %>%
      column_to_rownames("Visium_spot_id")

  print(paste0("There are ", sum(is.na(visium.seu@meta.data[[col_name]])), " empty spot/hexbins out of ",nrow(visium.seu@meta.data), " in-tissue spots. (", sum(is.na(visium.seu@meta.data[[col_name]]))/nrow(visium.seu@meta.data), ")"))

  return(visium.seu)
}

calculate_nCount_hexbin <- function(iST.seu, visium.seu, platform, name = "cell_count"){
  col_name = paste0(platform, "_", name)

  # sum nCount_RNA of all iST cells in the same hexbin
  ncount_sum <- iST.seu@meta.data %>%
    group_by(Visium_spot_id) %>%
    summarise(total_nCount_RNA = sum(nCount_RNA, na.rm = TRUE))

  rownames(ncount_sum) <- ncount_sum$Visium_spot_id # total read of all iST cells inside each hexbin

  # add sum parameter to visium metadata
  visium.seu@meta.data <- visium.seu@meta.data %>%
    rownames_to_column(var = "Visium_spot_id")%>%
    left_join(ncount_sum, by = "Visium_spot_id") %>%
    column_to_rownames("Visium_spot_id")

  name <- paste0(platform, "_total_nCount_hexbin")
  colnames(visium.seu@meta.data) <- gsub("total_nCount_RNA", name, colnames(visium.seu@meta.data))

  # calculate avg Count per iST cell in each hexbin
  avg_name <- paste0(platform, "_avg_nCount_hexbin")
  visium.seu@meta.data[,avg_name] <- visium.seu@meta.data[, name]/visium.seu@meta.data[, col_name]

  return(visium.seu)
}

calculate_count_hexbin_shared_genes <- function(iST, Visium, platform, name = "cell_count"){
  col_name = paste0(platform, "_", name)

  # new nCount on shared genes
  iST@meta.data$nCount_RNA <- colSums(iST[shared_genes,]@assays$RNA@layers$counts)

  # sum nCount_RNA of all iST cells in the same hexbin
  ncount_sum <- iST@meta.data %>%
    group_by(Visium_spot_id) %>%
    summarise(total_nCount_RNA = sum(nCount_RNA, na.rm = TRUE))

  rownames(ncount_sum) <- ncount_sum$Visium_spot_id # total read of all iST cells inside each hexbin

  # add sum parameter to visium metadata
  Visium@meta.data <- Visium@meta.data %>%
    rownames_to_column(var = "Visium_spot_id")%>%
    left_join(ncount_sum, by = "Visium_spot_id") %>%
    column_to_rownames("Visium_spot_id")

  name <- paste0(platform, "_total_nCount_hexbin_shared_genes")
  colnames(Visium@meta.data) <- gsub("total_nCount_RNA", name, colnames(Visium@meta.data))

  # calculate avg Count per iST cell in each hexbin
  avg_name <- paste0(platform, "_avg_nCount_hexbin_shared_genes")
  Visium@meta.data[,avg_name] <- Visium@meta.data[, name]/Visium@meta.data[, col_name]

  return(Visium)
}

### raw data
# Process raw data without QC, adding alignment details.

raw_data_list_AllSample <- vector("list", length(bench_samples))
names(raw_data_list_AllSample) <- bench_samples

for (sample in bench_samples) {
  platforms <- determine_platforms(sample)

  raw_data_list <- vector("list", 3)
  names(raw_data_list) <- c("Visium", "Xenium", "MERSCOPE")

  # visium
  print(paste0("============",sample, " - ","Visium", "============"))
  visium_raw_seu <- readRDS(file.path(visium_dir, sample, "so_raw.rds"))

  # iST
  for (platform in platforms) {
    print(paste0("============",sample, " - ",platform, "============"))

    # alignment results
    hexbin <- fread(file.path(hexbin_dir, paste0(sample, "_", platform, "_hexbin.csv.gz"))) %>% as.data.frame()

    # iST seu
    iST_raw_seu <- readRDS(file.path(ist_dir, platform, sample, "so_raw.rds"))

    if (platform == "Xenium") {
      rownames(hexbin) <- hexbin$cell_id
      iST_id <- "cell_id"
    } else {
      rownames(hexbin) <- hexbin$EntityID
      iST_id <- "EntityID"
    }

    # add QC metrics
    visium_raw_seu <- calculate_CellNumber_hexbin(hexbin = hexbin, visium_raw_seu, sample, platform) # calculate the number of cells in each hexbin corresponding to visium spot id => "Xenium_cell_count"

    iST_raw_seu@meta.data$Visium_spot_id <- hexbin[colnames(iST_raw_seu), "Visium_spot_id"] # hexbin id for each iST cell

    visium_raw_seu <- calculate_nCount_hexbin(iST_raw_seu, visium_raw_seu, platform) # calculate the total counts detected in all cells in each hexbin (separate panels) => "Xenium_total_nCount_hexbin" & "Xenium_avg_nCount_hexbin"

    visium_raw_seu <- calculate_count_hexbin_shared_genes(iST_raw_seu[shared_genes,], visium_raw_seu, platform) # calculate the total counts detected in all cells in each hexbin (shared panels) => "Xenium_total_nCount_hexbin_shared_genes" & "Xenium_avg_nCount_hexbin_shared_genes"

    visium_raw_seu$nCount_Spatial_shared_genes <- colSums(visium_raw_seu[shared_genes,]@assays$Spatial@layers$counts) # calculate the total counts for shared panel in each visium spot => "nCount_Spatial_shared_genes"

    # add aligned spatial map
    iST_raw_seu@meta.data$aligned_x <- hexbin$aligned_x
    iST_raw_seu@meta.data$aligned_y_reversed <- -hexbin$aligned_y # original coordinates follows image coord system which starts from left top corner

    sp_embeddings <- iST_raw_seu@meta.data[, c("aligned_x", "aligned_y_reversed")]
    colnames(sp_embeddings) <- c("aligned_Sp_1", "aligned_Sp_2")
    spatial <- CreateDimReducObject(embeddings = as.matrix(sp_embeddings), assay = "RNA", key = "aligned_Sp_")
    iST_raw_seu@reductions[["aligned_spatial"]] <- spatial

    raw_data_list[[platform]] <- iST_raw_seu
  }

  visium_raw_seu@meta.data[is.na(visium_raw_seu@meta.data)] <- 0
  raw_data_list$Visium <- visium_raw_seu

  raw_data_list_AllSample[[sample]] <- raw_data_list
}

saveRDS(raw_data_list_AllSample, file.path(processed_data_dir, "raw_data_list_AllSample.rds"))

### processed data

processed_data_list_AllSample <- vector("list", length(bench_samples))
names(processed_data_list_AllSample) <- bench_samples

for (sample in bench_samples) {
  platforms <- determine_platforms(sample)

  processed_data_list <- vector("list", 3)
  names(processed_data_list) <- c("Visium", "Xenium", "MERSCOPE")

  print(paste0("============",sample, " - ","Visium", "============"))
  visium_processed_seu <- readRDS(file.path(visium_dir, sample, "seu_ST_Filtered.rds"))

  for (platform in platforms) {
    print(paste0("============",sample, " - ",platform, "============"))

    hexbin <- fread(file.path(hexbin_dir, paste0(sample, "_", platform, "_hexbin.csv.gz"))) %>% as.data.frame()
    iST_processed_seu <- readRDS(file.path(ist_dir, platform, sample, "so.rds"))

    if (platform == "Xenium") {
      rownames(hexbin) <- hexbin$cell_id
      iST_id <- "cell_id"
    } else {
      rownames(hexbin) <- hexbin$EntityID
      iST_id <- "EntityID"
    }

    # remove hexbins' spot filtered out in visium processed data
    hexbin_filtered <- hexbin[colnames(iST_processed_seu),] %>%
      mutate(Visium_spot_id = case_when(!Visium_spot_id %in% colnames(visium_processed_seu) ~ "",
                                        .default = Visium_spot_id)) %>%
      group_by(Visium_spot_id) %>%
      mutate(iST_cell_count_filtered  = ifelse(Visium_spot_id == "", NA, n())) %>%
      ungroup() %>% column_to_rownames(iST_id) # re-count hexbin cell numbers based on filtered data

    # add QC metrics
    iST_processed_seu@meta.data$Visium_spot_id <- hexbin_filtered[colnames(iST_processed_seu), "Visium_spot_id"] # hexbin id for each iST cell

    visium_processed_seu <- calculate_CellNumber_hexbin(hexbin = hexbin_filtered, visium_processed_seu, sample, platform, name = "cell_count_filtered") # calculate the number of cells in each hexbin corresponding to visium spot id => "Xenium_cell_count_filtered"

    visium_processed_seu <- calculate_nCount_hexbin(iST_processed_seu, visium_processed_seu, platform, name = "cell_count_filtered") # calculate the total counts detected in all cells in each hexbin (separate panels) => "Xenium_total_nCount_hexbin" & "Xenium_avg_nCount_hexbin"

    visium_processed_seu <- calculate_count_hexbin_shared_genes(iST_processed_seu[shared_genes,], visium_processed_seu, platform, name = "cell_count_filtered") # calculate the total counts detected in all cells in each hexbin (shared panels) => "Xenium_total_nCount_hexbin_shared_genes" & "Xenium_avg_nCount_hexbin_shared_genes"

    visium_processed_seu$nCount_Spatial_shared_genes <- colSums(visium_processed_seu[shared_genes,]@assays$Spatial@layers$counts) # calculate the total counts for shared panel in each visium spot => "nCount_Spatial_shared_genes"

    processed_data_list[[platform]] <- iST_processed_seu
  }
  visium_processed_seu@meta.data[is.na(visium_processed_seu@meta.data)] <- 0
  processed_data_list$Visium <- visium_processed_seu

  processed_data_list_AllSample[[sample]] <- processed_data_list
}

## Get overlapped region
# Only use overlapped hexbins for downstream comparison -- i.e. hexbins with
# cells from both Xenium & MERSCOPE (where both platforms exist for a sample).

add_aligned_map <- function(seu, sample, platform) {
  hexbin <- fread(file.path(hexbin_dir, paste0(sample, "_", platform, "_hexbin.csv.gz"))) %>% as.data.frame()
  if (platform == "Xenium") {
    rownames(hexbin) <- hexbin$cell_id
    iST_id <- "cell_id"
  } else {
    rownames(hexbin) <- hexbin$EntityID
    iST_id <- "EntityID"
  }

  hexbin_filtered <- hexbin[colnames(seu),]
  seu$aligned_x <- hexbin_filtered$aligned_x
  seu$aligned_y <- hexbin_filtered$aligned_y
  seu$aligned_x_in_fullres <- hexbin_filtered$aligned_x_in_fullres
  seu$aligned_y_in_fullres <- hexbin_filtered$aligned_y_in_fullres

  # add spatial map
  sp_embeddings <- seu@meta.data[, c("aligned_x", "aligned_y")]
  sp_embeddings$aligned_y <- -sp_embeddings$aligned_y
  colnames(sp_embeddings) <- c("Aligned_1", "Aligned_2")
  spatial <- CreateDimReducObject(embeddings = as.matrix(sp_embeddings), assay = "RNA", key = "Aligned_")
  seu@reductions[["aligned"]] <- spatial

  return(seu)
}

overlapped_data_list_AllSample <- vector("list", length(bench_samples))
names(overlapped_data_list_AllSample) <- bench_samples

for (sample in bench_samples) {

  overlapped_data_list <- vector("list", 3)
  names(overlapped_data_list) <- c("Visium", "Xenium", "MERSCOPE")

  platforms <- determine_platforms(sample)

  print(paste0("============",sample, " - ","Visium", "============"))

  visium <- processed_data_list_AllSample[[sample]]$Visium

  if ("MERSCOPE" %in% platforms) {
    visium_overlapped <- visium[, visium$Xenium_cell_count_filtered != 0 & visium$MERSCOPE_cell_count_filtered != 0]
  } else {
    visium_overlapped <- visium[, visium$Xenium_cell_count_filtered != 0]
  }

  overlapped_data_list[["Visium"]] <- visium_overlapped

  for (platform in platforms) {
    print(paste0("============",sample, " - ",platform, "============"))
    iST <- processed_data_list_AllSample[[sample]][[platform]]

    # overlapping info
    iST_overlapped <- iST[, iST$Visium_spot_id %in% rownames(visium_overlapped@meta.data)]

    iST_overlapped <- add_aligned_map(iST_overlapped, sample = sample, platform = platform)

    overlapped_data_list[[platform]] <- iST_overlapped
  }

  overlapped_data_list_AllSample[[sample]] <- overlapped_data_list
}

saveRDS(overlapped_data_list_AllSample, file.path(processed_data_dir, "overlapped_data_list_AllSample.rds"))