# Purpose:  Proportion of Visium spots/bins with a significant BLISA ligand-receptor
#           interaction, broken down by pathologist-annotated tissue region, LR pair and
#           platform — i.e. where in the tissue each platform detects each interaction.
#           Spots in the "Unannotated" region are excluded. See README.md.
# Inputs:   <results_table_dir>/seu_v_annotated.rds                  from 01
#                     supplies the per-spot pathology_region labels
#           <results_table_dir>/BLISA_x.rds, BLISA_m.rds, BLISA_v.rds   from 02
# Outputs:  <results_table_dir>/prop_all.rds
#                     feeds figures/Fig5_celltype_cci.R panel d

# ---- Paths: edit config.R at the repository root, not this file ----
source(here::here("config.R"))
seu_v   <- readRDS(file.path(results_table_dir, "seu_v_annotated.rds"))
BLISA_x <- readRDS(file.path(results_table_dir, "BLISA_x.rds"))
BLISA_m <- readRDS(file.path(results_table_dir, "BLISA_m.rds"))
BLISA_v <- readRDS(file.path(results_table_dir, "BLISA_v.rds"))

spot_regions <- seu_v$pathology_region

all_lr      <- unique(c(rownames(BLISA_x$LR_out), rownames(BLISA_m$LR_out), rownames(BLISA_v$LR_out)))
all_regions <- setdiff(unique(spot_regions), "Unannotated")
total_per_region <- table(spot_regions)

calc_region_proportion <- function(LR_out, platform) {
  do.call(rbind, lapply(all_lr, function(lr) {
    if (!lr %in% rownames(LR_out)) {
      return(data.frame(platform = platform, LR_pair = lr, region = all_regions,
                         n_total = NA_integer_, n_sig = NA_integer_, stringsAsFactors = FALSE))
    }
    sig_idx <- LR_out[lr, "sig_index"][[1]]
    if (is.null(sig_idx) || length(sig_idx) == 0) {
      return(data.frame(platform = platform, LR_pair = lr, region = all_regions,
                         n_total = as.integer(total_per_region[all_regions]),
                         n_sig = 0L, stringsAsFactors = FALSE))
    }
    sig_regions    <- spot_regions[sig_idx]
    sig_per_region <- table(sig_regions)
    n_sig_vals     <- as.integer(sig_per_region[all_regions])
    n_sig_vals[is.na(n_sig_vals)] <- 0L
    data.frame(platform = platform, LR_pair = lr, region = all_regions,
               n_total = as.integer(total_per_region[all_regions]), n_sig = n_sig_vals,
               stringsAsFactors = FALSE)
  }))
}

prop_x <- calc_region_proportion(BLISA_x$LR_out, "Xenium")
prop_m <- calc_region_proportion(BLISA_m$LR_out, "MERSCOPE")
prop_v <- calc_region_proportion(BLISA_v$LR_out, "Visium")

prop_all <- rbind(prop_x, prop_m, prop_v)
prop_all$proportion <- ifelse(is.na(prop_all$n_total), NA_real_, prop_all$n_sig / prop_all$n_total)
prop_all$platform <- factor(prop_all$platform, levels = c("Xenium", "MERSCOPE", "Visium"))

saveRDS(prop_all, file.path(results_table_dir, "prop_all.rds"))
