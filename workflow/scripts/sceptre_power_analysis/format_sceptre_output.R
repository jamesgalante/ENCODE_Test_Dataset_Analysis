# Script: format_sceptre_output.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/format_sceptre_output")) { dir.create("RDA_objects/format_sceptre_output", recursive = TRUE) }
save.image(paste0("RDA_objects/format_sceptre_output/", snakemake@wildcards$sample, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
# Read power analysis results (multiple TSV files)
power_analysis_data <- lapply(snakemake@input$power_analysis_results, function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Extract effect size from filename
  list(effect_size = as.numeric(sub(".*es_(.+)\\.tsv", "\\1", file)), data = data)
})

# Read in the other input files
discovery_results <- readRDS(snakemake@input$discovery_results)
gene_gRNA_pairs <- readRDS(snakemake@input$gene_gRNA_group_pairs)
distances <- read.table(snakemake@input$distances, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
guide_targets <- read_tsv(snakemake@input$guide_targets)

# Remove any duplicates in "distances" - there shouldn't be any
if(any(duplicated(distances[, c("response_id", "grna_group")]))) {
  message("Found duplicates in distances file - removing them")
  distances <- distances[!duplicated(distances[, c("response_id", "grna_group")]),]
}

### MERGE FILES TOGETHER =====================================================

# Basically first merge the discovery_results, gene_gRNA_pairs, distances because the rows are all the same
merged_data <- discovery_results %>%
  left_join(gene_gRNA_pairs, by = c("response_id" = "response_id", "grna_target" = "grna_target")) %>%
  left_join(distances, by = c("response_id" = "response_id", "grna_target" = "grna_group"))


### COMBINE THE POWER ANALYSES ===============================================

combine_power_analysis <- function(power_analysis_data) {
  # Get all effect sizes and create column names
  effect_sizes <- sapply(power_analysis_data, function(x) x$effect_size)
  power_colnames <- paste0("PowerAtEffectSize", effect_sizes * 100)
  
  # Process each dataset and merge
  power_tables <- lapply(seq_along(power_analysis_data), function(i) {
    df <- power_analysis_data[[i]]$data
    # Rename power column using corresponding name from our vector
    names(df)[names(df) == "power"] <- power_colnames[i]
    # Keep only needed columns
    df[, c("grna_target", "response_id", power_colnames[i])]
  })
  
  # Merge all tables together
  combined_power <- Reduce(function(x, y) {
    merge(x, y, by = c("grna_target", "response_id"), all = TRUE)
  }, power_tables)
  
  # Order columns: first the IDs, then power columns in order of effect size
  col_order <- c("grna_target", "response_id", sort(power_colnames))
  combined_power <- combined_power[, col_order]
  
  return(combined_power)
}

# Use it like this:
combined_power_analysis_data <- combine_power_analysis(power_analysis_data)


### CONVERT POSITIVE CONTROLS TO TSSCtrl =====================================

# Convert specific target types to TSSCtrl
merged_data$target_type <- ifelse(
  merged_data$target_type %in% c("tss_pos", "tss_random"), 
  "TSSCtrl", 
  merged_data$target_type
)


### CREATE FINAL OUTPUT ======================================================

# Start with combined_power_analysis_data and merge in the merged_data info
final_merged <- merge(combined_power_analysis_data, 
                      merged_data,
                      by = c("grna_target", "response_id"),
                      all.x = TRUE)

# Get all unique guide targets (for smooth merging with final_merged)
unique_guide_targets <- guide_targets %>% 
  group_by(target_name) %>%
  dplyr::slice(1)

# Combine final_merged with guide_targets to get pert_start/end for TSSCtrls
final_merged <- final_merged %>%
  left_join(unique_guide_targets, by = c("grna_target" = "target_name"))

# Create the final output table
final_output <- data.frame(
  perturbation = final_merged$grna_target,
  gene = final_merged$response_id,
  logFC = final_merged$log_2_fold_change,
  ci_high = NA,
  ci_low = NA,
  pvalue = final_merged$p_value,
  pval_adj = p.adjust(final_merged$p_value, method = "BH"),
  pert_chr = final_merged$target_chr,
  pert_start = final_merged$target_start,
  pert_end = final_merged$target_end,
  gene_chr = final_merged$chr,
  gene_tss = final_merged$gene_tss,
  gene_strand = final_merged$strand.x,
  dist_to_tss = final_merged$distance,
  pert_level = "cre_perts",
  target_type = final_merged$target_type.x,
  cells = final_merged$n_nonzero_trt,
  avg_expr = NA,
  disp_outlier_deseq2 = FALSE
)

# Add power columns from final_merged
power_cols <- grep("PowerAtEffectSize", names(final_merged), value = TRUE)
final_output[power_cols] <- final_merged[power_cols]


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_output, snakemake@output$final_output)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)