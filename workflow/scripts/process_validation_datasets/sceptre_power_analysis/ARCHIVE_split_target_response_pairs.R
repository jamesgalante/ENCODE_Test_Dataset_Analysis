# Script: split_target_response_pairs.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/split_target_response_pairs")) { dir.create("RDA_objects/split_target_response_pairs", recursive = TRUE) }
save.image(paste0("RDA_objects/split_target_response_pairs/", snakemake@wildcards$sample, ".rda"))
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

message("Loading inputs")
gene_gRNA_group_pairs <- readRDS(snakemake@input$gene_gRNA_group_pairs)
batches <- snakemake@params$batches  # Get number of batches from snakemake params
output_files <- snakemake@output$splits


### SPLIT PAIRS FILE =======================================================

# Calculate the target number of unique grna_group for each split
message("Precomputations for splitting pairs")
total_unique_groups <- n_distinct(gene_gRNA_group_pairs$grna_target)
target_per_split <- ceiling(total_unique_groups / batches)

# Make sure that splitting the unique groups will work given the number of batches
if (total_unique_groups < batches) {
  stop(paste0("The number of batches (", batches, ") exceeds the number of unique grna_targets."))
}

# Initialize splits with empty data frames
splits <- vector("list", batches)
names(splits) <- paste0("split", seq_len(batches))
for (i in seq_len(batches)) {
  splits[[i]] <- data.frame(grna_target = character(0), response_id = character(0))
}

# Distribute grna_target to splits trying to even out the number of unique values
message("Distributing gene-gRNA group pairs")
unique_targets <- unique(gene_gRNA_group_pairs$grna_target)
for (i in seq_along(unique_targets)) {
  split_counts <- sapply(splits, function(x) n_distinct(x$grna_target))
  split_with_least <- which.min(split_counts)
  
  # Get the rows for the current grna_target
  current_rows <- gene_gRNA_group_pairs %>% 
    filter(grna_target == unique_targets[i])
  
  # Add the current rows to the appropriate split
  splits[[split_with_least]] <- bind_rows(splits[[split_with_least]], current_rows)
}

### SAVE OUTPUT ===============================================================

# Write each split to the corresponding output file
message("Saving splits to files")
for (i in seq_along(splits)) {
  write_tsv(splits[[i]], file=output_files[[i]], col_names = FALSE)
}


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)