# Script: # Script: combine_encode_features_with_enhancer_classes.R

### SETUP =====================================================================

# Get the dataset wildcard
dataset_name <- snakemake@wildcards$dataset

# Saving image for debugging
save.image(paste0("RDA_objects/combine_encode_features_with_enhancer_classes_", dataset_name, ".rda"))
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
enhancer_classes <- read_tsv(snakemake@input$enhancer_classes)
features_file <- read_tsv(snakemake@input$features_file)


### MERGING DATASETS ==========================================================

# Get common columns except pred_value
merge_cols <- setdiff(intersect(colnames(enhancer_classes), colnames(features_file)), "pred_value")

# Merge with suffixes to indicate source
message("Merging enhancer classes and features file")
merged_data <- enhancer_classes %>% 
  left_join(features_file, 
            by = merge_cols,
            suffix = c(".enhancer_classes", ".features_file"))

### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(merged_data, snakemake@output$combined_file)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)