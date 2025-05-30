# Script: combining_dc_tap_and_validation_dataset.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/combining_dc_tap_and_validation_dataset.rda"))
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
formatted_resized_and_merged_dc_tap_output <- read_tsv(snakemake@input$formatted_resized_and_merged_dc_tap_output)
validation_dataset <- read_tsv(snakemake@input$validation_dataset)


### COMBINE THE DATASETS ======================================================

# Combine the datasets
full_validation_dataset_w_dc_tap_seq <- rbind(validation_dataset, formatted_resized_and_merged_dc_tap_output)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(full_validation_dataset_w_dc_tap_seq, snakemake@output$full_validation_dataset_w_dc_tap_seq)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)