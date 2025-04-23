# Script: # Script: create_enhancer_subsets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_enhancer_subsets.rda"))
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
combined_validation <- read_tsv(snakemake@input$combined_validation)
combined_training <- read_tsv(snakemake@input$combined_training)


### CREATING SUBSETS ==========================================================

# Create each individual subset
message("Subsetting datasets")
full_validation_dataset <- combined_validation
k562_validation_dataset <- combined_validation %>% filter(ExperimentCellType == "K562")
k562_validation_wo_dc_tap <- combined_validation %>% filter(ExperimentCellType == "K562", category != "K562 DC TAP Seq")
training_dataset <- combined_training
other_cell_types <- combined_validation %>% filter(ExperimentCellType != "K562")
validation_without_dc_tap <- combined_validation %>% filter(!category %in% c("K562 DC TAP Seq", "WTC11 DC TAP Seq"))
k562_dc_tap <- combined_validation %>% filter(category == "K562 DC TAP Seq")
wtc11_dc_tap <- combined_validation %>% filter(category == "WTC11 DC TAP Seq")
wtc11_k562_validation <- combined_validation %>% filter(ExperimentCellType %in% c("K562", "WTC11"))
both_dc_tap_only <- combined_validation %>% filter(category %in% c("K562 DC TAP Seq", "WTC11 DC TAP Seq"))
other_cell_types_wo_dc_tap <- combined_validation %>% filter(!ExperimentCellType %in% c("WTC11", "K562"))
validation_wo_wtc11 <- combined_validation %>% filter(ExperimentCellType != "WTC11")
validation_wo_k562_dc_tap <- combined_validation %>% filter(category != "K562 DC TAP Seq")

# Create a subset list for saving the files
subsets <- list(
  full_validation_dataset = full_validation_dataset,
  k562_validation_dataset = k562_validation_dataset,
  k562_validation_wo_dc_tap = k562_validation_wo_dc_tap,
  training_dataset = training_dataset,
  other_cell_types = other_cell_types,
  validation_without_dc_tap = validation_without_dc_tap,
  k562_dc_tap = k562_dc_tap,
  wtc11_dc_tap = wtc11_dc_tap,
  wtc11_k562_validation = wtc11_k562_validation,
  both_dc_tap_only = both_dc_tap_only,
  other_cell_types_wo_dc_tap = other_cell_types_wo_dc_tap,
  validation_wo_wtc11 = validation_wo_wtc11,
  validation_wo_k562_dc_tap = validation_wo_k562_dc_tap
)


### SAVE OUTPUT ===============================================================

# Save each subset to its respective output file
message("Saving dataset subsets")
for (i in seq_along(subsets)) {
  saveRDS(subsets[[i]], file = snakemake@output[[i]])
  message(paste("Saved subset:", names(subsets)[i], "to file:", snakemake@output[[i]]))
}


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)