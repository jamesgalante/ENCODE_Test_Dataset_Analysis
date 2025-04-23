# Script: sample_training_dataset_for_class_balance.R

### SETUP =====================================================================

# Saving image for debugging
# Check if the directory exists, and create it only if it doesn't
if (!file.exists("RDA_objects/sample_training_dataset_for_class_balance")) { dir.create("RDA_objects/sample_training_dataset_for_class_balance", recursive = TRUE) }
save.image(paste0("RDA_objects/sample_training_dataset_for_class_balance/", snakemake@wildcards$subset, "_rep", snakemake@wildcards$replicate, ".rda"))
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
subset <- readRDS(snakemake@input$subset)
full_training <- read_tsv(snakemake@input$labelled_combined_training)


### RANDOMLY DOWNSAMPLE TRAINING ==============================================

# Count the number of positives and negatives in the subset
message("Calculating counts in subset")
positive_count <- sum(subset$Regulated == TRUE)
negative_count <- sum(subset$Regulated == FALSE)

# Filter full_training to get positives and negatives
full_training_positives <- full_training %>% filter(Regulated == TRUE)
full_training_negatives <- full_training %>% filter(Regulated == FALSE)

# Determine sample sizes (use all available if counts are insufficient)
positive_sample_size <- min(positive_count, nrow(full_training_positives))
negative_sample_size <- min(negative_count, nrow(full_training_negatives))

# Randomly sample from positives and negatives
message("Sampling from full training set to match these counts")
sampled_positives <- full_training_positives %>% sample_n(positive_sample_size, replace = FALSE)
sampled_negatives <- full_training_negatives %>% sample_n(negative_sample_size, replace = FALSE)

# Combine sampled positives and negatives
downsampled_training <- bind_rows(sampled_positives, sampled_negatives)

# Shuffle the combined dataset
downsampled_training <- downsampled_training %>% sample_frac(1)

### SAVE OUTPUT ===============================================================

# Save the downsampled training dataset
saveRDS(downsampled_training, file = snakemake@output$randomly_downsampled_training)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

