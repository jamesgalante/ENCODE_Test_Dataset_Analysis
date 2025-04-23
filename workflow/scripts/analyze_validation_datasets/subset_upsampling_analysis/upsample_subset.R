# Script: # Script: upsample_subset.R

### SETUP =====================================================================

# Saving image for debugging
# Check if the directory exists, and create it only if it doesn't
if (!file.exists("RDA_objects/upsample_subset")) { dir.create("RDA_objects/upsample_subset", recursive = TRUE) }
save.image(paste0("RDA_objects/upsample_subset/", snakemake@wildcards$subset, "_rep", snakemake@wildcards$replicate, ".rda"))
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


### UPSAMPLE SUBSET POSITIVES ================================================

# For Regulated == TRUE
regulated_true_subset <- subset %>% filter(Regulated == TRUE)
regulated_true_training <- full_training %>% filter(Regulated == TRUE)

# Initialize an empty list to store the balanced dataset
balanced_subset_true <- list()

for (cat in unique(regulated_true_training$Enhancer_Class)) {
  # Get counts for the current category in both subset and training
  subset_count <- regulated_true_subset %>% filter(Enhancer_Class == cat) %>% nrow()
  training_count <- regulated_true_training %>% filter(Enhancer_Class == cat) %>% nrow()
  
  # Always upsample subset to match training count if possible
  if (subset_count > 0) {
    sampled_data <- regulated_true_subset %>% 
      filter(Enhancer_Class == cat) %>%
      sample_n(training_count, replace = TRUE)
    
    # Add the sampled data to the list
    balanced_subset_true[[cat]] <- sampled_data
  }
}

# Combine all the sampled category data
balanced_subset_true <- bind_rows(balanced_subset_true)

# Check if we need to sample more to match total positives
remaining_true_samples <- nrow(regulated_true_training) - nrow(balanced_subset_true)

if (remaining_true_samples > 0) {
  # Sample from all positives to make up the difference
  extra_samples <- regulated_true_subset %>% 
    sample_n(remaining_true_samples, replace = TRUE)
  balanced_subset_true <- bind_rows(balanced_subset_true, extra_samples)
}


### UPSAMPLE SUBSET NEGATIVES ================================================

# For Regulated == FALSE
regulated_false_subset <- subset %>% filter(Regulated == FALSE)
regulated_false_training <- full_training %>% filter(Regulated == FALSE)

# Initialize an empty list to store the balanced dataset
balanced_subset_false <- list()

for (cat in unique(regulated_false_training$Enhancer_Class)) {
  # Get counts for the current category in both subset and training
  subset_count <- regulated_false_subset %>% filter(Enhancer_Class == cat) %>% nrow()
  training_count <- regulated_false_training %>% filter(Enhancer_Class == cat) %>% nrow()
  
  # Always upsample subset to match training count if possible
  if (subset_count > 0) {
    sampled_data <- regulated_false_subset %>% 
      filter(Enhancer_Class == cat) %>%
      sample_n(training_count, replace = TRUE)
    
    # Add the sampled data to the list
    balanced_subset_false[[cat]] <- sampled_data
  }
}

# Combine all the sampled category data
balanced_subset_false <- bind_rows(balanced_subset_false)

# Check if we need to sample more to match total negatives
remaining_false_samples <- nrow(regulated_false_training) - nrow(balanced_subset_false)

if (remaining_false_samples > 0) {
  # Sample from all negatives to make up the difference
  extra_samples <- regulated_false_subset %>% 
    sample_n(remaining_false_samples, replace = TRUE)
  balanced_subset_false <- bind_rows(balanced_subset_false, extra_samples)
}


### COMBINE BOTH DATAFRAMES ===================================================

# Final balanced subset data
balanced_subset <- bind_rows(balanced_subset_true, balanced_subset_false)

# Check the final result
balanced_subset %>% group_by(Regulated, Enhancer_Class) %>% summarize(n = n())

# Handle potential duplicates in names as in the original script
balanced_subset <- balanced_subset %>%
  group_by(name) %>%
  mutate(name = if(n() > 1) paste0(name, "_", row_number()) else name) %>%
  ungroup()


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(balanced_subset, snakemake@output$upsampled_subset)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)