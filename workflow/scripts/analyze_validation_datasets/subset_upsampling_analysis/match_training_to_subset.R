# Script: # Script: match_training_to_subset.R

### SETUP =====================================================================

# Saving image for debugging
# Check if the directory exists, and create it only if it doesn't
if (!file.exists("RDA_objects/match_training_to_subset")) { dir.create("RDA_objects/match_training_to_subset", recursive = TRUE) }
save.image(paste0("RDA_objects/match_training_to_subset/", snakemake@wildcards$subset, "_rep", snakemake@wildcards$replicate, ".rda"))
message("Saved Image")
#stop("Manually Stopped Program after Saving Image")

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


### FILTER TRAINING POSTIIVES =================================================

# For Regulated == TRUE
regulated_true_subset <- subset %>% filter(Regulated == TRUE)
regulated_true_training <- full_training %>% filter(Regulated == TRUE)

# Initialize an empty list to store the balanced dataset
balanced_training_true <- list()

for (cat in unique(regulated_true_training$Enhancer_Class)) {
  
  # Get counts for the current category in both subset and training
  subset_count <- regulated_true_subset %>% filter(Enhancer_Class == cat) %>% nrow()
  training_count <- regulated_true_training %>% filter(Enhancer_Class == cat) %>% nrow()
  
  if (training_count > subset_count) {
    # Downsample training
    sampled_data <- regulated_true_training %>% filter(Enhancer_Class == cat) %>%
      sample_n(subset_count)
  } else if (training_count < subset_count) {
    # Upsample training with replacement
    sampled_data <- regulated_true_training %>% filter(Enhancer_Class == cat) %>%
      sample_n(subset_count, replace = TRUE)
  } else {
    # If counts match, just take the training data as is
    sampled_data <- regulated_true_training %>% filter(Enhancer_Class == cat)
  }
  
  # Add the sampled data to the list
  balanced_training_true[[cat]] <- sampled_data
}

# Combine all the sampled data
balanced_training_true <- bind_rows(balanced_training_true)

# After processing all categories, check if there is still a difference in total counts
remaining_true_positives <- nrow(regulated_true_subset) - nrow(balanced_training_true)

# If there is a difference, sample from the entire training pool to make the counts equal
if (remaining_true_positives > 0) {
  extra_samples <- regulated_true_training %>% sample_n(remaining_true_positives, replace = TRUE)
  balanced_training_true <- bind_rows(balanced_training_true, extra_samples)
}


### FILTER TRAINING NEGATIVES =================================================

# For Regulated == FALSE
regulated_false_subset <- subset %>% filter(Regulated == FALSE)
regulated_false_training <- full_training %>% filter(Regulated == FALSE)

# Initialize an empty list to store the balanced dataset
balanced_training_false <- list()

for (cat in unique(regulated_false_training$Enhancer_Class)) {
  
  # Get counts for the current category in both subset and training
  subset_count <- regulated_false_subset %>% filter(Enhancer_Class == cat) %>% nrow()
  training_count <- regulated_false_training %>% filter(Enhancer_Class == cat) %>% nrow()
  
  if (training_count > subset_count) {
    # Downsample training
    sampled_data <- regulated_false_training %>% filter(Enhancer_Class == cat) %>%
      sample_n(subset_count)
  } else if (training_count < subset_count) {
    # Upsample training with replacement
    sampled_data <- regulated_false_training %>% filter(Enhancer_Class == cat) %>%
      sample_n(subset_count, replace = TRUE)
  } else {
    # If counts match, just take the training data as is
    sampled_data <- regulated_false_training %>% filter(Enhancer_Class == cat)
  }
  
  # Add the sampled data to the list
  balanced_training_false[[cat]] <- sampled_data
}

# Combine all the sampled data
balanced_training_false <- bind_rows(balanced_training_false)

# After processing all categories, check if there is still a difference in total counts
remaining_false_negatives <- nrow(regulated_false_subset) - nrow(balanced_training_false)

# If there is a difference, sample from the entire training pool to make the counts equal
if (remaining_false_negatives > 0) {
  extra_samples <- regulated_false_training %>% sample_n(remaining_false_negatives, replace = TRUE)
  balanced_training_false <- bind_rows(balanced_training_false, extra_samples)
}


### COMBINE BOTH DATAFRAMES ===================================================

# Final balanced training data
balanced_training <- bind_rows(balanced_training_true, balanced_training_false)

# Check the final result
balanced_training %>% group_by(Regulated, Enhancer_Class) %>% summarize(n = n())

# In `convertMergedForBootstrap`, there is a pivot_wider step that requires unique names in the dataset
# We can give any duplicates (upsampled pairs) a distinct name to abide by this
balanced_training <- balanced_training %>%
  group_by(name) %>%
  mutate(name = if(n() > 1) paste0(name, "_", row_number()) else name) %>%
  ungroup()


### SAVE OUTPUT ===============================================================

# Save the filtered training dataset
saveRDS(balanced_training, snakemake@output$balanced_training)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

