# Script: downsample_training_and_subset.R

### SETUP =====================================================================

# Saving image for debugging
# Check if the directory exists, and create it only if it doesn't
if (!file.exists("RDA_objects/downsample_training_and_subset")) { dir.create("RDA_objects/downsample_training_and_subset", recursive = TRUE) }
save.image(paste0("RDA_objects/downsample_training_and_subset/", snakemake@wildcards$subset, "_rep", snakemake@wildcards$replicate, ".rda"))
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


### DOWN SAMPLE TRAINING AND SUBSET ===========================================

# Function to downsample two datasets to the minimum count per category
downsample_to_min <- function(training_data, subset_data) {
  # Get unique categories
  categories <- unique(c(training_data$Enhancer_Class, subset_data$Enhancer_Class))
  
  # Initialize lists to store downsampled data
  downsampled_training_list <- list()
  downsampled_subset_list <- list()
  
  for (cat in categories) {
    # Filter data for the current category
    training_cat <- training_data %>% filter(Enhancer_Class == cat)
    subset_cat <- subset_data %>% filter(Enhancer_Class == cat)
    
    # Get counts
    training_count <- nrow(training_cat)
    subset_count <- nrow(subset_cat)
    
    # Determine the minimum count
    min_count <- min(training_count, subset_count)
    
    # Downsample both datasets to the minimum count if counts are greater than zero
    if (min_count > 0) {
      sampled_training <- training_cat %>% sample_n(min_count)
      sampled_subset <- subset_cat %>% sample_n(min_count)
      
      # Add to lists
      downsampled_training_list[[cat]] <- sampled_training
      downsampled_subset_list[[cat]] <- sampled_subset
    }
  }
  
  # Combine all categories
  downsampled_training <- bind_rows(downsampled_training_list)
  downsampled_subset <- bind_rows(downsampled_subset_list)
  
  return(list(training = downsampled_training, subset = downsampled_subset))
}

message("Processing Regulated == TRUE")

# Process positives (Regulated == TRUE)
regulated_true_training <- full_training %>% filter(Regulated == TRUE)
regulated_true_subset <- subset %>% filter(Regulated == TRUE)

downsampled_true <- downsample_to_min(regulated_true_training, regulated_true_subset)

message("Processing Regulated == FALSE")

# Process negatives (Regulated == FALSE)
regulated_false_training <- full_training %>% filter(Regulated == FALSE)
regulated_false_subset <- subset %>% filter(Regulated == FALSE)

downsampled_false <- downsample_to_min(regulated_false_training, regulated_false_subset)


### COMBINE BOTH DATAFRAMES ===================================================

# Combine positives and negatives for training and subset
downsampled_training <- bind_rows(downsampled_true$training, downsampled_false$training)
downsampled_subset <- bind_rows(downsampled_true$subset, downsampled_false$subset)


### SAVE OUTPUT ===============================================================

# Save the downsampled subsets
saveRDS(downsampled_training, file = snakemake@output$downsampled_training)
saveRDS(downsampled_subset, file = snakemake@output$downsampled_original)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

