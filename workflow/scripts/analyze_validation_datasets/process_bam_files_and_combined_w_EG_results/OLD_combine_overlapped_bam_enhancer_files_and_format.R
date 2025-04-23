# Script: combine_overlapped_bam_enhancer_files.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/combine_overlapped_bam_enhancer_files_and_format.rda"))
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

message("Loading overlapped files")
# Read in the validation overlaps
validation_files <- snakemake@input$validation_overlaps
validation_overlaps_list <- lapply(validation_files, read_tsv)

# Read in the training overlaps
training_files <- snakemake@input$training_overlaps
training_overlaps_list <- lapply(training_files, read_tsv)

# Read in the unfiltered k562 dc tap overlaps
unfiltered_k562_dc_tap_files <- snakemake@input$unfiltered_k562_dc_tap_overlaps
unfiltered_k562_dc_tap_overlaps_list <- lapply(unfiltered_k562_dc_tap_files, read_tsv)


### COMBINING DATASETS ========================================================

combine_overlaps <- function(overlaps_list) {
  # Initialize combined_df with the first dataframe
  combined_df <- overlaps_list[[1]]
  
  # Loop over the rest of the dataframes and merge
  if (length(overlaps_list) > 1) {
    for (i in 2:length(overlaps_list)) {
      df <- overlaps_list[[i]]
      # Merge the dataframes on common columns
      combined_df <- full_join(combined_df, df, by = intersect(names(combined_df), names(df)))
    }
  }
  return(combined_df)
}

# Combine training overlaps
message("Combining training data")
training_combined <- combine_overlaps(training_overlaps_list)

# Combine validation overlaps
message("Combining validation data")
validation_combined <- combine_overlaps(validation_overlaps_list)

# Combine unfiltered k562 dc tap overlaps
message("Combining unfiltered k562 dc tap data")
unfiltered_k562_dc_tap_combined <- combine_overlaps(unfiltered_k562_dc_tap_overlaps_list)

### FORMAT THE DATASETS ======================================================

# Let's create the `category` column for the combined_validation data
formatted_validation_combined <- validation_combined %>%
  mutate(category = case_when(
    Dataset == "Jurkat" ~ "Jurkat Nasser et al., 2021",
    Dataset == "GM12878" ~ "GM12878 Nasser et al., 2021",
    Dataset == "HCT116" ~ "HCT116 Guckelberger et al., inPrep",
    Dataset == "TAPSeq" ~ "K562 DC TAP Seq",
    Dataset == "WTC11" ~ "WTC11 DC TAP Seq",
    TRUE ~ Reference
  ))

# Let's create the `category` column for the combined_training data
formatted_training_combined <- training_combined %>%
  rename(Dataset = dataset) %>%
  mutate(category = case_when(
    Dataset == "TAPseq" ~ "Schraivogel et al., 2020",
    Dataset == "Gasperini2019" ~ "Gasperini et al., 2019",
    Dataset == "FlowFISH_K562" ~ "K562 FlowFISH"
  )) %>%
  # Training data has extra columns that we want to remove to match with the validation dataset
  select(colnames(formatted_validation_combined))

# Let's create the `category` column for the unfiltered_k562_dc_tap data
formatted_unfiltered_k562_dc_tap_combined <- unfiltered_k562_dc_tap_combined %>%
  mutate(category = "K562 DC TAP Seq")

### ADD CLASSIFICATION RESULTS AND FILTER ====================================

# Create compute_class function to compute classes of each formatted_combined_dataset file
compute_class <- function(formatted_combined_dataset) {
  
  # Add the threshold to compute the classification
  threshold <- snakemake@params$model_threshold
  
  formatted_combined_dataset <- formatted_combined_dataset %>%
    # We filter for only pairs that are evaluated in benchmarking - we're also only interested in the ENCODE model here
    filter(ValidConnection == TRUE, pred_id == "ENCODE_rE2G") %>%
    mutate(
      ClassificationResult = case_when(
        Regulated == TRUE & pred_value > threshold ~ "TP",
        Regulated == FALSE & pred_value > threshold ~ "FP",
        Regulated == FALSE & pred_value < threshold ~ "TN",
        Regulated == TRUE & pred_value < threshold ~ "FN"
      )
    )
  return(formatted_combined_dataset)
}

formatted_training_combined <- compute_class(formatted_training_combined)
formatted_validation_combined <- compute_class(formatted_validation_combined)
formatted_unfiltered_k562_dc_tap_combined <- compute_class(formatted_unfiltered_k562_dc_tap_combined)


### SAVE OUTPUT ===============================================================

# Save combined data
message("Saving combined datasets")
write_tsv(formatted_training_combined, snakemake@output$combined_training)
write_tsv(formatted_validation_combined, snakemake@output$combined_validation)
write_tsv(formatted_unfiltered_k562_dc_tap_combined, snakemake@output$combined_unfiltered_k562_dc_tap)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

