# Script: plot_auprc_results.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/plot_auprc_results.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
})

message("Loading input files")
datasets_list <- list(
  "Original" = lapply(snakemake@input$original_subsets, readRDS),
  "Balanced Training" = lapply(snakemake@input$balanced_training_subsets, readRDS),
  "Downsampled Training" = lapply(snakemake@input$downsampled_training_subsets, readRDS),
  "Downsampled Original" = lapply(snakemake@input$downsampled_original_subsets, readRDS),
  "Randomly Downsampled Training" = lapply(snakemake@input$randomly_downsampled_training_subsets, readRDS),
  "Upsampled Original" = lapply(snakemake@input$upsampled_original_subsets, readRDS)
  # Add more datasets as needed
)

# Extract the subset names from file paths
message("Grabbing the subset names")
subset_names <- snakemake@params$subset_names
n_replicates <- snakemake@params$n_replicates


### EXTRACT AND AVERAGE VALUES =================================================

# Helper function to process performance summaries
process_perf_summary <- function(data_obj, subset_name, dataset_type, replicate_id = NULL) {
  perf_summary <- data_obj$performance$perf_summary
  
  # Select only the last row or the single row
  perf_summary <- perf_summary[nrow(perf_summary), ]
  
  # Add Subset, DatasetType, and replicate info (if applicable)
  perf_summary$Subset <- subset_name
  perf_summary$DatasetType <- dataset_type
  if (!is.null(replicate_id)) {
    perf_summary$Replicate <- replicate_id
  }
  
  # Ensure required columns exist, adding NA if missing
  required_columns <- c("AUPRC", "AUPRC_lowerCi", "AUPRC_upperCi")
  for (col in required_columns) {
    if (!col %in% colnames(perf_summary)) {
      perf_summary[[col]] <- NA  # Assign NA if column is missing
    }
  }
  
  return(perf_summary)
}

# Initialize lists to store data frames
perf_data_list <- list()
replicate_points_list <- list()

# Iterate over dataset types
for (dataset_type in names(datasets_list)) {
  # Determine whether the dataset type has replicates
  is_replicated <- dataset_type != "Original"
  
  # Process each subset within the dataset type
  for (i in seq_along(subset_names)) {
    subset_name <- subset_names[i]
    
    if (is_replicated) {
      # Process each replicate for replicated dataset types
      replicates <- datasets_list[[dataset_type]][((i - 1) * n_replicates + 1):(i * n_replicates)]
      
      # Store replicate points
      replicate_summaries <- lapply(seq_along(replicates), function(rep_id) {
        process_perf_summary(replicates[[rep_id]], subset_name, dataset_type, replicate_id = rep_id)
      })
      
      # Add replicate data for plotting points
      replicate_points_list <- append(replicate_points_list, replicate_summaries)
      
      # Average replicates for bar plot
      replicate_summaries_df <- bind_rows(replicate_summaries)
      perf_summary <- replicate_summaries_df %>%
        group_by(Subset, DatasetType) %>%
        summarize(
          AUPRC = mean(AUPRC, na.rm = TRUE),
          AUPRC_lowerCi = mean(AUPRC_lowerCi, na.rm = TRUE),
          AUPRC_upperCi = mean(AUPRC_upperCi, na.rm = TRUE)
        ) %>%
        ungroup()
    } else {
      # For non-replicated dataset types (Original)
      data_obj <- datasets_list[[dataset_type]][[i]]
      perf_summary <- process_perf_summary(data_obj, subset_name, dataset_type)
    }
    
    # Add to perf_data_list for final bar plot data
    if (!is.null(perf_summary)) {
      perf_data_list <- append(perf_data_list, list(perf_summary))
    }
  }
}

# Combine all data frames into one
perf_data <- bind_rows(perf_data_list)
replicate_points <- bind_rows(replicate_points_list)

### PLOT VALUES ===============================================================

# Ensure numeric columns are properly formatted
perf_data <- perf_data %>% mutate(
  AUPRC = as.numeric(AUPRC),
  AUPRC_lowerCi = as.numeric(AUPRC_lowerCi),
  AUPRC_upperCi = as.numeric(AUPRC_upperCi)
)

# Create the plot 
p <- ggplot(perf_data, aes(x = Subset, y = AUPRC, fill = DatasetType)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", size = 0.3) +  # Adding borders to bars
  geom_errorbar(aes(ymin = AUPRC_lowerCi, ymax = AUPRC_upperCi),
                width = 0.2,
                position = position_dodge(width = 0.8)) +
  geom_point(data = replicate_points,  # Plot points for replicated datasets
             aes(x = Subset, y = AUPRC),
             color = "black",
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.7) +  # Adjust size as desired
  scale_fill_brewer(palette = "Set1") +  # Apply Set1 color palette from RColorBrewer
  theme_minimal(base_size = 14) +  # Using a minimal theme
  ylim(0, 1) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 12, color = "gray20"),  # Rotate x-axis labels for readability
    axis.text.y = element_text(size = 12, color = "gray20"),  # Adjust y-axis text size
    axis.title.x = element_text(size = 14, face = "bold", color = "gray20"),  # x-axis title
    axis.title.y = element_text(size = 14, face = "bold", color = "gray20"),  # y-axis title
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centered title
    legend.position = "top",  # Move legend to the top
    legend.title = element_blank(),  # Remove legend title
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Add gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  ) +
  labs(
    title = "AUPRC Comparison by Subset",
    x = "Subset",
    y = "AUPRC",
    fill = "Dataset"  # Legend label
  )

### SAVE OUTPUT ===============================================================

ggsave(plot = p, filename = snakemake@output$auprc_upsampling_summary, width = 10, height = 7)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

