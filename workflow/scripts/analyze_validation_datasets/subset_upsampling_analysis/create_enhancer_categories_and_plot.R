# Script: # Script: create_enhancer_categories_and_plot.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_enhancer_categories_and_plot.rda"))
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
combined_unfiltered_k562_dc_tap <- read_tsv(snakemake@input$combined_unfiltered_k562_dc_tap) %>% mutate(Dataset = "Unfilt. TAPSeq") %>% select(colnames(combined_validation))

### CREATE POSITIVES CATEGORIES ===============================================

# Define a function to apply the common operations with variable dhs/h3k27ac cutoffs
process_dataset <- function(data, category_h3k27ac, category_dhs) {
  data %>%
    mutate(DistToTSS = abs(((chromStart + chromEnd) / 2) - ((startTSS + endTSS) / 2))) %>%
    mutate(
      Enhancer_Class = case_when(
        `h3k27me3 overlap` == TRUE ~ "H3K27me3 Overlap",
        `ctcf overlap` == TRUE & `median_high_low_nonpromoter_h3k27ac` == "low" ~ "CTCF Overlap", # CTCF Overlap by median h3k27ac threshold always
        !!sym(category_h3k27ac) == "high" & !!sym(category_dhs) == "high" ~ "High H3K27ac High DNase",
        !!sym(category_h3k27ac) == "high" & !!sym(category_dhs) == "low" ~ "High H3K27ac Low DNase",
        !!sym(category_h3k27ac) == "low" & !!sym(category_dhs) == "low" ~ "Low H3K27ac Low DNase",
        !!sym(category_h3k27ac) == "low" & !!sym(category_dhs) == "high" ~ "Low H3K27ac High DNase",
        `ctcf overlap` == TRUE ~ "CTCF Overlap",
        is.na(`dhs_signal`) ~ "No peak overlap"
      )
    )
}

# Apply with 'high_low'
full_val_elbow <- process_dataset(combined_validation, "high_low_nonpromoter_h3k27ac", "high_low_nonpromoter_dhs")
full_train_elbow <- process_dataset(combined_training, "high_low_nonpromoter_h3k27ac", "high_low_nonpromoter_dhs")
unfilt_k562_dc_tap_elbow <- process_dataset(combined_unfiltered_k562_dc_tap, "high_low_nonpromoter_h3k27ac", "high_low_nonpromoter_dhs")

# Apply with 'cdf75_high_low'
full_val_cdf75 <- process_dataset(combined_validation, "cdf75_high_low_nonpromoter_h3k27ac", "cdf75_high_low_nonpromoter_dhs")
full_train_cdf75 <- process_dataset(combined_training, "cdf75_high_low_nonpromoter_h3k27ac", "cdf75_high_low_nonpromoter_dhs")
unfilt_k562_dc_tap_cdf75 <- process_dataset(combined_unfiltered_k562_dc_tap, "cdf75_high_low_nonpromoter_h3k27ac", "cdf75_high_low_nonpromoter_dhs")

# Apply with 'median_high_low'
full_val_median <- process_dataset(combined_validation, "median_high_low_nonpromoter_h3k27ac", "median_high_low_nonpromoter_dhs")
full_train_median <- process_dataset(combined_training, "median_high_low_nonpromoter_h3k27ac", "median_high_low_nonpromoter_dhs")
unfilt_k562_dc_tap_median <- process_dataset(combined_unfiltered_k562_dc_tap, "median_high_low_nonpromoter_h3k27ac", "median_high_low_nonpromoter_dhs")


### COMMON FUNCTION TO CREATE PLOTS ===========================================

# Function to create a plot based on combined data
create_plot <- function(combined_data, title_suffix) {
  
  # Calculate total counts for each dataset to include in the x-axis labels
  dataset_counts <- combined_data %>%
    group_by(Source, Dataset) %>%
    summarise(TotalCount = n())
  
  # Merge the counts back into the data for plotting
  prop_data <- combined_data %>%
    group_by(Source, Dataset, Enhancer_Class) %>%
    summarise(Count = n()) %>%
    left_join(dataset_counts, by = c("Source", "Dataset")) %>%
    mutate(
      Proportion = Count / sum(Count),  # Proportion of each enhancer class
      y_position = cumsum(Proportion) - Proportion / 2,  # Midpoint of each segment for label positioning
      DatasetLabel = paste0(Dataset)  # Dataset label
    )
  
  # Calculate combined proportions across all datasets within each source
  combined_source_data <- combined_data %>%
    group_by(Source, Enhancer_Class) %>%
    summarise(Count = n()) %>%
    mutate(
      Dataset = "All Datasets",
      TotalCount = sum(Count),
      Proportion = Count / sum(Count),
      DatasetLabel = paste0("All Datasets")
    )
  
  # Combine the original data with the combined source data
  prop_data <- bind_rows(prop_data, combined_source_data)
  
  # Create a stacked bar plot showing proportions with labels at the base of each class
  p <- ggplot(prop_data, aes(x = DatasetLabel, y = Proportion, fill = Enhancer_Class)) +
    geom_bar(stat = "identity") +  # Use identity to stack bars by proportion
    geom_text(aes(label = Count, y = Proportion), position = position_stack(vjust = 0.5), color = "black", size = 2.5) +  # Add count labels at the center of each segment
    geom_text(aes(x = DatasetLabel, y = -0.05, label = paste0("(", TotalCount, ")")), size = 3, vjust = 0) +  # Add total count label at the top of each bar
    labs(
      title = paste0("Proportion of Enhancer-Gene Classes Across Datasets ", title_suffix),
      x = "Dataset",
      y = "Proportion of EG Pairs",
      fill = "Enhancer Class"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Slant the x-axis labels
      strip.background = element_rect(fill = "lightgray"),
      panel.grid = element_blank()
    ) +
    facet_grid(. ~ Source, scales = "free", space = "free") +  # Add facet grid with both Source and combined datasets
    scale_fill_brewer(palette = "Set3")  # Optional: Change color palette
  
  # Create the heatmap plot
  q <- ggplot(prop_data, aes(x = Dataset, y = Enhancer_Class, fill = Proportion)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "lightgray", high = "steelblue") +
    labs(
      title = paste0("Enhancer Classes Heatmap ", title_suffix),
      x = "Dataset",
      y = "Enhancer Class",
      fill = "Proportion"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    facet_grid(. ~ Source, scales = "free_x", space = "free_x")
  
  return(list(p = p, q = q))
}

### CREATE PLOTS ==============================================================

# Combine the datasets for each version
combined_data_elbow <- bind_rows(
  full_val_elbow %>% mutate(Source = "Validation"),
  full_train_elbow %>% mutate(Source = "Training"),
  unfilt_k562_dc_tap_elbow %>% mutate(Source = "DC TAP")
)

combined_data_cdf75 <- bind_rows(
  full_val_cdf75 %>% mutate(Source = "Validation"),
  full_train_cdf75 %>% mutate(Source = "Training"),
  unfilt_k562_dc_tap_cdf75 %>% mutate(Source = "DC TAP")
)

combined_data_median <- bind_rows(
  full_val_median %>% mutate(Source = "Validation"),
  full_train_median %>% mutate(Source = "Training"),
  unfilt_k562_dc_tap_median %>% mutate(Source = "DC TAP")
)

# Generate plots for high_low
plots_elbow_positives <- create_plot(combined_data_elbow %>% filter(Regulated == TRUE), "(Elbow High Low - CRISPRi Positives)")
plots_elbow_negatives <- create_plot(combined_data_elbow %>% filter(Regulated == FALSE), "(Elbow High Low - CRISPRi Negatives)")
plots_elbow_all <- create_plot(combined_data_elbow, "(Elbow High Low - All Pairs)")

# Generate plots for cdf75_high_low
plots_cdf75_positives <- create_plot(combined_data_cdf75 %>% filter(Regulated == TRUE), "(CDF75 High Low - CRISPRi Positives)")
plots_cdf75_negatives <- create_plot(combined_data_cdf75 %>% filter(Regulated == FALSE), "(CDF75 High Low - CRISPRi Negatives)")
plots_cdf75_all <- create_plot(combined_data_cdf75, "(CDF75 High Low - All Pairs)")

# Generate plots for median_high_low
plots_median_positives <- create_plot(combined_data_median %>% filter(Regulated == TRUE), "(Median High Low - CRISPRi Positives)")
plots_median_negatives <- create_plot(combined_data_median %>% filter(Regulated == FALSE), "(Median High Low - CRISPRi Negatives)")
plots_median_all <- create_plot(combined_data_median, "(Median High Low - All Pairs)")


### SAVE OUTPUT ===============================================================

# Save the bar plots
ggsave(plot = plots_elbow_positives$p, filename = snakemake@output$enhancer_class_crispri_positives_plot,
       device = "pdf", height = 6, width = 11)
ggsave(plot = plots_elbow_negatives$p, filename = snakemake@output$enhancer_class_crispri_negatives_plot,
       device = "pdf", height = 6, width = 11)
ggsave(plot = plots_elbow_all$p, filename = snakemake@output$enhancer_class_all_pairs_plot, 
       device = "pdf", height = 6, width = 11)

# Save the heatmaps
ggsave(plot = plots_elbow_positives$q, filename = snakemake@output$enhancer_class_crispri_positives_heatmap_plot,
       device = "pdf", height = 6, width = 11)
ggsave(plot = plots_elbow_negatives$q, filename = snakemake@output$enhancer_class_crispri_negatives_heatmap_plot,
       device = "pdf", height = 6, width = 11)
ggsave(plot = plots_elbow_all$q, filename = snakemake@output$enhancer_class_all_pairs_heatmap_plot, 
       device = "pdf", height = 6, width = 11)

# Save each combined dataset with enhancer class labels
write_tsv(unfilt_k562_dc_tap_elbow, snakemake@output$labelled_combined_unfiltered_k562_dc_tap)
write_tsv(full_val_elbow, snakemake@output$labelled_combined_validation)
write_tsv(full_train_elbow, snakemake@output$labelled_combined_training)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

