# Script: enhancer_category_proportions.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/enhancer_category_proportions.rda"))
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
labelled_combined_validation <- read_tsv(snakemake@input$labelled_combined_validation) %>% mutate(Source = "Validation")
labelled_combined_training <- read_tsv(snakemake@input$labelled_combined_training) %>% mutate(Source = "Training")


### COMBINE THE DATASETS =====================================================

# Create combined version with category and Source (Basically so we can have one column with the full val or training dataset in the plot)
combined_labelled_combined_validation <- labelled_combined_validation %>% mutate(category = "All Validation Datasets", Source = "Combined")
combined_labelled_combined_training <- labelled_combined_training %>% mutate(category = "All Training Datasets", Source = "Combined")

# Combine all of these
all_datasets <- rbind(labelled_combined_training, labelled_combined_validation,
                      combined_labelled_combined_validation, combined_labelled_combined_training)


### CONSOLIDATE ENHANCER CATEGORIES ===========================================

# Add new column for consolidated enhancer categories
consolidate_enhancer_categories <- function(data) {
  data %>%
    mutate(plotting_enhancer_categories = case_when(
      Enhancer_Class %in% c("High H3K27ac High DNase", "High H3K27ac Low DNase") ~ "Strongly acetylated elements (High H3K27ac)",
      Enhancer_Class == "CTCF Overlap" ~ "CTCF elements (CTCF, low H3K27ac)",
      Enhancer_Class == "H3K27me3 Overlap" ~ "Polycomb elements (H3K27me3)",
      TRUE ~ "Weakly acetylated elements (Low H3K27ac, no CTCF)"
    ))
}
all_datasets <- consolidate_enhancer_categories(all_datasets)


### CREATE THE PROPORTIONS PLOT ==============================================

# Function to create plot with new enhancer categories, allowing optional labeling
create_plot <- function(combined_data, title_suffix, showCategoryNumbers = TRUE, showDatasetTotals = TRUE) {
  
  # First calculate total counts per dataset
  dataset_totals <- combined_data %>%
    group_by(Source, category) %>%
    summarise(TotalCount = n(), .groups = "drop")
  
  # Calculate proportions and prepare data for plotting
  prop_data <- combined_data %>%
    group_by(Source, category, plotting_enhancer_categories) %>%
    summarise(Count = n(), .groups = "drop") %>%
    left_join(dataset_totals, by = c("category", "Source")) %>%
    mutate(
      Proportion = Count / TotalCount,
      y_position = cumsum(Proportion) - Proportion / 2
    )
  
  text_angle <- ifelse(!showCategoryNumbers, 90, 45)
  
  # Create the base stacked bar plot
  p <- ggplot(prop_data, aes(x = category, y = Proportion, fill = plotting_enhancer_categories)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Proportion of Enhancer-Gene Classes Across Datasets",
      subtitle = title_suffix,
      x = "Dataset",
      y = "Proportion of EG Pairs",
      fill = "Element Category"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = text_angle, hjust = 1),
      strip.background = element_rect(fill = "lightgray"),
      panel.grid = element_blank()
    ) +
    facet_grid(. ~ Source, scales = "free", space = "free") +
    scale_fill_manual(values = c(
      "Strongly acetylated elements (High H3K27ac)" = "#8DD3C7",
      "CTCF elements (CTCF, low H3K27ac)" = "#FFFFB3",
      "Polycomb elements (H3K27me3)" = "#BEBADA",
      "Weakly acetylated elements (Low H3K27ac, no CTCF)" = "#FB8072"
    ))
  
  # Optionally add number of pairs in each enhancer class on the stacked bars
  if (showCategoryNumbers) {
    p <- p + 
      geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black", size = 2.5)
  }
  
  # Optionally add the total number of pairs under each bar
  if (showDatasetTotals) {
    p <- p + 
      geom_text(aes(x = category, y = -0.05, label = paste0("(", TotalCount, ")")), size = 3, vjust = 0)
  }
  
  return(p)
}

# Create the full and mini plots
all_proportions_plot_full <- create_plot(combined_data = all_datasets, title_suffix = "All Pairs", showCategoryNumbers = TRUE, showDatasetTotals = TRUE)
all_proportions_plot_mini <- create_plot(combined_data = all_datasets, title_suffix = "All Pairs", showCategoryNumbers = FALSE, showDatasetTotals = FALSE)
positives_proportions_plot_full <- create_plot(combined_data = all_datasets %>% filter(Regulated == TRUE), title_suffix = "CRISPRi Positives", showCategoryNumbers = TRUE, showDatasetTotals = TRUE)
positives_proportions_plot_mini <- create_plot(combined_data = all_datasets %>% filter(Regulated == TRUE), title_suffix = "CRISPRi Positives", showCategoryNumbers = FALSE, showDatasetTotals = FALSE)


### CREATE SELECTED PLOTS =====================================================

# (All training datasets, Both Combined datasets, Both DC TAP datasets, K562 Only, Other Cell Types Only)
# Need to create the K562 Only and other cell types only in addition to subsetting for the others
K562_Only <- labelled_combined_validation %>% filter(ExperimentCellType == "K562") %>% mutate(Source = "Validation", category = "K562 Only")
other_cell_types <- labelled_combined_validation %>% filter(ExperimentCellType != "K562") %>% mutate(Source = "Validation", category = "Other Cell Types Only")

# Grab the dc tap datasets
dc_tap_datasets <- labelled_combined_validation %>% filter(category %in% c("WTC11 DC TAP Seq", "K562 DC TAP Seq")) %>% mutate(Source = "Validation")

# Combine the selected categories together
selected_datasets <- rbind(labelled_combined_training, combined_labelled_combined_training, combined_labelled_combined_validation, K562_Only, other_cell_types, dc_tap_datasets)

# Create the plots
all_proportions_plot_select_datasets <- create_plot(combined_data = consolidate_enhancer_categories(selected_datasets), title_suffix = "All Pairs", showCategoryNumbers = TRUE, showDatasetTotals = TRUE)
positives_proportions_plot_select_datasets <- create_plot(combined_data = consolidate_enhancer_categories(selected_datasets) %>% filter(Regulated == TRUE), title_suffix = "CRISPRi Positives", showCategoryNumbers = TRUE, showDatasetTotals = TRUE)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = all_proportions_plot_full, file = snakemake@output$all_proportions_plot_full, device = "pdf", height = 6, width = 10)
ggsave(plot = all_proportions_plot_mini, file = snakemake@output$all_proportions_plot_mini, device = "pdf", height = 4.5, width = 7)
ggsave(plot = positives_proportions_plot_full, file = snakemake@output$positives_proportions_plot_full, device = "pdf", height = 6, width = 10)
ggsave(plot = positives_proportions_plot_mini, file = snakemake@output$positives_proportions_plot_mini, device = "pdf", height = 4.5, width = 7)
ggsave(plot = all_proportions_plot_select_datasets, file = snakemake@output$all_proportions_plot_select_datasets, device = "pdf", height = 6, width = 10)
ggsave(plot = positives_proportions_plot_select_datasets, file = snakemake@output$positives_proportions_plot_select_datasets, device = "pdf", height = 6, width = 10)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)