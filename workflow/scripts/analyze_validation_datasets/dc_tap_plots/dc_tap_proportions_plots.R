# Script: dc_tap_proportions_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/dc_tap_proportions_plots.rda"))
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
labelled_combined_validation <- read_tsv(snakemake@input$labelled_combined_validation)
labelled_combined_training <- read_tsv(snakemake@input$labelled_combined_training)
labelled_combined_unfiltered_k562_dc_tap <- read_tsv(snakemake@input$labelled_combined_unfiltered_k562_dc_tap) %>% mutate(Dataset = "Unfilt. TAPSeq") %>% select(colnames(labelled_combined_training))

# Filter for all "control" genes that weren't randomly chosen
biased_genes_to_remove <- read_tsv(snakemake@input$biased_genes_to_remove) %>% filter(type == "control") %>% pull(gene) %>% unique()


### FILTER FOR DATASETS =======================================================

# Select Gasperini and Fulco FlowFISH from the training data
training_datasets <- labelled_combined_training %>% filter(Dataset %in% c("FlowFISH_K562", "Gasperini2019"))

# Select the WTC11 from the combined_validation dataset
wtc11_dc_tap <- labelled_combined_validation %>% filter(category == "WTC11 DC TAP Seq")
# Remove genes that weren't randomly chosen
wtc11_dc_tap <- wtc11_dc_tap %>% filter(!measuredGeneSymbol %in% biased_genes_to_remove)

# Remove genes that weren't randomly chosen
labelled_combined_unfiltered_k562_dc_tap <- labelled_combined_unfiltered_k562_dc_tap %>% filter(!measuredGeneSymbol %in% biased_genes_to_remove)
# Combine all of these with the unfiltered_k562_dc_tap
dc_tap_and_training <- rbind(training_datasets, wtc11_dc_tap, labelled_combined_unfiltered_k562_dc_tap)


### CONSOLIDATE ENHANCER CATEGORIES ===========================================

# Add new column for consolidated enhancer categories
dc_tap_and_training <- dc_tap_and_training %>%
  mutate(plotting_enhancer_categories = case_when(
    Enhancer_Class %in% c("High H3K27ac High DNase", "High H3K27ac Low DNase") ~ 
      "Strongly acetylated elements (High H3K27ac)",
    Enhancer_Class == "CTCF Overlap" ~ 
      "CTCF elements (CTCF, low H3K27ac)",
    Enhancer_Class == "H3K27me3 Overlap" ~ 
      "Polycomb elements (H3K27me3)",
    TRUE ~ "Weakly acetylated elements (Low H3K27ac, no CTCF)"
  ))

dc_tap_and_training <- dc_tap_and_training %>%
  mutate(plotting_enhancer_categories = factor(
    plotting_enhancer_categories,
    levels = c(
      "Weakly acetylated elements (Low H3K27ac, no CTCF)",
      "Strongly acetylated elements (High H3K27ac)",
      "CTCF elements (CTCF, low H3K27ac)",
      "Polycomb elements (H3K27me3)"
    )
  ))


### CREATE THE PROPORTIONS PLOT ==============================================

# Function to create plot with new enhancer categories
create_plot <- function(combined_data, title_suffix) {
  # First calculate total counts per dataset
  dataset_totals <- combined_data %>%
    group_by(Source, Dataset) %>%
    summarise(TotalCount = n())
  
  # Calculate proportions and prepare data for plotting
  prop_data <- combined_data %>%
    group_by(Source, Dataset, plotting_enhancer_categories) %>%
    summarise(Count = n()) %>%
    left_join(dataset_totals, by = c("Source", "Dataset")) %>%
    mutate(
      Proportion = Count / TotalCount,
      y_position = cumsum(Proportion) - Proportion / 2
    )
  
  # Create stacked bar plot
  p <- ggplot(prop_data, aes(x = Dataset, y = Proportion, fill = plotting_enhancer_categories)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5),
              color = "black", size = 2.5) +
    geom_text(aes(x = Dataset, y = -0.05, label = paste0("(", TotalCount, ")")),
              size = 3, vjust = 0) +
    labs(
      title = "Proportion of Enhancer-Gene Classes Across Datasets",
      subtitle = title_suffix,
      x = "Dataset",
      y = "Proportion of EG Pairs",
      fill = "Element Category"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "lightgray"),
      panel.grid = element_blank()
    ) +
    scale_fill_manual(values = c(
      "Weakly acetylated elements (Low H3K27ac, no CTCF)" = "#F4631E",
      "Strongly acetylated elements (High H3K27ac)"        = "#CB041F",
      "CTCF elements (CTCF, low H3K27ac)"                  = "#309898",
      "Polycomb elements (H3K27me3)"                       = "#FF9F00"
    ))
  
  return(p)
}

# Process the data for plotting
dc_tap_paper_categories <- dc_tap_and_training %>%
  mutate(Source = "All") %>%
  mutate(Dataset = case_when(
    Dataset == "Unfilt. TAPSeq" ~ "K562 DC TAP Seq",
    Dataset == "WTC11_DC_TAP_Seq" ~ "WTC11 DC TAP Seq",
    Dataset == "Gasperini2019" ~ "Gasperini et al., 2019",
    Dataset == "FlowFISH_K562" ~ "K562 FlowFISH",
    TRUE ~ Dataset
  ))

# Set the desired order for datasets
dataset_order <- c("Gasperini et al., 2019", "K562 FlowFISH", "K562 DC TAP Seq", "WTC11 DC TAP Seq")

# Process the data and set factor levels for ordering
dc_tap_paper_categories <- dc_tap_paper_categories %>%
  mutate(Dataset = factor(Dataset, levels = dataset_order))

# Create the plot
proportions_plot_all <- create_plot(dc_tap_paper_categories, "All Pairs")

# Create the plot for Regulated == TRUE pairs
proportions_plot_pos <- create_plot(dc_tap_paper_categories %>% filter(Regulated == TRUE), "CRISPRi Positives")

# Save the proportion plots
ggsave(plot = proportions_plot_all, filename = snakemake@output$proportions_plot_all, device = "pdf", height = 5, width = 7)
ggsave(plot = proportions_plot_pos, filename = snakemake@output$proportions_plot_pos, device = "pdf", height = 5, width = 7)


### FOLD CHANGE PLOT =========================================================

# Calculate proportions and fold changes for each dataset and category
enhancer_percentages <- dc_tap_paper_categories %>%
  group_by(Dataset) %>%
  mutate(
    total_dataset_pairs = n(),
    total_dataset_positives = sum(Regulated == TRUE)
  ) %>%
  group_by(Dataset, plotting_enhancer_categories) %>%
  summarise(
    category_pairs = n(),
    category_positives = sum(Regulated == TRUE),
    total_dataset_pairs = dplyr::first(total_dataset_pairs),
    total_dataset_positives = dplyr::first(total_dataset_positives),
    .groups = 'drop'
  ) %>%
  mutate(
    # Calculate proportions
    proportion_of_all_pairs = category_pairs / total_dataset_pairs,
    proportion_of_positives = category_positives / total_dataset_positives,
    # Calculate fold change
    fold_change = (proportion_of_positives / proportion_of_all_pairs) * 100
  )

# Create the barplot
fold_change_plot <- ggplot(enhancer_percentages, 
                           aes(x = fold_change, 
                               y = plotting_enhancer_categories, 
                               fill = Dataset)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sprintf("%.1f", fold_change)), 
            position = position_dodge(width = 0.9),
            hjust = -0.1) +
  labs(
    title = "Enhancer Category Enrichment in Positive Pairs",
    x = "Fold-change (% positives / % tested)",
    y = "Element category",
    fill = "Dataset"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set2") +
  xlim(0, max(enhancer_percentages$fold_change) * 1.2) 

# Save the fold change plot
ggsave(plot = fold_change_plot, filename = snakemake@output$fold_change_plot, device = "pdf", height = 5, width = 9)


### SIGNIFICANCE TESTING ======================================================

# Perform chi-square test between datasets for each category
datasets <- c("Gasperini et al., 2019", "K562 FlowFISH", "K562 DC TAP Seq", "WTC11 DC TAP Seq")
combinations <- combn(datasets, 2, simplify = FALSE)
comparison_names <- c("Gas v Flow", "Gas v K562", "Gas v wt11", 
                      "Flow v K562", "Flow v wt11", "K562 v wt11")

# Function for chi-square test
chi_square_test <- function(data, d1, d2, cat) {
  data1 <- filter(data, Dataset == d1, plotting_enhancer_categories == cat)
  data2 <- filter(data, Dataset == d2, plotting_enhancer_categories == cat)
  
  matrix <- matrix(c(data1$category_positives, 
                     data1$category_pairs - data1$category_positives,
                     data2$category_positives, 
                     data2$category_pairs - data2$category_positives), 
                   nrow = 2)
  
  tryCatch(chisq.test(matrix)$p.value, error = function(e) NA)
}

# Calculate all chi-square tests
results <- crossing(
  category = unique(enhancer_percentages$plotting_enhancer_categories),
  pair = seq_along(combinations)
) %>%
  mutate(
    dataset1 = map_chr(pair, ~combinations[[.x]][1]),
    dataset2 = map_chr(pair, ~combinations[[.x]][2]),
    comparison = comparison_names[pair],
    p_value = pmap_dbl(list(dataset1, dataset2, category), 
                       ~chi_square_test(enhancer_percentages, ..1, ..2, ..3))
  ) %>%
  arrange(category, comparison)

# Add proportions to results
results_with_props <- results %>%
  left_join(
    select(enhancer_percentages, 
           dataset1 = "Dataset", 
           category = "plotting_enhancer_categories",
           prop1_pos = "proportion_of_positives",
           prop1_all = "proportion_of_all_pairs"),
    by = c("dataset1", "category")
  ) %>%
  left_join(
    select(enhancer_percentages, 
           dataset2 = "Dataset", 
           category = "plotting_enhancer_categories",
           prop2_pos = "proportion_of_positives",
           prop2_all = "proportion_of_all_pairs"),
    by = c("dataset2", "category")
  )

# Format results
formatted_results <- results_with_props %>%
  mutate(across(starts_with("prop"), ~sprintf("%.1f%%", . * 100)),
         p_value = sprintf("%.3e", p_value)) %>%
  arrange(category, comparison)

# Print significant results
cat("\nSignificant results (p < 0.05):\n")
print(filter(results_with_props, p_value < 0.05) %>%
        select(category, comparison, p_value))

cat("\nHighly significant results (p < 0.001):\n")
print(filter(results_with_props, p_value < 0.001) %>%
        select(category, comparison, p_value))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(formatted_results, snakemake@output$chi_squared_results)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)