# Script: # Script: chi_squared_enhancer_classes.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/chi_squared_enhancer_classes.rda"))
message("Saved Image")
stop("Manually Stopped Program after Saving Image")

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
expt_pred_merged_annot_files <- lapply(snakemake@input, read_tsv)
names(expt_pred_merged_annot_files) <- c("validation", "training", "unfiltered_k562_dc_tap")

val <- expt_pred_merged_annot_files[["validation"]]
train <- expt_pred_merged_annot_files[["training"]]


### SET UP CHI SQUARED TEST ===================================================

# Get counts of each enhancer class in both datasets
training_counts <- train %>%
  group_by(Enhancer_Class) %>%
  summarise(training_n = n()) %>%
  arrange(Enhancer_Class)

validation_counts <- val %>%
  group_by(Enhancer_Class) %>%
  summarise(validation_n = n()) %>%
  arrange(Enhancer_Class)

# Combine counts into one dataframe
combined_counts <- full_join(training_counts, validation_counts, by = "Enhancer_Class") %>%
  replace_na(list(training_n = 0, validation_n = 0))

# Create vectors for chi-squared test
observed <- combined_counts$validation_n
expected <- combined_counts$training_n

# Calculate proportions for expected (training) values
# Chi-squared test expects raw counts but scaled to match total observations
total_validation <- sum(observed)
total_training <- sum(expected)
expected_scaled <- expected * (total_validation/total_training)

# Perform chi-squared test
chi_squared_result <- chisq.test(observed, p = expected/sum(expected))

# Create summary dataframe with counts and proportions
summary_df <- combined_counts %>%
  mutate(
    training_prop = training_n/sum(training_n),
    validation_prop = validation_n/sum(validation_n),
    prop_difference = validation_prop - training_prop
  )

# Save results
results <- list(
  chi_squared_test = chi_squared_result,
  summary_table = summary_df,
  test_details = list(
    degrees_of_freedom = chi_squared_result$parameter,
    p_value = chi_squared_result$p.value,
    chi_squared_statistic = chi_squared_result$statistic
  )
)


### CALCULATE INDIVIDUAL CATEGORY SIGNIFICANCE ================================

# Function to perform prop test for each category
prop_test_by_category <- function(cat_data) {
  # Create contingency table for this category vs others
  cat_matrix <- matrix(
    c(cat_data$validation_n, 
      sum(validation_counts$validation_n) - cat_data$validation_n,
      cat_data$training_n, 
      sum(training_counts$training_n) - cat_data$training_n),
    nrow = 2
  )
  # Perform test
  test <- prop.test(cat_matrix)
  # Return p-value
  return(test$p.value)
}

# Add p-values to summary
summary_df <- summary_df %>%
  rowwise() %>%
  mutate(
    p_value = prop_test_by_category(cur_data()),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

### CREATE VISUALIZATION ====================================================

# Create plot
prop_plot <- ggplot(summary_df, aes(x = Enhancer_Class, y = prop_difference * 100)) +
  geom_bar(stat = "identity", fill = ifelse(summary_df$prop_difference > 0, "#F28E2B", "#4E79A7")) +
  geom_text(aes(label = significance, y = ifelse(prop_difference > 0, prop_difference * 100 + 1, prop_difference * 100 - 1)),
            size = 5) +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Enhancer Class",
    y = "Difference in Proportion (Validation - Training) %",
    title = "Differences in Enhancer Class Proportions",
    subtitle = paste("Chi-squared test p-value:", format.pval(chi_squared_result$p.value, digits = 3))
  ) +
  theme(
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(
  filename = snakemake@output$chi_squared_plot,
  plot = prop_plot,
  width = 6,
  height = 4
)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)