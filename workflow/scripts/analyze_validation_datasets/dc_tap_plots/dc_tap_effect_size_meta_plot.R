# Script: # Script: dc_tap_effect_size_meta_plot.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/dc_tap_effect_size_meta_plot.rda"))
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
combined_unfiltered_k562_dc_tap <- read.table(snakemake@input$combined_unfiltered_k562_dc_tap)
per_guide_effect_sizes_unfiltered_k562_dc_tap <- read.table(snakemake@input$per_guide_effect_sizes_unfiltered_k562_dc_tap)
guide_targets <- read.table(snakemake@input$guide_targets)


### x =========================================================

# Function to process data and create meta-plot
create_guide_meta_plot <- function(guide_targets, effect_sizes) {
  # Join guide information with effect sizes
  combined_data <- guide_targets %>%
    left_join(effect_sizes, by = c("name" = "grna_id")) %>%
    # Calculate center of target regions and guides and absolute effect size
    mutate(
      target_center = (target_start + target_end) / 2,
      guide_center = (start + end) / 2,
      # Calculate relative position from target center
      relative_position = guide_center - target_center,
      # Take absolute value of effect size
      abs_effect = abs(log_2_fold_change),
      # Convert significant to factor with better labels
      significant = factor(significant, 
                           levels = c(TRUE, FALSE),
                           labels = c("Significant", "Not Significant"))
    ) %>%
    # Remove NA values
    filter(!is.na(abs_effect))
  
  # Create binned averages to reduce data size
  binned_data <- combined_data %>%
    group_by(significant) %>%
    mutate(
      # Create bins of relative position
      pos_bin = cut(relative_position, breaks = 50)
    ) %>%
    group_by(significant, pos_bin) %>%
    summarise(
      mean_effect = mean(abs_effect),
      se_effect = sd(abs_effect) / sqrt(n()),
      mid_pos = mean(relative_position),
      .groups = 'drop'
    )
  
  # Create the meta-plot using binned data
  p <- ggplot(binned_data, aes(x = mid_pos, y = mean_effect, color = significant)) +
    # Add line
    geom_line() +
    # Add error bands
    geom_ribbon(aes(ymin = mean_effect - se_effect, 
                    ymax = mean_effect + se_effect,
                    fill = significant), 
                alpha = 0.2,
                color = NA) +
    # Customize appearance
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray50")) +
    scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "gray50")) +
    labs(
      title = "Guide Effect Size Meta-plot",
      x = "Distance from Target Center (bp)",
      y = "|Log2 Fold Change|",
      color = "Significance",
      fill = "Significance"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  return(list(
    plot = p,
    processed_data = combined_data,
    binned_data = binned_data
  ))
}

# Example usage:
result <- create_guide_meta_plot(guide_targets, effect_sizes)
print(result$plot)



### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)