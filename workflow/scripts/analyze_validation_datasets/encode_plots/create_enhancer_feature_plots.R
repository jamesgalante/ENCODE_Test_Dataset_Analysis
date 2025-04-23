# Script: # Script: create_enhancer_feature_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_enhancer_feature_plots.rda"))
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


### ORGANIZE DATA =============================================================

# Create combined versions but keep them in their respective Training/Validation sources
training_combined <- filt_combined_training %>%
  mutate(
    Source = "Training",  # Same source as original training data
    category = "Combined Training Data\n[All Training Datasets]"  # New category that will appear alongside others
  )

validation_combined <- filt_combined_validation %>%
  mutate(
    Source = "Validation",  # Same source as original validation data
    category = "Combined Heldout Data\n[All Heldout Datasets]"  # New category that will appear alongside others
  )

# Original datasets (no need to modify source)
training_split <- filt_combined_training %>% 
  mutate(Source = "Training")

validation_split <- filt_combined_validation %>%
  mutate(Source = "Validation")

# Bind all datasets together
all_pairs_with_combined <- rbind(
  training_split,
  training_combined,
  validation_split,
  validation_combined
)

# Convert DistToTSS from bp to Mb
all_pairs_with_combined <- all_pairs_with_combined %>%
  mutate(DistToTSS_Mb = DistToTSS / 1e6)

# Pivot to long format - now Source will only have "Training" and "Validation" levels
all_long <- all_pairs_with_combined %>%
  pivot_longer(
    cols = c(dhs_signal, h3k27ac_signal, DistToTSS_Mb), 
    names_to = "type", 
    values_to = "value"
  ) %>%
  mutate(Source = factor(Source, levels = c("Training", "Validation"))) %>%
  mutate(type = factor(
    type, 
    levels = c("dhs_signal", "h3k27ac_signal", "DistToTSS_Mb"), 
    labels = c("DNase-seq", "H3K27ac", "Distance to TSS (Mb)"))
  )


### REFORMATTING CATEGORY NAMES ===============================================

# For the paper figure, let's rename each of the categories:
# Reference (cell type)

# Create a named vector for category renaming
category_rename <- c(
  "K562 FlowFISH" = "Nasser et al., 2021 (K562)",
  "Gasperini et al., 2019" = "Gasperini et al., 2019 (K562)",
  "Schraivogel et al., 2020" = "Schraivogel et al., 2020 (K562)",
  "Combined Training Data\n[All Training Datasets]" = "Combined Training Data\n[All Training Datasets]",
  "GM12878 Nasser et al., 2021" = "Nasser et al., 2021 (GM12878)",
  "HCT116 Guckelberger et al., inPrep" = "Guckelberger et al., 2024 (HCT116)",
  "Jurkat Nasser et al., 2021" = "Nasser et al., 2021 (Jurkat)",
  "Klann et al., 2021" = "Klann et al., 2021 (K562)",  # Need cell type
  "Xie et al., 2019" = "Xie et al., 2019 (K562)",      # Need cell type
  "Morris et al., 2023" = "Morris et al., 2023 (K562)", # Need cell type
  "K562 DC TAP Seq" = "DC TAP-Seq in Prep (K562)",
  "Reilly et al., 2021" = "Reilly et al., 2021 (K562)", # Need cell type
  "WTC11 DC TAP Seq" = "DC TAP-Seq in Prep (WTC11)",
  "Combined Heldout Data\n[All Heldout Datasets]" = "Combined Heldout Data\n[All Heldout Datasets]"
)

# Apply the renaming in your data processing pipeline
all_long <- all_long %>%
  mutate(category = recode(category, !!!category_rename))


### PLOT DATA =================================================================

p_all <- ggplot(all_long, aes(x = category, y = value)) +
  geom_violin(aes(fill = type), trim = FALSE, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, fill = NA, color = "black") +
  coord_flip() +
  facet_grid(Source ~ type, scales = "free_y", space = "free_y") +
  scale_y_log10(
    breaks = c(1e-2, 1e0, 1e2), 
    labels = c("-2", "0", "2")
  ) +
  labs(y = "Value (log10)", x = NULL) +
  scale_fill_manual(values = c("DNase-seq" = "#4783B5",
                               "H3K27ac" = "#FAA51A",
                               "Distance to TSS (Mb)" = "grey")) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold", angle = 270),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = p_all, 
       filename = snakemake@output$properties_plot, 
       device = "pdf", 
       height = 6, 
       width = 8)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)