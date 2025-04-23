# Script: # Script: compute_power_from_simulations.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/compute_power_from_simulations")) { dir.create("RDA_objects/compute_power_from_simulations", recursive = TRUE) }
save.image(paste0("RDA_objects/compute_power_from_simulations/", snakemake@wildcards$sample, "_", snakemake@wildcards$effect_size, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# required packages and functions
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load input file
combined_power_analysis_output <- read_tsv(snakemake@input$combined_power_analysis_output)
discovery_results <- readRDS(snakemake@input$discovery_results)


### GET NOMINAL P VALUE ======================================================

# From the discovery results, retrieve the largest p value that is significant after FDR correction
max_nom_p_val <- discovery_results %>%  
  filter(significant == TRUE) %>% 
  pull(p_value) %>% 
  max()


### CALCULATING POWER ======================================================

power <- combined_power_analysis_output %>%
  filter(!is.na(log_2_fold_change)) %>%
  group_by(rep) %>%
  group_by(grna_target, response_id) %>%
  summarize(
    power = mean(p_value < max_nom_p_val & log_2_fold_change < 0), 
    mean_log_2_fold_change = mean(log_2_fold_change),
    mean_pert_cells = mean(num_pert_cells),
    .groups = "drop") %>%
  arrange(desc(power))
  

# Add the effect size (NOTE: if adding ability to have multiple effect sizes in one run, must change this to a wildcard)
power$effect_size <- as.numeric(snakemake@wildcards$effect_size)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(power, file = snakemake@output$power_analysis_results)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)