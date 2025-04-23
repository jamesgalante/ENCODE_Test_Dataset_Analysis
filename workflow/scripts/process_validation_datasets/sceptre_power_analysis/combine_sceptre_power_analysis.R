# Script: # Script: combine_sceptre_power_analysis.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/combine_sceptre_power_analysis")) { dir.create("RDA_objects/combine_sceptre_power_analysis", recursive = TRUE) }
save.image(paste0("RDA_objects/combine_sceptre_power_analysis/", snakemake@wildcards$sample, "_", snakemake@wildcards$effect_size, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(tidyverse)
})

### COMBINE SIMULATION OUTPUTS ==============================================

# Combine each file into one tsv and save it
message("Combining power analysis outputs...")

# Initialize an empty list to store data frames
data_list <- list()

# Loop over input files
message("Gathering inputs")
for (file_path in snakemake@input$splits) {
  # Read the current file
  current_data <- read_tsv(file_path, col_types = cols())
  
  # Append the data frame to the list
  data_list[[length(data_list) + 1]] <- current_data
}

# Combine all data frames into one
message("Combining all data frames into one")
combined_data <- do.call(rbind, data_list)

### SAVE COMBINE OUTPUTS ===================================================

# Write the combined data frame to the output file
write_tsv(combined_data, snakemake@output$combined_power_analysis_output)
message("Power analysis outputs combined successfully.")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)