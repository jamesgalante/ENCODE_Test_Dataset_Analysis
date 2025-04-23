# Script: # Script: sceptre_differential_expression.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/sceptre_differential_expression")) { dir.create("RDA_objects/sceptre_differential_expression", recursive = TRUE) }
save.image(paste0("RDA_objects/sceptre_differential_expression/", snakemake@wildcards$sample, ".rda"))
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
  library(sceptre)
})

message("Loading input files")
sceptre_object <- readRDS(snakemake@input$sceptre_diffex_input)


### RUN DIFFEX ================================================================

sceptre_object <- sceptre_object %>% 
  assign_grnas(method = "mixture") %>%
  run_qc() %>%
  run_discovery_analysis()


### SAVE OUTPUT ===============================================================

# Create "plots" subdirectory inside dirname(snakemake@output$discovery_results)
plots_dir <- file.path(dirname(snakemake@output$discovery_results), "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Define ggsave wrapper for consistent usage
save_plot <- function(plot, file_name) {
  file_path <- file.path(plots_dir, paste0(file_name, ".pdf"))
  ggsave(
    filename = file_path,
    plot = plot,
    device = "pdf",
    width = 5,
    height = 5
  )
}

# Save plots with their respective names
message("Saving plots")
save_plot(plot_grna_count_distributions(sceptre_object), "plot_grna_count_distributions")
save_plot(plot_assign_grnas(sceptre_object), "plot_assign_grnas")
# save_plot(plot_run_calibration_check(sceptre_object), "plot_run_calibration_check")
save_plot(plot_run_discovery_analysis(sceptre_object), "plot_run_discovery_analysis")
save_plot(plot_covariates(sceptre_object), "plot_covariates")
save_plot(plot_run_qc(sceptre_object), "plot_run_qc")

# Write all outputs to directory
message("Saving output files")
write_outputs_to_directory(sceptre_object, dirname(snakemake@output$discovery_results))

# Save the final sceptre object
saveRDS(sceptre_object, snakemake@output$final_sceptre_object)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)