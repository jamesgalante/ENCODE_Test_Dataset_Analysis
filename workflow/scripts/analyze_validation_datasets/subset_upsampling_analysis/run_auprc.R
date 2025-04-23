# Script: run_auprc.R

### SETUP =====================================================================

# Saving image for debugging
# Check if the directory exists, and create it only if it doesn't
if (!file.exists("RDA_objects/run_auprc")) { dir.create("RDA_objects/run_auprc", recursive = TRUE) }

# Check if replicate exists in wildcards
if ("replicate" %in% names(snakemake@wildcards)) {
  # Case for replicated analyses
  save.image(paste0("RDA_objects/run_auprc/", 
                    snakemake@wildcards$dataset_type, "_",
                    snakemake@wildcards$subset, "_rep",
                    snakemake@wildcards$replicate, ".rda"))
} else {
  # Case for original (non-replicated) analyses
  save.image(paste0("RDA_objects/run_auprc/original_",
                    snakemake@wildcards$subset, ".rda"))
}

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
  library(boot)
  library(ROCR)
  library(caTools)
  library(data.table)
  library(BiocParallel)
  library(ggpubr)
  library(ggcorrplot)
  source(file.path(snakemake@scriptdir, "auprc_helper_functions/crisprComparisonPlotFunctions.R"))
  source(file.path(snakemake@scriptdir, "auprc_helper_functions/crisprComparisonBootstrapFunctions.R"))
  source(file.path(snakemake@scriptdir, "auprc_helper_functions/functions_for_AUPRC.R"))
})

message("Loading input files")
dataset <- readRDS(snakemake@input$dataset)
pred_config <- fread(snakemake@params$pred_config, colClasses = c("alpha" = "numeric", "color" = "character"))


### RUN AUPRC ==========================================================

# Run AUPRC on the dataset
dataset_results <- auprc_calc(dataset, pred_config)


### SAVE OUTPUT ===============================================================

# Save the AUPRC results
saveRDS(dataset_results, snakemake@output$auprc_results)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

