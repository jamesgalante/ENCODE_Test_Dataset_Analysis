# Script: # Script: fetch_encode_features.R

### SETUP =====================================================================

# Get the dataset wildcard
dataset_name <- snakemake@wildcards$dataset

# Saving image for debugging
save.image(paste0("RDA_objects/fetch_encode_features_", dataset_name, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Home folder for crisprComparison helper functions in:
crispr_comparison_scripts <- paste0(sub("/scripts/.*", "/scripts/", snakemake@scriptdir), "benchmark_validation_datasets/crispr_benchmarking")

message("Loading in packages")
suppressPackageStartupMessages({
  source(file.path(crispr_comparison_scripts, "crisprComparisonLoadInputData.R"))
  source(file.path(crispr_comparison_scripts, "crisprComparisonMergeFunctions.R"))
  source(file.path(crispr_comparison_scripts, "crisprComparisonSimplePredictors.R"))
})


### LOADING FILES ============================================================

# config entry for this comparison is used to load named list of input files
config <- snakemake@config$benchmark_validation_datasets$crispr_benchmarking$comparisons[[snakemake@wildcards$dataset]]

# load pred_config file
include_col <- ifelse(is.null(snakemake@params$include_col), "include", snakemake@params$include_col)
pred_config <- importPredConfig(snakemake@input$pred_config,
                                expr = !is.null(snakemake@input$expressed_genes),
                                include_col = include_col,
                                filter = snakemake@params$filter_include_col)

# load experimental data
message("Reading experimental data in: ", snakemake@input$experiment)
expt <- fread(file = snakemake@input$experiment, showProgress = FALSE,
              colClasses = c("ValidConnection" = "character"))
message("\tLoaded experimental data with ", nrow(expt), " rows.\n")

# load all prediction files
pred_list <- loadPredictions(config$pred, show_progress = FALSE)

# load tss and gene universe files
tss_annot <- fread(snakemake@input$tss_universe, select = 1:6,
                   col.names = c("chrTSS", "startTSS", "endTSS", "gene", "score", "strandTSS"))
gene_annot <- fread(snakemake@input$gene_universe, select = 1:6,
                    col.names = c("chr", "start", "end", "gene", "score", "strand"))

# load optional cell mapping files if provided
ct_map_files <- config$cell_type_mapping
if (!is.null(ct_map_files)) {
  cell_mappings <- lapply(ct_map_files, FUN = fread)
  qcCellMapping(cell_mappings)
} else {
  cell_mappings <- list()
}

# load expressed genes files if provided
if (!is.null(snakemake@input$expressed_genes)) {
  expressed_genes <- loadGeneExpr(snakemake@input$expressed_genes)
} else {
  expressed_genes <- NULL
}

# QC pred_config file
qcPredConfig(pred_config, pred_list = pred_list)

# QC predictions and experimental data
pred_list <- qcPredictions(pred_list, pred_config = pred_config, one_tss = FALSE)
expt <- qcExperiment(expt, pos_col = snakemake@params$pos_col, remove_na_pos = TRUE)


### PROCESS INPUT DATA ========================================================

# base output directory for any output
outdir <- dirname(snakemake@output$merged)

# filter experimental data for genes in gene universe
missing_file <- file.path(outdir, "expt_missing_from_gene_universe.txt")
expt <- filterExptGeneUniverse(expt, genes = tss_annot, missing_file = missing_file)

# add expression information to experimental data if specified
if (!is.null(snakemake@input$expressed_genes)) {
  expt <- addGeneExpression(expt, expressed_genes = expressed_genes)
}

# cell type matching and filter predictions for cell types in experimental data
message("Mapping cell types in predictions to cell types in experimental data")
pred_list <- mapCellTypes(pred_list, cell_mappings = cell_mappings)
pred_list <- lapply(pred_list, FUN = function(p) p[p$ExperimentCellType %in% expt$CellType, ] )

# verify if bad cell matching resulted in no data for some predictions after matching
pred_rows <- vapply(pred_list, FUN = nrow, FUN.VALUE = integer(1))
if (any(pred_rows == 0)) {
  stop("No predictions left for ", paste(names(pred_rows[pred_rows == 0]), collapse = ", "),
       " after cell type matching. Check that cell type mapping files are correct.", call. = FALSE)
}


### OVERLAP EXPERIMENTAL DATA WITH PREDICTIONS ================================

# check if genes in experimental data are also found in predictions and write to file
# TODO: make this per cell type
genes_summary_file <- file.path(outdir, "experimental_genes_in_predictions.txt")
checkExistenceOfExperimentalGenesInPredictions(expt, pred_list, summary_file = genes_summary_file)

# merge experimental data with predictions
message("\nMerging experimentals data and predictions:")

# This is where the current file diverges from the normal `mergePredictionsWithExperiment.R` 
# I add an `extract_features` parameter to include the feature columns after merging
merged <- combineAllExptPred(expt = expt, 
                             pred_list = pred_list["ENCODE_rE2G"],
                             config = pred_config,
                             outdir = outdir,
                             fill_pred_na = TRUE,
                             extract_features = TRUE)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
readr::write_tsv(merged, file = snakemake@output$merged)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)