# Script: create_sce_object.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/create_sce_object")) { dir.create("RDA_objects/create_sce_object", recursive = TRUE) }
save.image(paste0("RDA_objects/create_sce_object/", snakemake@wildcards$sample, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
# required packages and functions
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SingleCellExperiment)
  library(rtracklayer)
  library(readr)
  library(DESeq2)
})
message("Loading input files")

# column types in guide targets file
guide_targets_cols <- cols(
  .default = col_character(),  # for target_type columns, which is optional as of now
  chr = col_character(),
  start = col_integer(),
  end = col_integer(),
  name = col_character(),
  strand = col_character(),
  spacer = col_character(),
  target_chr = col_character(),
  target_start = col_integer(),
  target_end = col_integer(),
  target_name = col_character(),
  target_strand = col_character()
)

# load guide targets
guide_targets <- read_tsv(snakemake@input$guide_targets, col_types = guide_targets_cols)

# load in the final sceptre object
final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object)


### CREATE SCE OBJECT ========================================================

expr <- final_sceptre_object@response_matrix[[1]]

if ("batch" %in% colnames(final_sceptre_object@covariate_data_frame)) {
  cell_metadata <- data.frame(
    cell_barcode = rownames(final_sceptre_object@covariate_data_frame),
    cell_batches = final_sceptre_object@covariate_data_frame$batch
  )
  sce <- SingleCellExperiment(assays = list(counts = expr), colData = cell_metadata)
  colnames(sce) <-  rownames(final_sceptre_object@covariate_data_frame)
} else {
  cell_metadata <- NULL
  sce <- SingleCellExperiment(assays = list(counts = expr))
  colnames(sce) <-  rownames(final_sceptre_object@covariate_data_frame)
}


### ADD SCEPTRE GRNA_PERTS ASSIGNMENTS TO SCE =================================

# Add the individual grna perts first
message("Adding the individual grna_perts")
get_grna_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@initial_grna_assignment_list)
}

individual_grna_assignments <- get_grna_assignments(final_sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_grna_assignments)
nCols <- length(colnames(sce))

# Initialize an empty sparse matrix
sparseMatrix <- Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)

# Assuming sparseMatrix is already initialized correctly
message("Looping through grna_assignments")
for (i in seq_along(individual_grna_assignments)) {
  indices <- individual_grna_assignments[[i]]
  if(length(indices) > 0) { # Check if there are any indices to assign
    sparseMatrix[i, indices] <- 1
  } else {
    # Optionally, handle the case where there are no indices
    # For example, by doing nothing, or logging this case
    cat(sprintf("No indices for row %d\n", i))
  }
}

# Set column names
colnames(sparseMatrix) <- colnames(sce)
rownames(sparseMatrix) <- names(individual_grna_assignments)

# Add to sce object
message("Adding grna_perts to sce object")
altExp(sce, e = "grna_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))


### ADD SCEPTRE CRE_PERTS ASSIGNMENTS TO SCE ==================================

# Now add the individual "cre_perts"
message("Adding the individual cre_perts")
get_cre_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@grna_assignments$grna_group_idxs)
}

individual_cre_assignments <- get_cre_assignments(final_sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_cre_assignments)
nCols <- length(colnames(sce))

# Initialize an empty sparse matrix
sparseMatrix <- Matrix(0, nrow = nRows, ncol = nCols, sparse = TRUE)

# Assuming sparseMatrix is already initialized correctly
message("Looping through cre_assignments")
for (i in seq_along(individual_cre_assignments)) {
  print(i)
  # Indices from individual_cre_assignments point to cells_in_use
  indices_to_cells_in_use <- individual_cre_assignments[[i]]
  
  # Now retrieve the actual indices from sceptre_object@cells_in_use
  indices <- final_sceptre_object@cells_in_use[indices_to_cells_in_use]
  
  if(length(indices) > 0) { # Check if there are any indices to assign
    sparseMatrix[i, indices] <- 1
  } else {
    # Optionally, handle the case where there are no indices
    # For example, by doing nothing, or logging this case
    cat(sprintf("No indices for row %d\n", i))
  }
}

# Set column names
colnames(sparseMatrix) <- colnames(sce)
rownames(sparseMatrix) <- names(individual_cre_assignments)

# Add to sce object
message("Adding cre_prets to sce object")
altExp(sce, e = "cre_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))


### ADD SCEPTRE DISPERION ESTIMATES ===========================================

dispersion_values <- lapply(final_sceptre_object@response_precomputations, function(x) 1/x$theta)

# Add to rowData
rowData(sce)$dispersion <- dispersion_values[rownames(sce)]
rowData(sce)$mean <- rowMeans(expr)


### CALCULATE SIZE FACTORS ====================================================

# Calculate total_umis and detected_genes for Deseq2 object creation
response_matrix <- assay(sce, "counts")
coldata <- data.frame(
  total_umis = colSums(response_matrix),
  detected_genes = colSums(response_matrix > 0)
)


# create DESeq2 object containing count data
dds <- DESeqDataSetFromMatrix(countData = response_matrix, 
                              colData = coldata, 
                              design = ~ 1)
  
# compute size factors (removed lib size option for simplicity)
dds <- estimateSizeFactors(dds, type = "poscounts")
colData(sce)[, "size_factors"] <- sizeFactors(dds)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(sce, snakemake@output$perturb_sce)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)