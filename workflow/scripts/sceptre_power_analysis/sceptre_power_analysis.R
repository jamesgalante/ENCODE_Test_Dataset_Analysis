# Script: sceptre_power_analysis.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/sceptre_power_analysis")) { dir.create("RDA_objects/sceptre_power_analysis", recursive = TRUE) }
save.image(paste0("RDA_objects/sceptre_power_analysis/", snakemake@wildcards$sample, "_", snakemake@wildcards$effect_size, "_", snakemake@wildcards$split, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING PACKAGES ==========================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(SingleCellExperiment)
  library(stringr)
  library(sceptre)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
})


### LOADING FILES =============================================================

# Read the input arguments into variables
message("Reading in snakemake variables")
# convert 'percentage decrease' effect size to 'relative expression level'
effect_size <- 1 - as.numeric(snakemake@wildcards$effect_size)
guide_sd <- snakemake@params$guide_sd
reps <- snakemake@params$reps
n_ctrl <- FALSE

# Read in the "guide_file"
final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object)
discovery_pairs_split <- read_tsv(snakemake@input$gene_gRNA_group_pairs_split, col_names = c("grna_group", "response_id"))
gRNA_groups_table <- readRDS(snakemake@input$gRNA_groups_table)
sce <- readRDS(snakemake@input$perturb_sce)


### SET UP STUFF ==============================================================

# Get all the perturbation targets in the subsetted guide_file that are represented in the sceptre_object and where at least one pair passes qc
discovery_pairs_which_pass_qc <- final_sceptre_object@discovery_pairs_with_info[final_sceptre_object@discovery_pairs_with_info$pass_qc == TRUE,]
perts <- intersect(discovery_pairs_split$grna_group, discovery_pairs_which_pass_qc$grna_group)

# Set up an empty datafame to store the results
discovery_results <- data.frame()


###  Doing Thing =========================================================

es_track <- list()
# Runs each power analysis rep for one perturbation
for (pert in perts){
  
  # Get a list of the relevant discovery pairs
  discovery_relevant_pairs_pert <- discovery_pairs_which_pass_qc[discovery_pairs_which_pass_qc$grna_group == pert,]
  
  # Get all the guides that target the current `pert`
  pert_guides <- gRNA_groups_table %>%
    filter(grna_target == pert) %>%
    pull(grna_id)  
  
  
  # Assigning appropriate sampling function based on status of n_ctrl
  message("Assigning pert_input_function with n_ctrl value")
  if (is.numeric(n_ctrl)) {
    pert_object <- pert_input_sampled(pert, sce = sce, pert_level = "cre_perts", n_ctrl = n_ctrl, cell_batches = "cell_batches")
  } else if (n_ctrl == FALSE) {
    pert_object <- pert_input(pert, sce = sce, pert_level = "cre_perts")
  } 
  
  # get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  # Convert to a sparse matrix, so the sampling function works in `create_guide_pert_status`
  grna_perts <- as(grna_perts, "CsparseMatrix") 
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts, pert_guides = pert_guides)
  
  # Subset the pert_object for only those genes which are tested against the current pert
  pert_genes <- discovery_relevant_pairs_pert$response_id
  pert_object <- pert_object[pert_genes,]
  
  # Create effect size matrix (sampled from negative binomial distribution around effect_size or 1)
  effect_sizes <- structure(rep(effect_size, nrow(pert_object)), names = rownames(pert_object))
  
  for (rep in seq(reps)) {
    
    # Create and center effect size matrices
    es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                        gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
    es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status, gene_effect_sizes = effect_sizes)
    es_mat_use <- es_mat[, colnames(assay(pert_object, "counts"))]
    
    # Simulate Counts
    message("Simulating Counts")
    sim_counts <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat_use)
    
    # c <- rowMeans(assay(sim_counts, "counts")[, pert_status == 0])
    # p <- rowMeans(assay(sim_counts, "counts")[, pert_status == 1])
    # m <- mean(p/c)
    # es_track <- append(es_track, m)
    
    # Now let's set up the discovery_analysis
    message("Setting up sceptre object for Disovery Analysis")
    full_response_matrix_sim_sparse <-  as(assay(sim_counts, "counts"), "RsparseMatrix")
    sceptre_object_use <- final_sceptre_object
    
    message("The class of what's being assigned to the response matrix")
    sceptre_object_use@response_matrix <- list(full_response_matrix_sim_sparse)
    sceptre_object_use@discovery_pairs_with_info <- discovery_relevant_pairs_pert
    
    
    # Fix the `cells_in_use` parameter for indexing when the n_ctrl is not FALSE
    if (is.numeric(n_ctrl)) {
      sceptre_object_use@cells_in_use <- seq(dim(full_response_matrix_sim_sparse)[[2]])
      
      # Only keep the covariate rows that are in the colnames of `full_response_matrix_sim_sparse`
      sceptre_object_use@covariate_data_frame <- sceptre_object_use@covariate_data_frame[rownames(sceptre_object_use@covariate_data_frame) %in% colnames(full_response_matrix_sim_sparse),]
      sceptre_object_use@covariate_matrix <- sceptre_object_use@covariate_matrix[rownames(sceptre_object_use@covariate_matrix) %in% colnames(full_response_matrix_sim_sparse),]
      
      sceptre_object_use@grna_assignments$grna_group_idxs[pert] <- list(seq(sum(colData(pert_object)$pert == 1)))
    }

    
    message("Running discovery analysis")
    sceptre_object_use <- run_discovery_analysis(
      sceptre_object = sceptre_object_use,
      parallel = FALSE
    )
    
    message("Returning discovery results")
    discovery_result <- get_result(
      sceptre_object = sceptre_object_use,
      analysis = "run_discovery_analysis"
    )
    
    # Add the number of perturbed cells and the rep to each pair
    n_pert_cells <- length(sceptre_object_use@grna_assignments$grna_group_idxs[[pert]])
    discovery_result$num_pert_cells <- n_pert_cells
    discovery_result$rep <- rep
    
    # Save the results
    discovery_results <- data.frame(rbind(discovery_results, discovery_result))
  }
}

# Add the average expression
average_expr = data.frame(
  response_id = rownames(sce), 
  average_expression_all_cells = rowMeans(counts(sce))
)
discovery_results <- left_join(discovery_results, average_expr, by = "response_id")


### SAVE OUTPUT ===============================================================


# save simulation output
message("Saving output to file.")
write_tsv(discovery_results, file = snakemake@output$power_analysis_output)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

