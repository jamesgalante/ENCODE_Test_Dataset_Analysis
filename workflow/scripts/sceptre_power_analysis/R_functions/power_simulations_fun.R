## Functions to perform power simulations using specified effect sizes

library(tidyverse)
library(Matrix)

#' Fit negative binomial distributions using DESeq2
#' 
#' Fit negative binomial distributions to estimate dispersion for every gene in a
#' SingleCellExperiment object.
fit_negbinom_deseq2 <- function(response_matrix,
                                simulated_sce,
                                coldata,
                                size_factors = c("ratio", "poscounts", "iterate", "libsize"),
                                fit_type = c("parametric", "local", "mean"),
                                disp_type = c("dispersion", "dispFit", "dispGeneEst", "dispMAP")) {
  
  # Load arguments
  size_factors <- match.arg(size_factors)
  fit_type <- match.arg(fit_type)
  disp_type <- match.arg(disp_type)
  
  # create DESeq2 object containing count data
  dds <- DESeqDataSetFromMatrix(countData = response_matrix,
                                colData = coldata,
                                design = ~ 1)
  
  # compute size factors (removed lib size option for simplicity)
  dds <- estimateSizeFactors(dds, type = size_factors)
  
  # estimate dispersion
  dds <- estimateDispersions(dds, fitType = fit_type)
  
  # add mean expression, dispersion and cell size factors to rowData and colData of simulated_sce
  rowData(simulated_sce)[, "mean"] <- rowData(dds)[, "baseMean"]
  rowData(simulated_sce)[, "dispersion"] <- rowData(dds)[, disp_type]
  rowData(simulated_sce)[, "disp_outlier_deseq2"] <- rowData(dds)[, "dispOutlier"]
  
  # Sample size factors with replacement
  calculated_size_factors <- sizeFactors(dds)
  num_samples <- ncol(simulated_sce)
  # Sample with replacement from the original calculated_size_factors
  new_calculated_size_factors <- sample(calculated_size_factors, size = num_samples, replace = TRUE)
  colData(simulated_sce)[, "size_factors"] <- new_calculated_size_factors
  
  # store dispersion function in simulated_sce
  metadata(simulated_sce)[["dispersionFunction"]] <- dispersionFunction(dds)
  
  return(simulated_sce)
  
}

## GENERALIZABLE SIMULATION FUNCTIONS ==============================================================

#' Simulate Perturb-seq data based on a SingleCellExperiment object with added mean expression,
#' dispersion and size factors from original real data
#' 
#' @param sce Perturb-seq SingleCellExperiment object
#' @param effect_size_mat Effect size matrix
sim_tapseq_sce <- function(sce, effect_size_mat) {
  
  # simulate Perturb-seq count data with parameters from SCE object
  sim_counts <- simulate_tapseq_counts(gene_means = rowData(sce)[, "mean"],
                                       gene_dispersions = rowData(sce)[, "dispersion"],
                                       cell_size_factors = colData(sce)[, "size_factors"],
                                       effect_size_mat = effect_size_mat)
  
  # convert to SingleCellExperiment object with colData and rowData from sce
  output <- SingleCellExperiment(assays = list(counts = sim_counts), rowData = rowData(sce),
                                 colData = colData(sce))
  
  return(output)
  
}

#' Simulate Perturb-seq data
#' 
#' Simulate Perturb-seq UMI counts for a given number of perturbed cells with a specified effect size on
#' specified genes.
#' 
#' @return A matrix with simulated Perturb-seq UMI counts
simulate_tapseq_counts <- function(gene_means, gene_dispersions, cell_size_factors, effect_size_mat,
                                   gene_ids = names(gene_means),
                                   cell_ids = names(cell_size_factors)) {
  
  # number of genes and cells
  n_genes <- length(gene_means)
  n_cells <- length(cell_size_factors)
  
  # make mu matrix for simulation
  mu <- matrix(rep(gene_means, n_cells), ncol = n_cells)
  shuffled_size_factors <- sample(cell_size_factors, n_cells, replace = FALSE)
  mu <- sweep(mu, 2, shuffled_size_factors, "*")  # add cell-to-cell variability based on size factors
  
  # inject perturbation effects by element-wise product of mu and effect_size_mat
  mu <- mu * effect_size_mat
  
  # simulate counts
  message("Simulating tapseq counts")
  Matrix(rnbinom(n_cells * n_genes, mu = mu, size = 1 / unlist(gene_dispersions)), ncol = n_cells,
         dimnames = list(gene_ids, cell_ids))
  
}

# guide-guide variability specific functions -------------------------------------------------------

# function to randomly pick 1 expressed guide per cell
sample_guide <- function(pert_status) {
  
  # Get the number of colums and initialize a return vector
  num_cols <- dim(pert_status)[[2]]
  return_vector <- numeric(0)
  
  # Iterate over columns
  for (col_idx in 1:num_cols) {
    start_idx <- pert_status@p[col_idx] + 1
    end_idx <- pert_status@p[col_idx + 1]
    column_non_zeros <- end_idx - (start_idx - 1)
    
    if (column_non_zeros > 0) {
      selected_i <- pert_status@i[start_idx + sample(0:(column_non_zeros-1), 1)]
      return_vector <- c(return_vector, selected_i)
    } else {
      # Add -1 because the vector gets 1 added in convert_pert_mat_to_vector
      return_vector <- c(return_vector, -1)
    }
  }
  
  # Return the sampled vector
  return(return_vector)
}

# create vector with the gRNA perturbation status for each cell
create_guide_pert_status <- function(pert_status, grna_perts, pert_guides) {
  
  # get gRNA perturbations for perturbed cells
  grnas_pert_cells <- grna_perts[pert_guides, pert_status == 1]
  
  # if >1 guides target the given perturbation convert guide perturbation matrix into vector with
  # unique guide-level perturbation status for each cell. If a cell expresses multiple guides, one
  # is randomly selected
  if (!is.null(nrow(grnas_pert_cells))) {
    grnas_pert_cells <- convert_pert_mat_to_vector(grnas_pert_cells)
  }
  
  # get gRNA perturbations for control cells and also convert these into a vector if needed
  grnas_ctrl_cells <- grna_perts[!rownames(grna_perts) %in% pert_guides, pert_status == 0]
  if (!is.null(nrow(grnas_ctrl_cells))) {
    grnas_ctrl_cells <- convert_pert_mat_to_vector(grnas_ctrl_cells)
  }
  
  # adjust control gRNA status so that they 'come after' targeting gRNAs
  ctrl_perts <- grnas_ctrl_cells > 0
  grnas_ctrl_cells[ctrl_perts] <- grnas_ctrl_cells[ctrl_perts] + max(grnas_pert_cells)
  
  # combine perturbed and control perturbation status vectors
  c(grnas_pert_cells, grnas_ctrl_cells)
  
}

# convert a perturbation status matrix to a vector with a unique status for every perturbation. If a
# cell has >1 perturbations, one is chosen randomly
convert_pert_mat_to_vector <- function(pert_mat) {
  
  # randomly pick one perturbation if a cell has > 1
  pert_mat_sampled <- sample_guide(pert_mat)
  
  # create unique perturbation status for each perturbation and transform to vector
  row_of_pert <- pert_mat_sampled + 1
  names(row_of_pert) <- colnames(pert_mat)
  return(row_of_pert)
  
}

# function to create an effect size matrix with added guide-guide variability in effect size
create_effect_size_matrix <- function(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd) {
  
  # randomly draw effect size variation of guides on every gene
  n_pert_guides <- length(pert_guides)
  n_ctrl_guides <- max(grna_pert_status) - n_pert_guides
  guide_effect_sizes_pert <- vapply(gene_effect_sizes, FUN = rnorm, n = n_pert_guides,
                                    sd = guide_sd, FUN.VALUE = numeric(n_pert_guides))
  guide_effect_sizes_ctrl <- vapply(rlang::rep_along(gene_effect_sizes, 1), FUN = rnorm,
                                    n = n_ctrl_guides, sd = guide_sd,
                                    FUN.VALUE = numeric(n_ctrl_guides))
  guide_effect_sizes <- rbind(guide_effect_sizes_pert, guide_effect_sizes_ctrl)
  
  # set negative guide effect sizes to 0
  guide_effect_sizes[guide_effect_sizes < 0] <- 0
  
  # add row with no effect for non-perturbed cells
  guide_effect_sizes <- rbind(1, guide_effect_sizes)
  
  # pick correct effect sizes for every cell based on it's gRNA perturbation status
  es_mat <- t(guide_effect_sizes[grna_pert_status + 1, ])
  colnames(es_mat) <- names(grna_pert_status)
  
  return(es_mat)
  
}

# function to center effect size matrix with guide variability so that the average effect size per
# gene corresponds to a specified effect size
center_effect_size_matrix <- function(effect_size_mat, pert_status, gene_effect_sizes) {
  
  # get mean effect size for every gene for perturbed and control cells
  mean_es_pert <- rowMeans(effect_size_mat[, pert_status == 1, drop = FALSE])
  mean_es_ctrl <- rowMeans(effect_size_mat[, pert_status == 0, drop = FALSE])
  
  # compute required shift to center guide effect sizes on the specified effect sizes
  pert_shift <- gene_effect_sizes - mean_es_pert
  ctrl_shift <- 1 - mean_es_ctrl
  
  # center guide-level effect sizes on specified effect sizes
  effect_size_mat[, pert_status == 1] <- effect_size_mat[, pert_status == 1] + pert_shift
  effect_size_mat[, pert_status == 0] <- effect_size_mat[, pert_status == 0] + ctrl_shift
  
  # set negative guide effect sizes due to shift to 0
  effect_size_mat[effect_size_mat < 0] <- 0
  
  return(effect_size_mat)
  
}



