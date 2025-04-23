# Script: overlap_dnase_and_h3k27ac_with_enhancers.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/overlap_dnase_and_h3k27ac_with_enhancers.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================
message("Loading packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(data.table)
})

# Load in the datasets
message("Loading training and validation data")
training <- read_tsv(snakemake@input$training_expt_pred_merged_annot)
validation <- read_tsv(snakemake@input$validation_expt_pred_merged_annot)
unfiltered_k562_dc_tap <- read_tsv(snakemake@input$unfiltered_k562_dc_tap_expt_pred_merged_annot)

# Load enhancer lists
message("Loading enhancer lists")
enhancer_files <- snakemake@input$enhancer_list
enhancer_data <- map(enhancer_files, fread)

### FUNCTIONS ================================================================

# Function to calculate the elbow point using maximum curvature method
calculate_elbow <- function(sorted_data, cdf) {
  # Normalize the data
  normalized_data <- (sorted_data - min(sorted_data)) / (max(sorted_data) - min(sorted_data))
  normalized_cdf <- (cdf - min(cdf)) / (max(cdf) - min(cdf))
  
  # Calculate distances from line connecting start to end point (for maximum curvature)
  start_point <- c(normalized_data[1], normalized_cdf[1])
  end_point <- c(normalized_data[length(normalized_data)], normalized_cdf[length(normalized_cdf)])
  
  distances <- sapply(seq_along(normalized_data), function(i) {
    point <- c(normalized_data[i], normalized_cdf[i])
    dist <- abs((end_point[2] - start_point[2]) * point[1] - 
                  (end_point[1] - start_point[1]) * point[2] +
                  end_point[1] * start_point[2] - end_point[2] * start_point[1]) /
      sqrt((end_point[2] - start_point[2])^2 + (end_point[1] - start_point[1])^2)
    return(dist)
  })
  
  # Return the elbow value
  elbow_value <- sorted_data[which.max(distances)]
  return(elbow_value)
}

# Function to create GRanges from perturbation data with size adjustment
create_perturbation_granges <- function(df) {
  # Convert positions to numeric
  starts <- as.numeric(df$chromStart)
  ends <- as.numeric(df$chromEnd)
  
  # Ensure chromStart <= chromEnd
  if (any(starts > ends, na.rm = TRUE)) {
    stop("Error: Found chromStart greater than chromEnd in the data.")
  }
  
  # Adjust perturbations to be at least 300bp in size
  widths <- ends - starts + 1
  small_indices <- which(widths < 300)

  if (length(small_indices) > 0) {
    centers <- (starts[small_indices] + ends[small_indices]) / 2
    new_starts <- floor(centers - 150)
    new_ends <- ceiling(centers + 149)
    starts[small_indices] <- new_starts
    ends[small_indices] <- new_ends
  }
  
  # Create GRanges object
  GRanges(
    seqnames = df$chrom,
    ranges = IRanges(start = starts, end = ends)
  )
}

# Function to create GRanges from enhancer data (no size adjustment needed)
create_enhancer_granges <- function(df) {
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end),
    normalized_dhs = df$normalized_dhs,
    normalized_h3k27ac = df$normalized_h3K27ac
  )
}

# Function to process signal values for one cell type
process_cell_type_signals <- function(perturbation_data, enhancers) {
  # Create GRanges for perturbations with 300bp minimum
  perturbation_gr <- create_perturbation_granges(perturbation_data)
  
  # Create GRanges for enhancers (with signal values)
  enhancer_gr <- create_enhancer_granges(enhancers)
  
  # Find overlaps
  overlaps <- findOverlaps(perturbation_gr, enhancer_gr)
  
  # Create a results dataframe starting with perturbation data
  result_df <- perturbation_data
  
  # Initialize columns for signal values
  result_df$dhs_signal <- NA_real_
  result_df$h3k27ac_signal <- NA_real_
  
  # Convert overlaps to data frame for easier aggregation
  overlap_df <- as.data.frame(overlaps)
  
  # Aggregate overlapping signals for each query
  # Summing values of dhs and h3k27ac for queries with multiple hits
  dhs_aggregated <- tapply(
    enhancers$normalized_dhs[overlap_df$subjectHits],
    overlap_df$queryHits,
    sum,
    na.rm = TRUE
  )
  
  h3k27ac_aggregated <- tapply(
    enhancers$normalized_h3K27ac[overlap_df$subjectHits],
    overlap_df$queryHits,
    sum,
    na.rm = TRUE
  )
  
  # Assign aggregated values back to result_df based on queryHits
  result_df$dhs_signal[as.numeric(names(dhs_aggregated))] <- dhs_aggregated
  result_df$h3k27ac_signal[as.numeric(names(h3k27ac_aggregated))] <- h3k27ac_aggregated
  
  # Calculate non-promoter medians from enhancer data
  non_promoter_enhancers <- enhancers %>%
    filter(!isPromoterElement)
  
  # Calculate elbow thresholds for DHS and H3K27ac
  dhs_sorted <- sort(na.omit(non_promoter_enhancers$normalized_dhs))
  h3k27ac_sorted <- sort(na.omit(non_promoter_enhancers$normalized_h3K27ac))
  dhs_cdf <- ecdf(dhs_sorted)(dhs_sorted)
  h3k27ac_cdf <- ecdf(h3k27ac_sorted)(h3k27ac_sorted)
  
  elbow_dhs <- calculate_elbow(dhs_sorted, dhs_cdf)
  elbow_h3k27ac <- calculate_elbow(h3k27ac_sorted, h3k27ac_cdf)
  
  # Classify as high or low using elbow thresholds
  result_df$high_low_nonpromoter_dhs <- ifelse(result_df$dhs_signal > elbow_dhs, "high", "low")
  result_df$high_low_nonpromoter_h3k27ac <- ifelse(result_df$h3k27ac_signal > elbow_h3k27ac, "high", "low")
  
  # Classify as high or low based on 0.75 threshold on CDF
  dhs_threshold_75 <- quantile(non_promoter_enhancers$normalized_dhs, 0.75, na.rm = TRUE)
  h3k27ac_threshold_75 <- quantile(non_promoter_enhancers$normalized_h3K27ac, 0.75, na.rm = TRUE)
  
  result_df$cdf75_high_low_nonpromoter_dhs <- ifelse(result_df$dhs_signal > dhs_threshold_75, "high", "low")
  result_df$cdf75_high_low_nonpromoter_h3k27ac <- ifelse(result_df$h3k27ac_signal > h3k27ac_threshold_75, "high", "low")
  
  # Calculate high v low dhs/h3k27ac sites based on the median value without promoters
  result_df$median_high_low_nonpromoter_dhs <- ifelse(result_df$dhs_signal > median(non_promoter_enhancers$normalized_dhs, na.rm = TRUE), "high", "low")
  result_df$median_high_low_nonpromoter_h3k27ac <- ifelse(result_df$h3k27ac_signal > median(non_promoter_enhancers$normalized_h3K27ac, na.rm = TRUE), "high", "low")
  
  return(result_df)
}

# Function to process dataset
process_dataset <- function(data) {
  # Convert ExperimentCellType to lower case
  data$cell_type_lower <- tolower(data$ExperimentCellType)
  
  # Process each cell type
  result <- data %>%
    group_by(cell_type_lower) %>%
    group_modify(~{
      cell_type <- .y$cell_type_lower
      
      # Find corresponding enhancer file for this cell type
      enhancer_file_idx <- grep(cell_type, enhancer_files, value = FALSE)
      
      if (length(enhancer_file_idx) == 0) {
        warning(sprintf("No enhancer data found for cell type: %s", cell_type))
        return(.x)
      }
      
      # Get enhancer data
      cell_enhancers <- enhancer_data[[enhancer_file_idx]]
      
      # Process signals
      process_cell_type_signals(.x, cell_enhancers)
      
    }) %>%
    ungroup() %>%
    select(-cell_type_lower)
  
  return(result)
}


### PROCESS DATA =============================================================

message("Processing signal values for each dataset")
# Process all datasets
training_with_signals <- process_dataset(training)
validation_with_signals <- process_dataset(validation)
unfiltered_k562_dc_tap_with_signals <- process_dataset(unfiltered_k562_dc_tap)


### SAVE OUTPUT ==============================================================

message("Saving output files")
write_tsv(training_with_signals, snakemake@output$training_overlap)
write_tsv(validation_with_signals, snakemake@output$validation_overlap)
write_tsv(unfiltered_k562_dc_tap_with_signals, snakemake@output$unfiltered_k562_dc_tap_overlap)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

