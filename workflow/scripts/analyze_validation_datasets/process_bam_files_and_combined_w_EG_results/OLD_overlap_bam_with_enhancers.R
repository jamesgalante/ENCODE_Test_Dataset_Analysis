# Script to overlap the bam read counts overlapping a peak with an enhancer perturbation

### SETUP =====================================================================

# Get the bam_type from Snakemake wildcards
bam_type <- snakemake@wildcards$bam_type

# Saving image for debugging
save.image(paste0("RDA_objects/overlap_bam_with_enhancers_", bam_type, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Opening log file to collect all messages, warnings, and errors
message("Opening log file")
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

# Load in the expt_pred_merged_annot_data
message("Loading training and validation data")
training <- read_tsv(snakemake@input$training_expt_pred_merged_annot)
validation <- read_tsv(snakemake@input$validation_expt_pred_merged_annot)
unfiltered_k562_dc_tap <- read_tsv(snakemake@input$unfiltered_k562_dc_tap_expt_pred_merged_annot)


### LOADING AND FORMATTING PEAK COUNTS ========================================

# Define the column headers for the peak_counts files
peak_counts_header <- c(
  "Chromosome",                 
  "Start_Position",             
  "End_Position",               
  "Peak_Name",                  
  "Score",                      
  "Strand",                     
  "Signal_Value",               
  "p_value",                    
  "q_value",                    
  "Peak_Point",                 
  "Number_of_Overlapping_Reads",
  "Number_of_Bases_Covered",    
  "Peak_Length",                
  "Fraction_of_Peak_Covered"    
)

# Extract cell types from peak_counts file names
message("Extracting cell types from peak_counts file names")
cell_types <- tolower(sapply(strsplit(basename(snakemake@input$peak_counts_files), "_"), `[`, 1))

# Read in the peak_counts files into a named list and compute CPM
message("Loading and normalizing peak_counts files")
peak_counts_list <- setNames(lapply(snakemake@input$peak_counts_files, function(x) {
  df <- read_tsv(x, col_names = peak_counts_header)
  total_counts <- sum(df$Number_of_Overlapping_Reads)
  df$CPM <- (df$Number_of_Overlapping_Reads / total_counts) * 1e6
  return(df)
}), cell_types)

# Compute median CPM values for each bam type and cell type
message("Computing median CPM values for each bam type and cell type")
median_cpm_list <- list()
for(cell_type in names(peak_counts_list)) {
  peaks_df <- peak_counts_list[[cell_type]]
  median_cpm <- median(peaks_df$CPM, na.rm = TRUE)
  median_cpm_list[[cell_type]] <- median_cpm
}


### PROCESS DATA ==============================================================

# Function to create GRanges with optional extension
create_granges <- function(df, seq_col, start_col, end_col, extend=0) {
  # Ensure that start and end positions are numeric
  starts <- as.numeric(df[[start_col]])
  ends <- as.numeric(df[[end_col]])
  
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
  
  # Apply extension
  starts <- pmax(starts - extend, 0)  # Ensure starts are not negative
  ends <- ends + extend
  
  # Create GRanges object
  gr <- GRanges(
    seqnames = df[[seq_col]],
    ranges = IRanges(start = starts, end = ends)
  )
  return(gr)
}

# Function to determine extension based on bam_type
get_extension <- function(bam_type) {
  if (bam_type == 'h3k27ac') {
    return(175)
  } else if (bam_type == "h3k27me3") {
    return(175)
  } else if (bam_type == 'ctcf') {
    return(0)
  } else if (bam_type == "dnase") {
    return(0)
  } else {
    # Throw an error if bam_type is unsupported
    stop(paste("Error: Unsupported bam_type:", bam_type))
  }
}

process_cell_type <- function(cell_type, df, peaks_df, bam_type, median_cpm) {
  message("Processing cell type: ", cell_type)
  
  # Overlap perturbations with peak file
  extend_bp <- get_extension(bam_type)
  peaks_gr <- create_granges(peaks_df, "Chromosome", "Start_Position", "End_Position")
  df_gr <- create_granges(df, "chrom", "chromStart", "chromEnd", extend=extend_bp)
  overlaps <- findOverlaps(df_gr, peaks_gr)
  
  # Initialize overlap and counts columns
  overlap_col_name <- paste0(bam_type, " overlap")
  counts_col_name <- paste0(bam_type, " overlap counts")
  cpm_col_name <- paste0(bam_type, " CPM")
  high_low_col_name <- paste0(bam_type, " High_Low")
  
  df[[overlap_col_name]] <- FALSE
  df[[counts_col_name]] <- 0
  df[[cpm_col_name]] <- 0  # To store CPM values
  
  if (length(overlaps) > 0) {
    # For each query (enhancer), sum the counts of overlapping peaks
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Raw counts and CPM (normalized counts per million)
    raw_counts_per_enhancer <- tapply(peaks_df$Number_of_Overlapping_Reads[subject_hits], query_hits, sum)
    cpms_per_enhancer <- tapply(peaks_df$CPM[subject_hits], query_hits, sum) 
    
    # Update overlap column to TRUE for overlapping enhancers
    df[[overlap_col_name]][unique(query_hits)] <- TRUE
    
    # Assign summed raw counts and CPM to the corresponding enhancers
    df[[counts_col_name]][as.numeric(names(raw_counts_per_enhancer))] <- raw_counts_per_enhancer
    df[[cpm_col_name]][as.numeric(names(cpms_per_enhancer))] <- cpms_per_enhancer
  }
  # Calculate DistToTSS and classify bam counts as High or Low based on median CPM
  df$DistToTSS <- abs(((df$chromStart + df$chromEnd) / 2) - ((df$startTSS + df$endTSS) / 2))
  df[[high_low_col_name]] <- ifelse(df[[cpm_col_name]] >= median_cpm, "High", "Low")
  
  return(df)
}

# Function to perform overlap and add overlap and counts columns
perform_overlap <- function(data, peak_counts_list, bam_type, median_cpm_list) {
  # Convert ExperimentCellType to lower case for matching
  data$cell_type_lower <- tolower(data$ExperimentCellType)
  
  # Split data by cell type
  data_list <- split(data, data$cell_type_lower)
  
  # Process each cell type using the process_cell_type function
  result_list <- lapply(names(data_list), function(cell_type) {
    df <- data_list[[cell_type]]
    
    if (cell_type %in% names(peak_counts_list)) {
      peaks_df <- peak_counts_list[[cell_type]]
      df <- process_cell_type(cell_type, df, peaks_df, bam_type, median_cpm_list[[cell_type]])
    } else {
      message("No peak counts data for cell type: ", cell_type)
      # Initialize overlap and counts columns with NA
      overlap_col_name <- paste0(bam_type, " overlap")
      counts_col_name <- paste0(bam_type, " overlap counts")
      df[[overlap_col_name]] <- NA
      df[[counts_col_name]] <- NA
    }
    return(df)
  })
  # Combine the list back into a data frame and remove the temporary column
  result <- bind_rows(result_list)
  result <- result %>% select(-cell_type_lower)
  
  return(result)
}

# Process the datasets
message("Processing datasets")
training_overlap <- perform_overlap(training, peak_counts_list, bam_type, median_cpm_list)
validation_overlap <- perform_overlap(validation, peak_counts_list, bam_type, median_cpm_list)
unfiltered_k562_dc_tap_overlap <- perform_overlap(unfiltered_k562_dc_tap, peak_counts_list, bam_type, median_cpm_list)


### SAVE OUTPUT ===============================================================

# Save the results
message("Saving output to files")
write_tsv(training_overlap, snakemake@output$training_overlap)
write_tsv(validation_overlap, snakemake@output$validation_overlap)
write_tsv(unfiltered_k562_dc_tap_overlap, snakemake@output$unfiltered_k562_dc_tap_overlap)


### CLEAN UP ==================================================================

# Close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)

