# Script: overlap_h3k27me3_and_ctcf_with_enhancers.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/overlap_h3k27me3_and_ctcf_with_enhancers.rda"))
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
})

# Load in the datasets
message("Loading training and validation data")
training <- read_tsv(snakemake@input$training_expt_pred_merged_annot)
validation <- read_tsv(snakemake@input$validation_expt_pred_merged_annot)
unfiltered_k562_dc_tap <- read_tsv(snakemake@input$unfiltered_k562_dc_tap_expt_pred_merged_annot)

# Get list of peak files and extract cell types
peak_files <- snakemake@input$peak_bed_file
cell_types <- tolower(unique(sapply(strsplit(basename(peak_files), "\\."), `[`, 1)))

# Separate CTCF and H3K27me3 peak files
ctcf_files <- peak_files[grep("ctcf", peak_files)]
h3k27me3_files <- peak_files[grep("h3k27me3", peak_files)]


### FUNCTIONS ================================================================

# Function to read bed file and create GRanges object
read_bed_to_granges <- function(file) {
  peaks <- read_tsv(file, col_names = c("chrom", "start", "end"))
  gr <- GRanges(
    seqnames = peaks$chrom,
    ranges = IRanges(start = peaks$start + 1, end = peaks$end)
  )
  return(gr)
}

# Function to create GRanges from enhancer data with size adjustment and extension
create_enhancer_granges <- function(df, mark_type) {
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
  
  # Apply mark-specific extension
  extend_bp <- if(mark_type == "h3k27me3") 175 else 0  # 175bp for H3K27me3, 0 for CTCF
  if(extend_bp > 0) {
    starts <- pmax(starts - extend_bp, 0)  # Ensure starts are not negative
    ends <- ends + extend_bp
  }
  
  # Create GRanges object
  GRanges(
    seqnames = df$chrom,
    ranges = IRanges(start = starts, end = ends)
  )
}

# Function to process overlaps for one cell type
process_cell_type_overlaps <- function(df, ctcf_gr, h3k27me3_gr) {
  # Create GRanges for enhancers with appropriate extensions for each mark type
  enhancer_gr_ctcf <- create_enhancer_granges(df, "ctcf")
  enhancer_gr_h3k27me3 <- create_enhancer_granges(df, "h3k27me3")
  
  # Find overlaps
  ctcf_overlaps <- as.logical(countOverlaps(enhancer_gr_ctcf, ctcf_gr))
  h3k27me3_overlaps <- as.logical(countOverlaps(enhancer_gr_h3k27me3, h3k27me3_gr))
  
  # Add overlap columns
  df$`ctcf overlap` <- ctcf_overlaps
  df$`h3k27me3 overlap` <- h3k27me3_overlaps
  
  return(df)
}


### PROCESS DATA =============================================================

message("Processing overlaps for each cell type")
process_dataset <- function(data) {
  # Convert ExperimentCellType to lower case
  data$cell_type_lower <- tolower(data$ExperimentCellType)
  
  # Process each cell type
  result <- data %>%
    group_by(cell_type_lower) %>%
    group_modify(~{
      cell_type <- .y$cell_type_lower
      
      # Get corresponding peak files for this cell type
      ctcf_file <- ctcf_files[grep(cell_type, ctcf_files)]
      h3k27me3_file <- h3k27me3_files[grep(cell_type, h3k27me3_files)]
      
      # Read peak files
      ctcf_gr <- read_bed_to_granges(ctcf_file)
      h3k27me3_gr <- read_bed_to_granges(h3k27me3_file)
      
      # Process overlaps
      process_cell_type_overlaps(.x, ctcf_gr, h3k27me3_gr)
      
    }) %>%
    ungroup() %>%
    select(-cell_type_lower)
  
  return(result)
}

# Process all datasets
message("Processing all datasets")
training_overlap <- process_dataset(training)
validation_overlap <- process_dataset(validation)
unfiltered_k562_dc_tap_overlap <- process_dataset(unfiltered_k562_dc_tap)


### SAVE OUTPUT ==============================================================

message("Saving output files")
write_tsv(training_overlap, snakemake@output$training_overlap)
write_tsv(validation_overlap, snakemake@output$validation_overlap)
write_tsv(unfiltered_k562_dc_tap_overlap, snakemake@output$unfiltered_k562_dc_tap_overlap)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

