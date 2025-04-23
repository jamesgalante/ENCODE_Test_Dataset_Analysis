# Script: split_target_response_pairs_checkpoint.R

### SETUP =====================================================================

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to extract values from arguments
get_arg_value <- function(args, flag) {
  i <- which(grepl(paste0("^--", flag, "="), args))
  if (length(i) > 0) {
    # Extract value after the equals sign
    value <- sub(paste0("^--", flag, "="), "", args[i])
    return(value)
  }
  return(NULL)
}

# Get input, output, and parameters
input_file <- get_arg_value(args, "input")
batches <- as.numeric(get_arg_value(args, "params.batches"))
log_filename <- get_arg_value(args, "log")

# Extract the first output file (the one prefixed with --output=)
first_output <- get_arg_value(args, "output")

# Find all positional arguments (not prefixed with --) - these are additional output files
positional_args <- args[!grepl("^--", args)]

# Combine all output files
output_files <- c(first_output, positional_args)
message("Found ", length(output_files), " output files in total")

# Extract sample name from input path
sample <- basename(dirname(input_file))

# Create directory for saving debugging info if needed
if (!file.exists("RDA_objects/split_target_response_pairs")) { 
  dir.create("RDA_objects/split_target_response_pairs", recursive = TRUE) 
}
save.image(paste0("RDA_objects/split_target_response_pairs/", sample, ".rda"))
message("Saved Image")

# Open log file to collect messages, warnings, and errors
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading inputs")
gene_gRNA_group_pairs <- readRDS(input_file)

# Get the output directory from the first file
output_dir <- dirname(output_files[1])
message("Output directory: ", output_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


### SPLIT PAIRS FILE =======================================================

# Calculate the target number of unique grna_group for each split
message("Precomputations for splitting pairs")
total_unique_groups <- n_distinct(gene_gRNA_group_pairs$grna_target)
target_per_split <- ceiling(total_unique_groups / batches)

# Make sure that splitting the unique groups will work given the number of batches
if (total_unique_groups < batches) {
  stop(paste0("The number of batches (", batches, ") exceeds the number of unique grna_targets."))
}

# Initialize splits with empty data frames
splits <- vector("list", batches)
names(splits) <- paste0("split", seq_len(batches))
for (i in seq_len(batches)) {
  splits[[i]] <- data.frame(grna_target = character(0), response_id = character(0))
}

# Distribute grna_target to splits trying to even out the number of unique values
message("Distributing gene-gRNA group pairs")
unique_targets <- unique(gene_gRNA_group_pairs$grna_target)
for (i in seq_along(unique_targets)) {
  split_counts <- sapply(splits, function(x) n_distinct(x$grna_target))
  split_with_least <- which.min(split_counts)
  
  # Get the rows for the current grna_target
  current_rows <- gene_gRNA_group_pairs %>% 
    filter(grna_target == unique_targets[i])
  
  # Add the current rows to the appropriate split
  splits[[split_with_least]] <- bind_rows(splits[[split_with_least]], current_rows)
}

### SAVE OUTPUT ===============================================================

# Write each split to the corresponding output file
message("Saving splits to files")
for (i in seq_along(splits)) {
  write_tsv(splits[[i]], file=output_files[i], col_names = FALSE)
}


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)