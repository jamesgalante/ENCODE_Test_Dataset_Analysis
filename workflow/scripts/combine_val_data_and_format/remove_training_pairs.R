# Script: overlap_h3k27me3_and_ctcf_with_enhancers.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/remove_training_pairs.rda"))
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
val_pairs <- read_tsv(snakemake@input$full_validation_dataset_w_dc_tap_seq)
train_pairs <- read_tsv(snakemake@input$training_data)


### CREATE GRANGES OBJECTS ====================================================

# Subset validation pairs for K562 because "train_pairs" are K562 only
val_pairs_k562 <- val_pairs %>% filter(CellType == "K562")

# Create GRanges objects
message("Creating GRanges objects")
val_gr <- GRanges(
  seqnames = val_pairs_k562$chrom,
  ranges = IRanges(start = val_pairs_k562$chromStart, end = val_pairs_k562$chromEnd),
  names = val_pairs_k562$name
)

train_gr <- GRanges(
  seqnames = train_pairs$chrom,
  ranges = IRanges(start = train_pairs$chromStart, end = train_pairs$chromEnd),
  names = train_pairs$name
)


### OVERLAP TRAIN & VAL =======================================================

# Find overlapping ranges
message("Finding Overlaps")
overlaps <- findOverlaps(val_gr, train_gr)

# Extract overlapping ranges
val_gr_overlaps <- val_gr[queryHits(overlaps)]
train_gr_overlaps <- train_gr[subjectHits(overlaps)]

# Extract gene names from the `names` metadata
val_gr_names <- sapply(strsplit(mcols(val_gr_overlaps)$names, "\\|"), `[`, 1)
train_gr_names <- sapply(strsplit(mcols(train_gr_overlaps)$names, "\\|"), `[`, 1)

# Identify matches
matches <- val_gr_names == train_gr_names

# Get names of validation pairs that overlap with training pairs
val_names_in_training_set <- val_gr_overlaps[matches,]$names

# Subset validation pairs to exclude those that overlap with training
filtered_val_pairs_k562 <- val_pairs_k562 %>%
  filter(!name %in% val_names_in_training_set)

# Combine with non-K562 validation pairs
filtered_val_pairs <- val_pairs %>%
  filter(CellType != "K562") %>%  # Keep all non-K562 pairs
  bind_rows(filtered_val_pairs_k562) %>%  # Add filtered K562 pairs back
  distinct()  # Ensure no duplicates


### SUMMARIZE FILTERING =======================================================

# Calculate original stats
message("Summarizing overlaps")
original_total <- nrow(val_pairs)
original_positives <- sum(val_pairs$Regulated, na.rm = TRUE)
original_negatives <- original_total - original_positives

# Calculate filtered stats
filtered_total <- nrow(filtered_val_pairs)
filtered_positives <- sum(filtered_val_pairs$Regulated, na.rm = TRUE)
filtered_negatives <- filtered_total - filtered_positives

# Compute statistics for K562-specific filtering
removed_k562 <- nrow(val_pairs_k562) - nrow(filtered_val_pairs_k562)
removed_k562_positives <- sum(val_pairs_k562$Regulated, na.rm = TRUE) - sum(filtered_val_pairs_k562$Regulated, na.rm = TRUE)
removed_k562_negatives <- removed_k562 - removed_k562_positives

# Step 11: Report summary
message("Original Validation Dataset:")
message("- Total pairs: ", original_total)
message("- Positive pairs: ", original_positives)
message("- Negative pairs: ", original_negatives)

message("\nFiltered Validation Dataset:")
message("- Total pairs: ", filtered_total)
message("- Positive pairs: ", filtered_positives)
message("- Negative pairs: ", filtered_negatives)

message("\nChanges Due to K562 Filtering:")
message("- Removed K562 pairs: ", removed_k562)
message("- Removed K562 positives: ", removed_k562_positives)
message("- Removed K562 negatives: ", removed_k562_negatives)


### SAVE OUTPUT ==============================================================

message("Saving output files")
write_tsv(filtered_val_pairs, snakemake@output$full_validation_dataset_w_dc_tap_seq_wo_training_pairs)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)

