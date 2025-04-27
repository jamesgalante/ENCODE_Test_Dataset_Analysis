# Script: remove_k562_duplicates_between_datasets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/remove_k562_duplicates_between_datasets.rda"))
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
  library(GenomicRanges)
})

message("Loading input files")
full_validation_dataset_w_dc_tap_seq_wo_training_pairs <- read_tsv(snakemake@input$full_validation_dataset_w_dc_tap_seq_wo_training_pairs)


### SPLIT DC TAP DATA BY CELL TYPE ===========================================

message("Splitting DC TAP data by cell type")
# Extract DC-TAP K562 data
dc_tap_k562 <- full_validation_dataset_w_dc_tap_seq_wo_training_pairs %>% 
  filter(Dataset == "K562_DC_TAP_Seq")

# Extract non-DC-TAP K562 validation data
validation_k562 <- full_validation_dataset_w_dc_tap_seq_wo_training_pairs %>%
  filter(CellType == "K562" & Dataset != "K562_DC_TAP_Seq")


### FIND OVERLAPS BETWEEN K562 DC TAP AND VALIDATION K562 ====================

message("Finding overlaps between K562 datasets")
# Create GRanges objects for overlap detection
dc_tap_k562_gr <- with(dc_tap_k562,
                       GRanges(seqnames = chrom,
                               ranges = IRanges(chromStart, chromEnd),
                               gene = measuredGeneSymbol,
                               idx = 1:nrow(dc_tap_k562)))

validation_k562_gr <- with(validation_k562,
                           GRanges(seqnames = chrom,
                                   ranges = IRanges(chromStart, chromEnd),
                                   gene = measuredGeneSymbol,
                                   idx = 1:nrow(validation_k562)))

# Find overlaps between datasets
overlaps <- findOverlaps(dc_tap_k562_gr, validation_k562_gr)

# Only keep overlaps where gene symbols match
query_hits <- queryHits(overlaps)
subject_hits <- subjectHits(overlaps)
gene_matches <- dc_tap_k562_gr$gene[query_hits] == validation_k562_gr$gene[subject_hits]

overlaps_filtered <- overlaps[gene_matches]
query_hits <- queryHits(overlaps_filtered)
subject_hits <- subjectHits(overlaps_filtered)

# Create mapping of overlapping regions with unique groups
overlap_mapping <- data.frame(
  dc_tap_idx = dc_tap_k562_gr$idx[query_hits],
  validation_idx = validation_k562_gr$idx[subject_hits]
)

# Group overlapping regions to handle cases where multiple regions overlap
overlap_mapping <- overlap_mapping %>%
  group_by(dc_tap_idx) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup()

# Get unique groups
unique_groups <- unique(overlap_mapping$group_id)

# Function to process overlap groups
process_overlap_group <- function(group_id) {
  # Get all pairs in this group
  pairs <- overlap_mapping %>%
    filter(group_id == !!group_id)
  
  # Get corresponding rows
  dc_tap_entries <- dc_tap_k562[pairs$dc_tap_idx, ]
  validation_entries <- validation_k562[pairs$validation_idx, ]
  
  # Check significance
  dc_tap_significant <- any(dc_tap_entries$Significant)
  validation_significant <- any(validation_entries$Significant)
  
  # If either source has significant entries
  if (dc_tap_significant || validation_significant) {
    # If both are significant, compare power values
    if (dc_tap_significant && validation_significant) {
      dc_max_power <- max(dc_tap_entries$PowerAtEffectSize15, na.rm = TRUE)
      val_max_power <- max(validation_entries$PowerAtEffectSize15, na.rm = TRUE)
      
      # Keep the source with higher power
      if (!is.na(dc_max_power) && (is.na(val_max_power) || dc_max_power > val_max_power)) {
        return(list(keep_source = "dc_tap", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
      } else {
        return(list(keep_source = "validation", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
      }
    }
    # If only one source is significant, keep that one
    else if (dc_tap_significant) {
      return(list(keep_source = "dc_tap", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
    } else {
      return(list(keep_source = "validation", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
    }
  }
  # If neither is significant, compare power values
  else {
    dc_max_power <- max(dc_tap_entries$PowerAtEffectSize15, na.rm = TRUE)
    val_max_power <- max(validation_entries$PowerAtEffectSize15, na.rm = TRUE)
    
    # Keep the source with higher power
    if (!is.na(dc_max_power) && (is.na(val_max_power) || dc_max_power > val_max_power)) {
      return(list(keep_source = "dc_tap", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
    } else {
      return(list(keep_source = "validation", dc_tap_indices = pairs$dc_tap_idx, validation_indices = pairs$validation_idx))
    }
  }
}

# Process each overlap group
overlap_decisions <- lapply(unique_groups, process_overlap_group)

# Extract indices to keep or remove
dc_tap_indices_to_remove <- unique(unlist(lapply(overlap_decisions, function(x) {
  if (x$keep_source == "validation") x$dc_tap_indices else NULL
})))

validation_indices_to_remove <- unique(unlist(lapply(overlap_decisions, function(x) {
  if (x$keep_source == "dc_tap") x$validation_indices else NULL
})))

# Remove overlapping rows based on decisions
if (length(dc_tap_indices_to_remove) > 0) {
  dc_tap_k562 <- dc_tap_k562[-dc_tap_indices_to_remove, ]
}

if (length(validation_indices_to_remove) > 0) {
  validation_k562 <- validation_k562[-validation_indices_to_remove, ]
}


### COMBINE ALL DATASETS =====================================================

message("Combining all datasets")
# Combine all datasets
Final_Validation_Dataset <- bind_rows(
  # Non-K562 validation data
  full_validation_dataset_w_dc_tap_seq_wo_training_pairs %>% filter(CellType != "K562"),
  # Filtered K562 validation data (after removing overlaps we don't want)
  validation_k562,
  # DC TAP K562 data (after removing overlaps we don't want)
  dc_tap_k562
)

# Sort according to genomic coordinates
Final_Validation_Dataset <- Final_Validation_Dataset %>%
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(Final_Validation_Dataset, snakemake@output$Final_Validation_Dataset)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)