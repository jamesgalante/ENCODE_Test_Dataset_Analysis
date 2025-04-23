# Script: add_DC_TAP_Seq_pairs.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/add_DC_TAP_Seq_pairs.rda"))
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
validation_data <- read_tsv(snakemake@input$validation)
dc_tap_data <- read_tsv(snakemake@input$dc_tap_data)

# Format the dc_tap_data to match the validation_data format
dc_tap_data_formatted <- dc_tap_data %>%
  filter(Validation_DistalElement_Gene == TRUE) %>%
  mutate(
    name = paste0(gene_symbol, "|", targeting_chr_hg38, ":", targeting_start_hg38, "-", targeting_end_hg38, ":."),
    EffectSize = pct_change_effect_size / 100,
    ValidConnection = TRUE,
    Reference = paste0(cell_type, "_DC_TAP_Seq"),
    Regulated = ifelse(significant == TRUE & EffectSize < 0, TRUE, FALSE),
    Dataset = Reference
  ) %>%
  dplyr::rename(
    chrom = targeting_chr_hg38,
    chromStart = targeting_start_hg38,
    chromEnd = targeting_end_hg38,
    chrTSS = chrTSS_hg38,
    startTSS = startTSS_hg38,
    endTSS = endTSS_hg38,
    measuredGeneSymbol = gene_symbol,
    Significant = significant,
    pValueAdjusted = sceptre_adj_p_value,
    PowerAtEffectSize10 = power_at_effect_size_10,
    PowerAtEffectSize15 = power_at_effect_size_15,
    PowerAtEffectSize20 = power_at_effect_size_20,
    PowerAtEffectSize25 = power_at_effect_size_25,
    PowerAtEffectSize50 = power_at_effect_size_50,
    CellType = cell_type
  ) %>%
  select(colnames(validation_data))


### SPLIT DC TAP DATA BY CELL TYPE ===========================================

message("Splitting DC TAP data by cell type")
dc_tap_k562 <- dc_tap_data_formatted %>% 
  filter(CellType == "K562")

dc_tap_wtc11 <- dc_tap_data_formatted %>% 
  filter(CellType == "WTC11")

# Filter validation data for K562 cells
validation_k562 <- validation_data %>%
  filter(CellType == "K562")


### FIND OVERLAPS BETWEEN K562 DC TAP AND VALIDATION K562 ====================

message("Finding overlaps between K562 datasets")
# Create GRanges objects for overlap detection - fixed the n() error
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


### PROCESS OVERLAPPING REGIONS ==============================================

message("Processing overlapping regions")
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
combined_data <- bind_rows(
  # Validation data (excluding K562)
  validation_data %>% filter(CellType != "K562"),
  # Filtered K562 validation data (after removing overlaps we don't want)
  validation_k562,
  # DC TAP K562 data (after removing overlaps we don't want)
  dc_tap_k562,
  # DC TAP WTC11 data (all included without overlap checking)
  dc_tap_wtc11
)

# Sort according to genomic coordinates
combined_data <- combined_data %>%
  arrange(chrom, chromStart, chromEnd, measuredGeneSymbol)


### SAVE OUTPUT ===============================================================

message("Saving output files")
write_tsv(combined_data, file = snakemake@output$full_validation_dataset)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)