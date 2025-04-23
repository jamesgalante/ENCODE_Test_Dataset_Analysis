# Script: resize_crispr_elements.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/resize_crispr_elements.rda"))
message("Saved Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING PACKAGES ==========================================================

message("Loading required packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})


### FUNCTIONS =================================================================

message("Defining functions")

# Function to resize and merge elements in crispr data
resize_crispr_elements <- function(crispr, size = 500, filter_valid_connections = TRUE) {
  
  # Filter for valid E-G pairs only if specified
  if (filter_valid_connections == TRUE) {
    crispr <- filter(crispr, ValidConnection == "TRUE")
  }
  
  # Extract unique elements and create GRanges object
  crispr <- mutate(crispr, element_uid = paste0(chrom, ":", chromStart, "-", chromEnd))
  elements <- distinct(select(crispr, chrom, chromStart, chromEnd, element_uid))
  elements <- makeGRangesFromDataFrame(elements, seqnames.field = "chrom",
                                       start.field = "chromStart", end.field = "chromEnd",
                                       keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # NEW: Keep a copy of the original elements for debugging
  elements_before <- elements
  
  # Resize elements and merge overlaps
  resize_elements <- width(elements) < size
  elements[resize_elements] <- resize(elements[resize_elements], width = size, fix = "center")
  
  # NEW: Keep a copy of the resized elements for debugging
  elements_after <- elements
  
  # Get merged element uids for each resized element
  elements_merged <- reduce(elements, with.revmap = TRUE)
  
  # NEW: Debug step - Show detailed merging information
  debug_element_merging(elements_before, elements_after, elements_merged)
  
  merged_uids <- lapply(elements_merged$revmap, FUN = function(x) elements[x]$element_uid)
  elements_merged$merged_uid <- merged_uids
  
  # Convert to table in long format
  elements_merged_df <- elements_merged %>% 
    as.data.frame() %>% 
    select(chrom = seqnames, chromStart = start, chromEnd = end, merged_uid) %>% 
    unnest(cols = merged_uid)
  
  # Add resized element coordinates to crispr results based on original merged element uids
  crispr_merged <- crispr %>% 
    select(-c(chrom, chromStart, chromEnd)) %>% 
    left_join(elements_merged_df, ., by = c("merged_uid" = "element_uid"))
  
  # For every gene - merged element pair, summarize the Regulated column
  crispr_merged <- crispr_merged %>% 
    mutate(pair_uid = paste0(measuredGeneSymbol, "|", chrom, ":", chromStart, "-", chromEnd)) %>%
    group_by(pair_uid) %>% 
    group_split() %>% 
    map_dfr(resolve_merge)
  
  # Reformat for output 
  crispr_merged <- crispr_merged %>% 
    select(-pair_uid, -merged_uid)
  
  return(crispr_merged)
}

# Function to resolve cases where elements of a pair were merged
resolve_merge <- function(pairs) {
  if (nrow(pairs) > 1) {
    pairs %>% 
      arrange(desc(Regulated), desc(pValueAdjusted)) %>% 
      slice_head(n = 1) %>% 
      mutate(name = pair_uid, merged_uid = paste(pairs$merged_uid, collapse = ",")) %>% 
      mutate(EffectSize = NA_real_, pValueAdjusted = NA_real_)
  } else {
    pairs %>% 
      mutate(name = pair_uid)
  }
}

# Add this new function to the FUNCTIONS section
debug_element_merging <- function(elements_before, elements_after, elements_merged) {
  message("DEBUG - Element merging details:")
  
  # For each merged element, show the elements that were merged
  merge_count <- 0
  
  for (i in 1:length(elements_merged)) {
    merged_element <- elements_merged[i]
    original_indices <- merged_element$revmap[[1]]
    
    # Only print debug info if there was actually merging (more than 1 element)
    if (length(original_indices) > 1) {
      merge_count <- merge_count + 1
      
      # Get the original elements that were merged
      original_elements <- elements_before[original_indices]
      original_df <- original_elements %>%
        as.data.frame() %>%
        select(chrom = seqnames, chromStart = start, chromEnd = end, element_uid)
      
      # Get the resized elements before merging
      resized_elements <- elements_after[original_indices]
      resized_df <- resized_elements %>%
        as.data.frame() %>%
        select(chrom = seqnames, chromStart = start, chromEnd = end, element_uid)
      
      # Get the merged element
      merged_df <- merged_element %>%
        as.data.frame() %>%
        select(chrom = seqnames, chromStart = start, chromEnd = end)
      
      message(paste0("Merged group ", merge_count, " (", merged_df$chrom, ":", merged_df$chromStart, "-", merged_df$chromEnd, ")"))
      message("Original elements:")
      print(original_df)
      message("After resizing:")
      print(resized_df)
      message("Final merged element:")
      print(merged_df)
      message("------------------------------")
    }
  }
  
  message(paste0("Total merged groups: ", merge_count))
}


### DATA PROCESSING ===========================================================

message("Loading input CRISPR data")
# Load combined CRISPR data from snakemake input
crispr_file <- snakemake@input$crispr_data
crispr <- read_tsv(crispr_file, show_col_types = FALSE,
                   col_types = c(ValidConnection = col_character()))

message("Extracting datasets")
# Extract DC TAP-seq data for K562 and WTC11
k562_dc_tap <- filter(crispr, Dataset == "K562_DC_TAP_Seq")
wtc11_dc_tap <- filter(crispr, Dataset == "WTC11_DC_TAP_Seq")
other_crispr <- filter(crispr, !Dataset %in% c("K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"))

message("Resizing and merging elements")
# Merge elements for DC-TAP-seq datasets
k562_dc_tap <- resize_crispr_elements(k562_dc_tap, size = 500)
wtc11_dc_tap <- resize_crispr_elements(wtc11_dc_tap, size = 500)

# Combine with other CRISPR data to create new combined file
merged_crispr <- bind_rows(other_crispr, k562_dc_tap, wtc11_dc_tap)


### SAVE OUTPUTS ==============================================================

message("Saving output files")
# Save merged CRISPR data to output files specified by snakemake
write_tsv(merged_crispr, file = snakemake@output$combined_output)
message("Processing complete")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)