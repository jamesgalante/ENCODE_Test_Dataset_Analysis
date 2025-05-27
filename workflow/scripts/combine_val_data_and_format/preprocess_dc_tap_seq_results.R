# Script: preprocess_dc_tap_seq_results.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/preprocess_dc_tap_seq_results.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading required packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
})

message("Loading input files")
dc_tap_data <- read_tsv(snakemake@input$dc_tap_data)


### FORMAT DC TAP FOR ENCODE ==================================================

# Subset DistalElement_Gene for pairs that have a gencode distance <1Mb -> Create a Validation_DistalElement_Gene column
dc_tap_data_formatted <- dc_tap_data %>% 
  mutate(Validation_DistalElement_Gene = (abs(distance_to_gencode_gene_TSS) < 1e6) & (DistalElement_Gene))

# Format the dc_tap_data to match the validation_data format
dc_tap_data_formatted <- dc_tap_data_formatted %>%
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
  select(chrom, chromStart, chromEnd, name, EffectSize, chrTSS, startTSS, endTSS, measuredGeneSymbol, 
         Significant, pValueAdjusted, PowerAtEffectSize10, PowerAtEffectSize15, PowerAtEffectSize20, 
         PowerAtEffectSize25, PowerAtEffectSize50, ValidConnection, CellType, Reference, Regulated, Dataset)


### RESIZE AND MERGE FUNCTIONS ================================================

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
  
  # Resize elements and merge overlaps
  resize_elements <- width(elements) < size
  elements[resize_elements] <- resize(elements[resize_elements], width = size, fix = "center")
  
  # Get merged element uids for each resized element
  elements_merged <- reduce(elements, with.revmap = TRUE)
  
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


### DATA PROCESSING ===========================================================

message("Resizing and merging elements")
# Merge elements for DC-TAP-seq datasets
k562_dc_tap <- resize_crispr_elements(filter(dc_tap_data_formatted, Dataset == "K562_DC_TAP_Seq"), size = 500)
wtc11_dc_tap <- resize_crispr_elements(filter(dc_tap_data_formatted, Dataset == "WTC11_DC_TAP_Seq"), size = 500)

# Combine with other CRISPR data to create new combined file
formatted_resized_and_merged_dc_tap_output <- bind_rows(k562_dc_tap, wtc11_dc_tap)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(formatted_resized_and_merged_dc_tap_output, snakemake@output$formatted_resized_and_merged_dc_tap_output)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)