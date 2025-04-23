# Script: create_sceptre_diffex_input_Morrisv2.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_sceptre_diffex_input_Morrisv2.rda"))
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
  library(sceptre)
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(R.utils)
  source(file.path(snakemake@scriptdir, "gene_target_pairing_functions.R"))
})

message("Loading input files")
# Import main counts data matrix
dge <- fread(snakemake@input$dge) %>% column_to_rownames("VECTOR") %>% as.matrix() %>% as("dgCMatrix")
# Import perturbation status data frame
perturb_status <- fread(snakemake@input$perturb_status) %>% column_to_rownames("VECTOR") %>% as.matrix() %>% as("dgCMatrix")
# Import the annotation file
annot <- import(snakemake@input$annot)
# Import the guide_targets file
guide_targets <- read_tsv(snakemake@input$guide_targets)
guide_targets <- guide_targets %>% filter(target_type != "TSSCtrl")


### CREATE GENE_GRNA_GROUP_PAIRS ==============================================

# Create the file
target_search_results <- find_genes_near_targets(guide_targets, annot, rownames(dge), max_distance=1e6)
gene_gRNA_group_pairs <- target_search_results[[1]]
errors <- target_search_results[[2]]

# Modify for Sceptre Input specifications
gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>%
  select(grna_group, response_id) %>%
  dplyr::rename(grna_target = grna_group)


### CREATE GRNA_GROUPS_TABLE ================================================== 

# Create the grna_groups_table which is just the guide name matched to the target name
gRNA_groups_table <- guide_targets %>%
  select(name, target_name) %>%
  dplyr::rename(grna_id = name, grna_target = target_name)


### CREATE THE METADATA FILE ==================================================

# Create the metadata indicating the batch for each cell barcode
cell_barcode <- colnames(dge)
cell_batches <- as.factor(sub("^([A-Z])_.*", "\\1", cell_barcode))

# Create the metadata dataframe
metadata <- data.frame(row.names = cell_barcode, batch = cell_batches)


### CREATE SCEPTRE OBJECT =====================================================

# Create sceptre_object
sceptre_object <- import_data(
  response_matrix = dge,
  grna_matrix = perturb_status,
  grna_target_data_frame = gRNA_groups_table,
  moi = "high",
  extra_covariates = metadata
)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = gene_gRNA_group_pairs,
  side = "both",
  grna_integration_strategy = "union",
)

print(sceptre_object)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(gene_gRNA_group_pairs, snakemake@output$gene_gRNA_group_pairs)
saveRDS(gRNA_groups_table, snakemake@output$gRNA_groups_table)
saveRDS(metadata, snakemake@output$metadata)
saveRDS(sceptre_object, snakemake@output$sceptre_diffex_input)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)