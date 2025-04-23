
# This file contains the majority of code which pairs up genes and targets within a certain distance
# Used for `create_sceptre_diffex_input_....R`

# Filter and prepare gene annotations
prepare_gene_coords <- function(annotation_file, gene_ids) {
  # Filter the annot file
  genes <- annotation_file[annotation_file$type == "gene"]
  genes <- genes[!grepl("_PAR_Y", genes$gene_id), ]
  names(genes) <- sub("\\..*", "", genes$gene_id)
  genes <- genes[names(genes) %in% gene_ids]
  
  # calculate TSS coordinates
  gene_tss <- resize(genes, width = 1, fix = "start")
  
  # Create gene coordinates dataframe
  # Note: In genomic coordinates (like in annotation files):
  # - 'start' is ALWAYS the lower number and 'end' is ALWAYS the higher number
  # - This is true for BOTH + and - strand genes
  # - The strand (+/-) indicates transcription direction, not coordinate order
  # For example:
  #   chr1 5000-7000 +  : Transcribed left to right, start=5000, end=7000
  #   chr1 5000-7000 -  : Transcribed right to left, start=5000, end=7000
  gene_coords <- data.frame(
    gene_chr = as.character(seqnames(gene_tss)),
    gene_start = start(gene_tss),
    gene_end = end(gene_tss),
    gene_strand = as.character(strand(gene_tss)),
    gene = names(gene_tss),
    gene_name = gene_tss$gene_name,
    stringsAsFactors = FALSE
  )
  return(gene_coords)
}

# Prepare perturbation coordinates
prepare_pert_coords <- function(guide_targets) {
  unique_targets <- unique(guide_targets[, c("target_chr", "target_start", "target_end", "target_name", "target_type")])
  pert_coords <- data.frame(
    pert_chr = unique_targets$target_chr,
    pert_start = unique_targets$target_start,
    pert_end = unique_targets$target_end,
    perturbation = unique_targets$target_name,
    stringsAsFactors = FALSE
  )
  return(list(pert_coords = pert_coords, unique_targets = unique_targets))
}

# Vectorized version of distance calculation for multiple perturbations
calculate_strand_aware_distances_vectorized <- function(pert_centers, gene) {
  # Distance calculation considers both genomic coordinates and strand:
  # For + strand genes:
  #   - Negative distance = upstream (pert_center < gene_start)
  #   - Positive distance = downstream (pert_center > gene_end)
  # For - strand genes:
  #   - Positive distance = upstream (pert_center < gene_start)
  #   - Negative distance = downstream (pert_center > gene_end)
  if(gene$gene_strand == "-") {
    distances <- ifelse(pert_centers < gene$gene_start,
                        gene$gene_start - pert_centers,  # positive = upstream
                        ifelse(pert_centers > gene$gene_end,
                               gene$gene_end - pert_centers,  # negative = downstream
                               0))
  } else {
    distances <- ifelse(pert_centers < gene$gene_start,
                        pert_centers - gene$gene_start,  # negative = upstream
                        ifelse(pert_centers > gene$gene_end,
                               pert_centers - gene$gene_end,  # positive = downstream
                               0))
  }
  return(distances)
}

# Report missing pairs
# report_missing_pairs <- function(unique_targets, pairs, gene_coords, max_distance) {
#   targets_without_genes <- setdiff(unique_targets$target_name, unique(pairs$grna_group))
#   if(length(targets_without_genes) > 0) {
#     message("Found ", length(targets_without_genes), " targets without any genes within ", max_distance/1e6, "MB:")
#     targets_without_genes_df <- unique_targets[unique_targets$target_name %in% targets_without_genes, c("target_name", "target_type")]
#     message("\nBreakdown by target type:")
#     print(table(targets_without_genes_df$target_type))
#     message("\nFull list of targets without genes:")
#     print(as.data.frame(targets_without_genes_df))
#   }
#   
#   # Fixed this line to use gene_coords$gene instead of names(genes)
#   genes_without_targets <- setdiff(gene_coords$gene, unique(pairs$response_id))
#   if(length(genes_without_targets) > 0) {
#     message("\nFound ", length(genes_without_targets), " genes without any targets within ", max_distance/1e6, "MB:")
#     print(genes_without_targets)
#   }
# }

report_missing_pairs <- function(unique_targets, pairs, gene_coords, max_distance) {
  # Report targets without genes
  targets_without_genes <- setdiff(unique_targets$target_name, unique(pairs$grna_group))
  if(length(targets_without_genes) > 0) {
    message("Found ", length(targets_without_genes), " targets without any genes within ", max_distance/1e6, "MB:")
    targets_without_genes_df <- unique_targets[unique_targets$target_name %in% targets_without_genes, 
                                               c("target_name", "target_type")]
    message("\nBreakdown by target type:")
    print(table(targets_without_genes_df$target_type))
    message("\nFull list of targets without genes:")
    print(as.data.frame(targets_without_genes_df))
  } else {
    targets_without_genes_df <- data.frame(target_name = character(), target_type = character(), stringsAsFactors = FALSE)
  }
  
  # Report genes without targets: now include gene_name along with gene ID
  genes_without_targets <- setdiff(gene_coords$gene, unique(pairs$response_id))
  if(length(genes_without_targets) > 0) {
    message("\nFound ", length(genes_without_targets), " genes without any targets within ", max_distance/1e6, "MB:")
    missing_genes <- gene_coords %>%
      filter(gene %in% genes_without_targets) %>%
      select(gene, gene_name, gene_chr, gene_start, gene_end)
    print(missing_genes)
  } else {
    missing_genes <- data.frame(gene = character(), gene_name = character(), gene_chr = character(), gene_start = numeric(), gene_end = numeric(), stringsAsFactors = FALSE)
  }
  
  return(list(targets_without_genes_df, missing_genes))
}

# Main function that uses all the modular components
find_genes_near_targets <- function(guide_targets, annotation_file, gene_ids, max_distance=2e6) {
  message("Preparing gene and perturbation coordinates")
  gene_coords <- prepare_gene_coords(annotation_file, gene_ids)
  pert_data <- prepare_pert_coords(guide_targets)
  pert_coords <- pert_data$pert_coords
  unique_targets <- pert_data$unique_targets
  
  message("Calculating distances efficiently by chromosome")
  pairs_list <- list()
  
  # Process chromosome by chromosome
  unique_chrs <- unique(pert_coords$pert_chr)
  
  for(chr in unique_chrs) {
    message(paste("Processing chromosome:", chr))
    
    # Get perturbations and genes on this chromosome
    chr_perts <- pert_coords[pert_coords$pert_chr == chr,]
    chr_genes <- gene_coords[gene_coords$gene_chr == chr,]
    
    if(nrow(chr_perts) > 0 && nrow(chr_genes) > 0) {
      # Calculate all perturbation centers at once
      pert_centers <- floor((chr_perts$pert_start + chr_perts$pert_end) / 2)
      
      # For each gene, calculate distances to all perturbations on same chromosome
      for(j in seq_len(nrow(chr_genes))) {
        gene <- chr_genes[j,]
        
        # Calculate TSS based on strand
        gene_tss <- ifelse(gene$gene_strand == "+", 
                           gene$gene_start, 
                           gene$gene_end)
        
        # Use vectorized distance calculation
        distances <- calculate_strand_aware_distances_vectorized(pert_centers, gene)
        
        # Keep only pairs within distance threshold
        keep_perts <- abs(distances) <= max_distance
        if(any(keep_perts)) {
          pairs_list[[length(pairs_list) + 1]] <- data.frame(
            response_id = gene$gene,
            grna_group = chr_perts$perturbation[keep_perts],
            distance = distances[keep_perts],
            strand = gene$gene_strand,
            gene_tss = gene_tss
          )
        }
      }
    }
  }
  
  message("Combining all pairs")
  # Create the pairs output with distances
  pairs <- do.call(rbind, pairs_list)
  pairs <- pairs[!duplicated(pairs[,c("response_id", "grna_group")]),]
  
  # Add target_type information from guide_targets
  message("Adding target type information")
  # Get unique target info
  target_info <- unique(guide_targets[, c("target_name", "target_type")])
  
  # Merge with pairs
  pairs <- merge(pairs, 
                 target_info, 
                 by.x = "grna_group", 
                 by.y = "target_name", 
                 all.x = TRUE)
  
  # Save the filtered pairs with distances
  message("Saving distance information for filtered pairs")
  write_tsv(pairs, snakemake@output$distances)
  
  # Return original output (without distance column)
  pairs_original <- pairs[,c("response_id", "grna_group", "strand")]
  
  message("Reporting all missing pairs")
  errors <- report_missing_pairs(unique_targets, pairs_original, gene_coords, max_distance)
  
  return(list(pairs_original, errors))
}