
# Function to load and preprocess the data
load_and_process_data <- function(merged, pred_config) {
  
  pred_config <- processPredConfig(pred_config, merged = merged)
  
  pred_config <- pred_config %>%
    filter(pred_id != "baseline" | (pred_id == "baseline" & pred_col %in% c("distToTSS", "nearestExprGene", "randomExprGene")))
  
  merged <- processMergedData(merged, pred_config = pred_config, filter_valid_connections = TRUE,
                              include_missing_predictions = TRUE,
                              distToTSS_as_kb = TRUE)
  
  merged <- merged[merged$measuredGeneSymbol != "BCYRN1", ]
  merged <- merged[merged$name != "SEC22B|chr1:119920115-119920791:.", ]
  
  list(merged = merged, pred_config = pred_config)
}

# Function to calculate AUPRC and performance summaries
calculate_auprc_and_performance <- function(merged, pred_config) {
  pos_col <- "Regulated"
  
  pr <- applyCellTypes(merged, .fun = calcPRCurves, pred_config = pred_config, pos_col = pos_col)
  pr_table <- rbindlist(pr, idcol = "ExperimentCellType")
  
  perf_summary <- applyCellTypes(merged, .fun = makePRSummaryTableBS, pred_config = pred_config,
                                 pos_col = pos_col)
  
  perf_summary <- perf_summary %>% 
    bind_rows(.id = "cell_type") %>% 
    left_join(select(pred_config, pred_uid, pred_name_long), by = "pred_uid") %>% 
    relocate(pred_name_long, .after = pred_uid)
  
  list(pr = pr, pr_table = pr_table, perf_summary = perf_summary)
}

# Function to create PRC plots
create_prc_plots <- function(merged, pr, pred_config) {
  pos_col <- "Regulated"
  
  if (all(is.na(pred_config$color))) {
    pred_colors <- NULL
  } else {
    pred_colors <- deframe(select(pred_config, pred_name_long, color))
  }
  
  n_pos <- applyCellTypes(merged, .fun = calcNPos, pos_col = pos_col)
  pct_pos <- applyCellTypes(merged, .fun = calcPctPos, pos_col = pos_col)
  
  n_pairs <- applyCellTypes(merged, .fun = function(df) n_distinct(df$name))
  pr_title <- paste0("Every_Validation_Dataset_All_Cell_Types", " (", unlist(n_pairs[names(pr)]), " pairs)")
  
  n_pos <- n_pos[names(pr)]
  pct_pos <- pct_pos[names(pr)]
  pr_plots <- mapply(FUN = makePRCurvePlot, pr_df = pr, n_pos = n_pos, pct_pos = pct_pos, plot_name = pr_title,
                     MoreArgs = list(pred_config = pred_config,
                                     min_sensitivity = 0.7,
                                     line_width = 0.8, point_size = 3,
                                     text_size = 15, colors = pred_colors,
                                     plot_thresholds = FALSE),
                     SIMPLIFY = FALSE)
  
  return(pr_plots)
}

# Outer function to run the entire process
auprc_calc <- function(merged, pred_config) {
  # Step 1: Load and process data
  data <- load_and_process_data(merged, pred_config)
  
  # Step 2: Calculate AUPRC and performance summaries
  performance <- calculate_auprc_and_performance(data$merged, data$pred_config)
  
  # Step 3: Create PRC plots
  #prc_plot <- create_prc_plots(data$merged, performance$pr, data$pred_config)
  
  list(data = data, performance = performance)
  #list(pr_table = performance$pr_table, perf_summary = performance$perf_summary, prc_plot = prc_plot)
}

# # Extract the desired resulting AUPRC values from the resulting auprc list
# extract_auprc_results <- function(auprc_results, results_cell_types, dataset_names, condition_type) {
#   auprc_df <- mapply(function(x, cell_types) {
#     
#     # Apply the appropriate filtering based on condition_type
#     if (condition_type == "ENCODE-rE2G") {
#       filtered_df <- x$perf_summary %>%
#         filter(pred_name_long == "ENCODE-rE2G", cell_type == cell_types)
#     } else if (condition_type == "ABC-DNase") {
#       filtered_df <- x$perf_summary %>%
#         filter(pred_uid == "ABCdnase.ABC.Score", cell_type == cell_types)
#     } else {
#       stop("Invalid condition_type provided. Use 'ENCODE-rE2G' or 'ABC-DNase'.")
#     }
#     
#     # Select the relevant columns
#     filtered_df %>%
#       select(AUPRC, AUPRC_lowerCi, AUPRC_upperCi)
#     
#   }, auprc_results, results_cell_types, SIMPLIFY = FALSE)
#   
#   # Combine the results into a single dataframe
#   auprc_df <- do.call(rbind, auprc_df)
#   
#   # Add the dataset name as a new column
#   auprc_df$Name <- dataset_names
#   
#   return(auprc_df)
# }
















