# This rule is to visualize the DC TAP dataset per guide effect sizes and compare with duplicate pairs in the training datasets
rule dc_tap_effect_size_meta_plot:
  input:
    combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt",
    per_guide_effect_sizes_unfiltered_k562_dc_tap = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/k562_dc_tap_per_guide_effect_sizes.txt",
    # Note that there are two of these in resources. one for encode and one for sceptre setup - figure out which one is final
    guide_targets = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/guide_targets.tsv"
  output:
    "results/analyze_validation_datasets/dc_tap_plots/dc_tap_meta_plot.pdf"
  conda: 
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "48G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/dc_tap_effect_size_meta_plot.R"
    
# Rule to create the proportion dc tap proportion plots without the "control" genes
rule dc_tap_proportions_plots:
  input:
    labelled_combined_validation = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    labelled_combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt",
    biased_genes_to_remove = "resources/analyze_validation_datasets/dc_tap_plots/train_list_w_gene_info.txt"
  output:
    proportions_plot_all = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/proportions_plot_all.pdf",
    proportions_plot_pos = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/proportions_plot_pos.pdf",
    fold_change_plot = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/fold_change_plot.pdf",
    chi_squared_results = "results/analyze_validation_datasets/dc_tap_plots/proportions_plots/chi_squared_results.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/dc_tap_proportions_plots.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/dc_tap_proportions_plots.R"
  

  
  

