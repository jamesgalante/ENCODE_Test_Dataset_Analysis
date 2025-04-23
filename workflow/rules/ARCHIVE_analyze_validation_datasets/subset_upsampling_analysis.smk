# This rule file analyzes the processed expt_pred_merged_annot subset file with a feature and upsampling analysis

# Generate the output paths based on the subsets
config_subsets = config["analyze_validation_datasets"]["subset_upsampling_analysis"]["create_enhancer_subsets"]["subsets"]
subset_outputs = [
    f"results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
    for subset in config_subsets
]

    
# Visualize features of positives in the dataset
rule create_enhancer_categories_and_plot:
  input:
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  output:
    enhancer_class_crispri_positives_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_crispri_positives_plot.pdf",
    enhancer_class_crispri_negatives_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_crispri_negatives_plot.pdf",
    enhancer_class_all_pairs_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_all_pairs_plot.pdf",
    enhancer_class_crispri_positives_heatmap_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_crispri_positives_heatmap_plot.pdf",
    enhancer_class_crispri_negatives_heatmap_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_crispri_negatives_heatmap_plot.pdf",
    enhancer_class_all_pairs_heatmap_plot = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/enhancer_class_all_pairs_heatmap_plot.pdf",
    labelled_combined_validation = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    labelled_combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/create_enhancer_categories_and_plot.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/create_enhancer_categories_and_plot.R"

    
# Let's divide the expt_pred_merged_annot files into different subset so that we can compare their performances
rule create_enhancer_subsets:
  input:
    combined_validation =  "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt"
  output:
    subset_outputs
  params:
    subsets = config_subsets
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/create_enhancer_subsets.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/create_enhancer_subsets.R"
    
    
# Upsample / Downsample Training CRISPRi Positives and Negatives to match the Subset Proportions
rule match_training_to_subset:
  input:
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    subset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
  output:
    balanced_training = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/balanced_training_{subset}_rep{replicate}.rds"
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/match_training_to_subset_{subset}_rep{replicate}.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/match_training_to_subset.R"
    

# Randomly sample training dataset to same class balance as validation subsets
rule sample_training_dataset_for_class_balance:
  input: 
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    subset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
  output:
    randomly_downsampled_training = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/randomly_downsampled_training_{subset}_rep{replicate}.rds"
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/sample_training_dataset_for_class_balance_{subset}_rep{replicate}.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/sample_training_dataset_for_class_balance.R"
    

# Only Downsample - don't duplicate any points
rule downsample_training_and_subset:
  input:
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    subset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
  output:
    downsampled_training = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/downsampled_training_{subset}_rep{replicate}.rds",
    downsampled_original = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/downsampled_original_{subset}_rep{replicate}.rds"
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/downsample_training_and_subset_{subset}_rep{replicate}.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/downsample_training_and_subset.R"
    
    
# Randomly upsample the validation dataset to match the training dataset proportions
rule upsample_subset:
  input:
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    subset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
  output:
    upsampled_subset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/upsampled_original_{subset}_rep{replicate}.rds"
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/upsampled_subset_{subset}_rep{replicate}.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/upsample_subset.R"
  

# Separate rule for running AUPRC on original (non-replicated) subsets
rule run_auprc_original:
  input:
    dataset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/original_{subset}.rds"
  output:
    auprc_results = "results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/original_{subset}.rds"
  params:
    pred_config = config["analyze_validation_datasets"]["subset_upsampling_analysis"]["auprc"]["pred_config"]
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/run_auprc_original_{subset}.log"
  conda:
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/run_auprc.R"

# Modified rule for running AUPRC on replicated analyses
rule run_auprc_replicated:
  input:
    dataset = "results/analyze_validation_datasets/subset_upsampling_analysis/subsets/{dataset_type}_{subset}_rep{replicate}.rds"
  output:
    auprc_results = "results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/{dataset_type}_{subset}_rep{replicate}.rds"
  params:
    pred_config = config["analyze_validation_datasets"]["subset_upsampling_analysis"]["auprc"]["pred_config"]
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/run_auprc_{dataset_type}_{subset}_rep{replicate}.log"
  conda:
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/run_auprc.R"


replicates = config["analyze_validation_datasets"]["subset_upsampling_analysis"]["auprc"]["replicates"]
rule plot_auprc_results:
  input:
    # Original subsets (no replication)
    original_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/original_{subset}.rds",  subset=config_subsets),
    # Replicated analyses
    balanced_training_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/balanced_training_{subset}_rep{replicate}.rds",  subset=config_subsets, replicate=range(1, replicates + 1)),
    downsampled_training_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/downsampled_training_{subset}_rep{replicate}.rds",  subset=config_subsets, replicate=range(1, replicates + 1)),
    downsampled_original_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/downsampled_original_{subset}_rep{replicate}.rds",  subset=config_subsets, replicate=range(1, replicates + 1)),
    randomly_downsampled_training_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/randomly_downsampled_training_{subset}_rep{replicate}.rds",  subset=config_subsets, replicate=range(1, replicates + 1)),
    upsampled_original_subsets = expand("results/analyze_validation_datasets/subset_upsampling_analysis/auprc_results/upsampled_original_{subset}_rep{replicate}.rds",  subset=config_subsets, replicate=range(1, replicates + 1))
  output: 
    auprc_upsampling_summary = "results/analyze_validation_datasets/subset_upsampling_analysis/plots/auprc_summary_plot.pdf"
  params:
    subset_names = config["analyze_validation_datasets"]["subset_upsampling_analysis"]["create_enhancer_subsets"]["subset_names"],
    n_replicates = replicates
  log: "results/analyze_validation_datasets/subset_upsampling_analysis/logs/plot_auprc_results.log"
  conda:
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/subset_upsampling_analysis/plot_auprc_results.R"


