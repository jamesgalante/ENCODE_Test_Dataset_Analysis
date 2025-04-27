# Script: combine_val_data_and_format.smk

# Preprocess the DC TAP Seq results to prepare for merging with Validation Dataset
rule preprocess_dc_tap_seq_results:
  input:
    dc_tap_data = "resources/combine_val_data_and_format/DC_TAP_Seq_data.tsv"
  output:
    formatted_resized_and_merged_dc_tap_output = "results/combine_val_data_and_format/formatted_resized_and_merged_dc_tap_output.tsv"
  log: 
    "results/combine_val_data_and_format/logs/preprocess_dc_tap_seq_results.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/preprocess_dc_tap_seq_results.R" 
    
# Combine dc tap seq and other validation dataset
rule combining_dc_tap_and_validation_dataset:
  input:
    formatted_resized_and_merged_dc_tap_output = "results/combine_val_data_and_format/formatted_resized_and_merged_dc_tap_output.tsv",
    validation_dataset = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  output:
    full_validation_dataset_w_dc_tap_seq = "results/combine_val_data_and_format/full_validation_dataset_w_dc_tap_seq.tsv"
  log: 
    "results/combine_val_data_and_format/logs/combining_dc_tap_and_validation_dataset.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/combining_dc_tap_and_validation_dataset.R"
    
# Remove pairs from the full validation dataset that overlap with the training dataset
rule remove_training_pairs:
  input:
    full_validation_dataset_w_dc_tap_seq = "results/combine_val_data_and_format/full_validation_dataset_w_dc_tap_seq.tsv",
    training_data = "resources/combine_val_data_and_format/EPCrisprBenchmark_ensemble_data_GRCh38.tsv"
  output:
    full_validation_dataset_w_dc_tap_seq_wo_training_pairs = "results/combine_val_data_and_format/full_validation_dataset_w_dc_tap_seq_wo_training_pairs.tsv.gz"
  log:
    "results/combine_val_data_and_format/logs/remove_training_pairs.log"
  conda:
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/remove_training_pairs.R"

# Remove duplicates between the dc tap k562 dataset and the other k562 validation datasets
# Note: this was done in create_ensemble in create_encode_output.smk, but since we added the dc tap seq datasets after that step, we have to perform this deduplication again
rule remove_k562_duplicates_between_datasets:
  input:
    full_validation_dataset_w_dc_tap_seq_wo_training_pairs = "results/combine_val_data_and_format/full_validation_dataset_w_dc_tap_seq_wo_training_pairs.tsv.gz"
  output:
    Final_Validation_Dataset = "results/combine_val_data_and_format/Final_Validation_Dataset.tsv.gz"
  log:
    "results/combine_val_data_and_format/logs/remove_k562_duplicates_between_datasets.log"
  conda:
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/remove_k562_duplicates_between_datasets.R" 
