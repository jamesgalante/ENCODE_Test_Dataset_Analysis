# Script: combine_val_data_and_format.smk

# Rule to remove the training pairs from all K562 Validation Datasets & create unfiltered_k562_dc_tap file for analysis
rule remove_training_pairs:
  input:
    combined_validation_w_train_pairs = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    training_data = "resources/combine_val_data_and_format/EPCrisprBenchmark_ensemble_data_GRCh38.tsv"
  output:
    validation = "results/combine_val_data_and_format/validation_wo_training_pairs_no_dc_tap.tsv.gz"
  log: 
    "results/combine_val_data_and_format/logs/remove_training_pairs.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/remove_training_pairs.R"
    
# Add the results from the DC TAP Seq experiments to the other perturb-seq / flowFISH results
rule add_DC_TAP_Seq_pairs:
  input:  
    validation = "results/combine_val_data_and_format/validation_wo_training_pairs_no_dc_tap.tsv.gz",
    dc_tap_data = "resources/combine_val_data_and_format/DC_TAP_Seq_data.tsv"
  output:
    full_validation_dataset = "results/combine_val_data_and_format/full_validation.tsv.gz"
  log: 
    "results/combine_val_data_and_format/logs/add_DC_TAP_Seq_pairs.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/add_DC_TAP_Seq_pairs.R"  

# Rule to resize and merge CRISPR elements in DC-TAP-seq datasets
rule resize_crispr_elements:
  input:
    crispr_data = "results/combine_val_data_and_format/full_validation.tsv.gz"
  output:
    combined_output = "results/combine_val_data_and_format/full_validation_resized_elements.tsv.gz"
  log: "results/combine_val_data_and_format/logs/resize_crispr_elements.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/combine_val_data_and_format/resize_crispr_elements.R"
