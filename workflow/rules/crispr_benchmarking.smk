# Script: crispr_benchmarking.smk

# Rule to remove the training pairs from all K562 Validation Datasets & create unfiltered_k562_dc_tap file for analysis
rule remove_training_pairs:
  input:
    combined_validation_w_train_pairs = "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    training_data = "resources/benchmark_validation_datasets/crispr_benchmarking/training_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv"
  output:
    validation = "results/benchmark_validation_datasets/crispr_benchmarking/validation_wo_training_pairs.tsv.gz"
  log: 
    "results/benchmark_validation_datasets/crispr_benchmarking/logs/remove_training_pairs.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/crispr_benchmarking/remove_training_pairs.R"
    
# Add the results from the DC TAP Seq experiments to the other perturb-seq / flowFISH results
rule add_DC_TAP_Seq_pairs:
  input:  
    validation = "results/benchmark_validation_datasets/crispr_benchmarking/validation_wo_training_pairs.tsv.gz",
    dc_tap_data = "resources/benchmark_validation_datasets/crispr_benchmarking/Data/DC_TAP_Seq_data.tsv"
  output:
    full_validation_dataset = "results/benchmark_validation_datasets/crispr_benchmarking/validation.tsv.gz"
  log: 
    "results/benchmark_validation_datasets/crispr_benchmarking/logs/add_DC_TAP_Seq_pairs.log"
  conda: 
    "../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/crispr_benchmarking/add_DC_TAP_Seq_pairs.R"  

# Rule to resize and merge CRISPR elements in DC-TAP-seq datasets
rule resize_crispr_elements:
  input:
    crispr_data = "results/benchmark_validation_datasets/crispr_benchmarking/validation.tsv.gz"
  output:
    combined_output = "results/benchmark_validation_datasets/crispr_benchmarking/validation_resized_elements.tsv.gz",
    k562_output = "results/benchmark_validation_datasets/crispr_benchmarking/K562_DC_TAP_Seq_resized_elements.tsv.gz",
    wtc11_output = "results/benchmark_validation_datasets/crispr_benchmarking/WTC11_DC_TAP_Seq_resized_elements.tsv.gz"
  log:
    "results/benchmark_validation_datasets/crispr_benchmarking/logs/resize_crispr_elements.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/crispr_benchmarking/resize_crispr_elements.R"
    
    
    
    
    
    
    
    
    
    

