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
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/benchmark_validation_datasets/crispr_benchmarking/remove_training_pairs.R"
    
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
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/benchmark_validation_datasets/crispr_benchmarking/add_DC_TAP_Seq_pairs.R"  

# merge predictions with experimental data
# rule mergePredictionsWithExperiment:
#   input:
#     predictions = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred"].values(),
#     experiment = "results/benchmark_validation_datasets/crispr_benchmarking/{dataset}.tsv.gz",
#     tss_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["tss_universe"],
#     gene_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["gene_universe"],
#     pred_config = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred_config"],
#     cell_type_mapping = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["cell_type_mapping"].values(),
#     expressed_genes = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["expressed_genes"]
#   output:
#     merged = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/{dataset}_expt_pred_merged_annot.txt.gz"
#   params:
#     pos_col = "Regulated",
#     include_col = "include",
#     filter_include_col = False
#   log: 
#     "results/benchmark_validation_datasets/crispr_benchmarking/logs/mergePredictionsWithExperiment_{dataset}.log"
#   conda: 
#     "../../envs/r_crispr_comparison.yml"
#   resources:
#     mem_mb = 72000
#   script:
#     "../../scripts/benchmark_validation_datasets/crispr_benchmarking/mergePredictionsWithExperiment.R"


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
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/benchmark_validation_datasets/crispr_benchmarking/resize_crispr_elements.R"
    
    
    
    
    
    
    
    
    
    

