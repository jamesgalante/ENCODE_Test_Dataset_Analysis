
# Let's divide the expt_pred_merged_annot files into different subset so that we can compare their performances
rule create_enhancer_feature_plots:
  input:
    combined_validation = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt"
  output:
    properties_plot = "results/analyze_validation_datasets/encode_plots/properties_plot.pdf" # Temp - change appropriately
  params:
  log: "results/analyze_validation_datasets/encode_plots/logs/create_encode_feature_plots.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/create_enhancer_feature_plots.R"
    
    
# merge predictions with experimental data
rule fetch_encode_features:
  input:
    predictions = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred"].values(),
    experiment = lambda wildcards: (
      "/oak/stanford/groups/engreitz/Users/jgalante/ENCODE_Sceptre_Analysis/resources/analyze_validation_datasets/encode_plots/EPCrisprBenchmark_training_data_GRCh38.tsv.gz"  
      if wildcards.dataset == "training"  
      else "results/benchmark_validation_datasets/crispr_benchmarking/{dataset}.tsv.gz"
    ),    
    tss_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["tss_universe"],
    gene_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["gene_universe"],
    pred_config = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred_config"],
    cell_type_mapping = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["cell_type_mapping"].values(),
    expressed_genes = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["expressed_genes"]
  output:
    merged = "results/analyze_validation_datasets/encode_plots/fetch_encode_features/{dataset}_expt_pred_merged_annot.txt.gz"
  params:
    pos_col = "Regulated",
    include_col = "include",
    filter_include_col = False
  log: 
    "results/analyze_validation_datasets/encode_plots/logs/fetch_encode_features_{dataset}.log"
  conda: 
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/fetch_encode_features.R"
  
    
# Combine the feature files with the constructed enhancer classes
rule combine_encode_features_with_enhancer_classes:
  input:
    enhancer_classes = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_{dataset}_expt_pred_merged_annot.txt",
    features_file = "results/analyze_validation_datasets/encode_plots/fetch_encode_features/{dataset}_expt_pred_merged_annot.txt.gz"
  output:
    combined_file = "results/analyze_validation_datasets/encode_plots/combined_features_and_classes/{dataset}_expt_pred_merged_annot.txt"
  params:
  log: "results/analyze_validation_datasets/encode_plots/logs/combine_encode_features_with_enhancer_classes_{dataset}.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/combine_encode_features_with_enhancer_classes.R"
    
    
# Compare the encode features between the training and validation datasets
rule compare_encode_features:
  input:
    expand("results/analyze_validation_datasets/encode_plots/combined_features_and_classes/{dataset}_expt_pred_merged_annot.txt", dataset = ["validation", "training", "unfiltered_k562_dc_tap"])
  output:
    "results/analyze_validation_datasets/encode_plots/compare_encode_features.html"
  conda: 
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "48G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/compare_encode_features.Rmd"
    

# Compute significance between enhancer class proportions in training versus validation
rule chi_squared_enhancer_classes:
  input: 
    expand("results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_{dataset}_expt_pred_merged_annot.txt", dataset = ["validation", "training", "unfiltered_k562_dc_tap"])
  output:
    chi_squared_plot = "results/analyze_validation_datasets/encode_plots/chi_squared_plot.pdf"
  log: "results/analyze_validation_datasets/encode_plots/logs/chi_squared_enhancer_classes.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/chi_squared_enhancer_classes.R"


rule enhancer_category_proportions:
  input:
    labelled_combined_training = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_training_expt_pred_merged_annot.txt",
    labelled_combined_validation = "results/analyze_validation_datasets/subset_upsampling_analysis/expt_pred_merged_annot/labelled_combined_validation_expt_pred_merged_annot.txt"
  output:
    all_proportions_plot_mini = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/all_proportions_plot_mini.pdf",
    all_proportions_plot_full = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/all_proportions_plot_full.pdf",
    positives_proportions_plot_mini = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/positives_proportions_plot_mini.pdf",
    positives_proportions_plot_full = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/positives_proportions_plot_full.pdf",
    all_proportions_plot_select_datasets = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/all_proportions_plot_select_datasets.pdf",
    positives_proportions_plot_select_datasets = "results/analyze_validation_datasets/encode_plots/enhancer_category_proportions/positives_proportions_plot_select_datasets.pdf"
  log: "results/analyze_validation_datasets/encode_plots/logs/enhancer_category_proportions.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/encode_plots/enhancer_category_proportions.R"
    
    
    
    
    
    
    
    
    
    
    
    
    


