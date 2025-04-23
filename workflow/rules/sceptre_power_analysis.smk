# Snakemake rules to run power analysis using Sceptre

# Run sceptre differential expression with "union"
rule sceptre_differential_expression:
  input:
    sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds"
  output:
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/sceptre_differential_expression_{sample}.log"
  conda:
    "../envs/sceptre_env.yml"
  resources:
    mem = "32G",
    time = "12:00:00"
  script:
    "../scripts/sceptre_power_analysis/sceptre_differential_expression.R"

rule create_sce:
  input:
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds",
    guide_targets =  "resources/sceptre_setup/{sample}/guide_targets.tsv"
  output:
    perturb_sce = "results/process_validation_datasets/{sample}/perturb_sce.rds"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/create_sce_{sample}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "196G",
    time = "4:00:00"
  script:
     "../scripts/sceptre_power_analysis/create_sce_object.R"
    

# Define a function to get the number of batches for a given sample
def get_n_batches(wildcards):
    return config["process_validation_datasets"]["power_simulations"]["n_batches"].get(wildcards.sample, 50)

# Use a checkpoint instead of a regular rule for splitting
checkpoint split_target_response_pairs:
  input:
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds"
  output:
    # Use directory as output to avoid specifying exact files
    directory("results/process_validation_datasets/{sample}/pair_splits/")
  params:
    batches = lambda wildcards: get_n_batches(wildcards)
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/split_target_response_pairs_{sample}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  shell:
    """
    # Create the output directory
    mkdir -p {output}
    
    # Create a temporary outputs file listing all expected split files
    for i in $(seq 1 {params.batches}); do
      echo "{output}/gene_gRNA_group_pairs_$i.txt" >> {output}/split_files_list.txt
    done
    
    # Run the R script, passing the list of files as output
    Rscript workflow/scripts/sceptre_power_analysis/split_target_response_pairs_checkpoint.R \
      --input={input.gene_gRNA_group_pairs} \
      --output=$(cat {output}/split_files_list.txt | tr '\\n' ' ') \
      --params.batches={params.batches} \
      --log={log}
    
    # Clean up the temporary file
    # rm {output}/split_files_list.txt
    """

# Function to get files created by the checkpoint
def get_split_files(wildcards):
    # Wait for the checkpoint to complete
    checkpoint_output = checkpoints.split_target_response_pairs.get(**wildcards).output[0]
    # Get all files in the directory matching the pattern
    pattern = os.path.join(checkpoint_output, "gene_gRNA_group_pairs_{split}.txt")
    split_nums = glob_wildcards(pattern).split
    return expand(pattern, split=split_nums)

    
# Run the power simulation with sceptre for each split
rule sceptre_power_analysis:
  input:
    gene_gRNA_group_pairs_split = "results/process_validation_datasets/{sample}/pair_splits/gene_gRNA_group_pairs_{split}.txt",
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds",
    gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
    perturb_sce = "results/process_validation_datasets/{sample}/perturb_sce.rds"
  output:
    power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/effect_size_{effect_size}/power_analysis_output_{split}.tsv"
  params:
    reps = config["process_validation_datasets"]["power_simulations"]["n_reps"],
    guide_sd = config["process_validation_datasets"]["power_simulations"]["guide_sd"]
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/sceptre_power_analysis_es{effect_size}_split{split}_{sample}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_power_analysis/sceptre_power_analysis.R"


# Function to get power analysis output files
def get_power_analysis_outputs(wildcards):
    # Wait for the checkpoint to complete
    checkpoint_output = checkpoints.split_target_response_pairs.get(sample=wildcards.sample).output[0]
    # Get all split files created by the checkpoint
    pattern = os.path.join(checkpoint_output, "gene_gRNA_group_pairs_{split}.txt")
    split_nums = glob_wildcards(pattern).split
    return expand("results/process_validation_datasets/{sample}/power_analysis/effect_size_{effect_size}/power_analysis_output_{split}.tsv", sample=wildcards.sample, effect_size=wildcards.effect_size, split=split_nums)

# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis:
  input:
    splits = get_power_analysis_outputs
  output:
    combined_power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/combine_sceptre_power_analysis_es{effect_size}_{sample}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_power_analysis/combine_sceptre_power_analysis.R"


# Compute the power from the power simulations
rule compute_power_from_simulations:
  input:
    combined_power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds"
  output:
    power_analysis_results = "results/process_validation_datasets/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/compute_power_from_simulations_es{effect_size}_{sample}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "24G",
    time = "1:00:00"
  script:
    "../scripts/sceptre_power_analysis/compute_power_from_simulations.R"

# format sceptre output for compatibility with ENCODE pipelines
rule format_sceptre_output:
  input:
    power_analysis_results = expand("results/process_validation_datasets/{{sample}}/power_analysis/power_analysis_results_es_{effect_size}.tsv", effect_size = [0.15, 0.25]),
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds",
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
    distances = "results/process_validation_datasets/{sample}/distances.tsv",
    guide_targets = "resources/sceptre_setup/{sample}/guide_targets.tsv"
  output:
    final_output = "results/process_validation_datasets/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/format_sceptre_output_{sample}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/sceptre_power_analysis/format_sceptre_output.R"
