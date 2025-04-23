
# download gencode annotations
rule download_gencode_v32lift37_K562_DC_TAP:
  output: "resources/process_validation_datasets/sceptre_setup/genome_annotation_files/gencode.v32lift37.annotation.gtf.gz"
  params:
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

rule download_gencode_v29:
  output: "resources/sceptre_setup/genome_annotation_files/gencode.v29.annotation.gtf.gz"
  params:
    url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"
    
rule create_sceptre_diffex_input_Morrisv1:
  input:
    dge = "resources/sceptre_setup/Morrisv1/dge.txt.gz",
    perturb_status = "resources/sceptre_setup/Morrisv1/perturb_status.txt.gz",
    guide_targets = "resources/sceptre_setup/Morrisv1/guide_targets.tsv",
    annot = "resources/sceptre_setup/genome_annotation_files/gencode.v29.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/Morrisv1/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/Morrisv1/gRNA_groups_table.rds",
    sceptre_diffex_input = "results/process_validation_datasets/Morrisv1/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/Morrisv1/distances.tsv"
  log: "results/process_validation_datasets/sceptre_setup/logs/create_sceptre_diffex_input_Morrisv1.log"
  conda:
    "../envs/sceptre_env.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_setup/create_sceptre_diffex_input_Morrisv1.R"
    
rule create_sceptre_diffex_input_Morrisv2:
  input:
    dge = "resources/sceptre_setup/Morrisv2/dge.txt.gz",
    perturb_status = "resources/sceptre_setup/Morrisv2/perturb_status.txt.gz",
    guide_targets = "resources/sceptre_setup/Morrisv2/guide_targets.tsv",
    annot = "resources/sceptre_setup/genome_annotation_files/gencode.v29.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/Morrisv2/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/Morrisv2/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/Morrisv2/metadata.rds",
    sceptre_diffex_input = "results/process_validation_datasets/Morrisv2/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/Morrisv2/distances.tsv"
  log: "results/process_validation_datasets/sceptre_setup/logs/create_sceptre_diffex_input_Morrisv2.log"
  conda:
    "../envs/sceptre_env.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_setup/create_sceptre_diffex_input_Morrisv2.R"
  
rule create_sceptre_diffex_input_Xie:
  input:
    dge = "resources/sceptre_setup/Xie/dge.txt.gz",
    perturb_status = "resources/sceptre_setup/Xie/perturb_status.txt.gz",
    guide_targets = "resources/sceptre_setup/Xie/guide_targets.tsv",
    annot = "resources/sceptre_setup/genome_annotation_files/gencode.v29.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/Xie/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/Xie/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/Xie/metadata.rds",
    sceptre_diffex_input = "results/process_validation_datasets/Xie/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/Xie/distances.tsv"
  log: "results/process_validation_datasets/sceptre_setup/logs/create_sceptre_diffex_input_Xie.log"
  conda:
    "../envs/sceptre_env.yml"
  resources:
    mem = "108G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_setup/create_sceptre_diffex_input_Xie.R"
  
rule create_sceptre_diffex_input_Klann:
  input:
    dge = "resources/sceptre_setup/Klann/dge.txt.gz",
    perturb_status = "resources/sceptre_setup/Klann/perturb_status.txt.gz",
    guide_targets = "resources/targets.tsv",
    annot = "resources/sceptre_setup/genome_annotation_files/gencode.v29.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/Klann/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/Klann/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/Klann/metadata.rds",
    sceptre_diffex_input = "results/process_validation_datasets/Klann/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/Klann/distances.tsv"
  log: "results/process_validation_datasets/sceptre_setup/logs/create_sceptre_diffex_input_Klann.log"
  conda:
    "../envs/sceptre_env.yml"
  resources:
    mem = "108G",
    time = "2:00:00"
  script:
    "../scripts/sceptre_setup/create_sceptre_diffex_input_Klann.R"
  
  
  
  
  
  
  


