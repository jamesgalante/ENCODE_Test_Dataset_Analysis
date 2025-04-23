# ENCODE_Validation_Dataset_Analysis

This repository collects and processes validation datasets to test ENCODE_rE2G performance on predicting held-out enhancer-gene pairs. It first runs differential expression using SCEPTRE. A power analysis is run on the negative pairs, selecting for pairs that pass a power threshold of 80% with a 15% simulated decrease in expression. The results are formatted, passed through filters for ENCODE, and formatted for EPBenchmarking. Finally, two datasets from the DC TAP experiments are added to the resulting validation dataset.

## Inputs
- DC TAP Seq results: 
  - resources/combine_val_data_and_format/DC_TAP_Seq_data.tsv
- Training Dataset results: 
  - resources/combine_val_data_and_format/EPCrisprBenchmark_ensemble_data_GRCh38.tsv
- All the validation dataset raw counts matrices, guide matrices, metadata and guide targets files:
  - resources/sceptre_setup/*
- The flowFISH ENCODE validation datasets processed by Andreas:
  - resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_HCRFlowFISH_perCRE_filt.tsv.gz
  - resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_HCT116_FlowFISH_GRCh38.tsv.gz
  - resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Nasser2021_GM12878_GRCh38.tsv.gz
  - resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Nasser2021_Jurkat_GRCh38.tsv.gz

## All intermediate files are generated dynamically

## Outputs
- Full validation dataset formatted for benchmarking:
  - results/combine_val_data_and_format/full_validation.tsv.gz
- Full validation dataset with DC TAP Datasets resized:
  - results/combine_val_data_and_format/full_validation_resized_elements.tsv.gz

## NOTE
Not all files needed to rerun the pipeline are on github - as some are too large (>50Mb)
- Specific files in process_validation_datasets
  - results/process_validation_datasets/-/power_analysis/combined_power*
  - results/process_validation_datasets/-/perturb_sce.rds
  - results/process_validation_datasets/-/differential_expression/final_sceptre_object.rds
  - results/process_validation_datasets/-/differential_expression/sceptre_diffex_input.rds
- Folders to ignore
  - results/genome_annotation_files/
- Specific files in sceptre_setup
  - resources/sceptre_setup/-/dge*
  - resources/sceptre_setup/-/metatdata*
  - resources/sceptre_setup/-/perturb*


