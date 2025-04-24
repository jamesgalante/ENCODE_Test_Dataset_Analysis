# ENCODE Validation Dataset Analysis

## Overview
This repository processes validation datasets to test ENCODE_rE2G performance on predicting enhancer-gene pairs. The pipeline:
1. Runs differential expression analysis using SCEPTRE
2. Performs power analysis on negative pairs (80% threshold, 15% expression decrease)
3. Formats results and applies ENCODE filters
4. Formats data for EPBenchmarking
5. Incorporates two DC TAP experiment datasets

## COPY DETAILED METHODS FROM GOOGLE DOC

## Inputs
- DC TAP Seq Data: `resources/combine_val_data_and_format/DC_TAP_Seq_data.tsv`
- Training Dataset: `resources/combine_val_data_and_format/EPCrisprBenchmark_ensemble_data_GRCh38.tsv`
- Validation Raw Data: `resources/sceptre_setup/*` (contains count matrices, guide matrices, metadata, guide targets)
- ENCODE Validation Datasets:
  - `resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_HCRFlowFISH_perCRE_filt.tsv.gz`
  - `resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_HCT116_FlowFISH_GRCh38.tsv.gz`
  - `resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Nasser2021_GM12878_GRCh38.tsv.gz`
  - `resources/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Nasser2021_Jurkat_GRCh38.tsv.gz`

## Intermediate Files
All intermediate files are generated dynamically during pipeline execution.

## Outputs
- Full validation dataset formatted for benchmarking:
  `results/combine_val_data_and_format/full_validation.tsv.gz`
- Validation dataset with resized DC TAP elements:
  `results/combine_val_data_and_format/full_validation_resized_elements.tsv.gz`

## Important Notes
Some large files (>50MB) are not included in the GitHub repository. These files are listed here, and can be found in the following directory on oak if needed: `/oak/stanford/groups/engreitz/Users/jgalante/ENCODE_Sceptre_Analysis`:

- Excluded files in `process_validation_datasets/`:
  - `results/process_validation_datasets/-/power_analysis/combined_power*`
  - `results/process_validation_datasets/-/perturb_sce.rds`
  - `results/process_validation_datasets/-/differential_expression/final_sceptre_object.rds`
  - `results/process_validation_datasets/-/differential_expression/sceptre_diffex_input.rds`
- Excluded folder: `results/genome_annotation_files/`
- Excluded files in `sceptre_setup/`:
  - `resources/sceptre_setup/-/dge*`
  - `resources/sceptre_setup/-/metadata*`
  - `resources/sceptre_setup/-/perturb*`
