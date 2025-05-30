# Processing held-out CRISPR enhancer screens for benchmarking ENCODE-rE2G

## Overview

This pipeline processes validation datasets to test ENCODE_rE2G performance on predicting enhancer-gene pairs through:

1. **Differential Expression Analysis**: Uses SCEPTRE to identify significant enhancer-gene interactions
2. **Power Analysis**: Simulates perturbation effects at various effect sizes
3. **Data Integration**: Incorporates various other datasets and applies filtering
4. **Benchmarking Output**: Produces final validation datasets for enhancer-gene prediction benchmarking

## Repository Structure

```
ENCODE_Sceptre_Analysis/
├── config/
│   └── config.yml                 # Pipeline configuration parameters
├── workflow/
│   ├── Snakefile                  # Main workflow definition
│   ├── rules/                     # Snakemake rule definitions
│   │   ├── sceptre_setup.smk
│   │   ├── sceptre_power_analysis.smk
│   │   ├── create_encode_output.smk
│   │   └── combine_val_data_and_format.smk
│   ├── scripts/                   # Analysis scripts
│   │   ├── sceptre_setup/
│   │   ├── sceptre_power_analysis/
│   │   ├── encode_datasets/
│   │   └── combine_val_data_and_format/
│   └── envs/                      # Conda environment definitions
├── resources/                     # Input data and reference files
├── results/                       # Pipeline outputs
├── Perturb_Seq_Test_Set_Preprocessing/  # Pre-processing to create pipeline inputs
└── README.md
```

## Requirements

### Software Dependencies
- **Snakemake**: 7.3.2
- **Conda**: 24.11.3

### Environment Management
Dependencies are managed automatically through conda environments defined in `workflow/envs/`. The pipeline will create and activate the appropriate environments for each step.

Key environments include:
- `sceptre_env.yml`: SCEPTRE differential expression analysis
- `analyze_crispr_screen.yml`: Single-cell analysis tools
- `r_process_crispr_data.yml`: Data processing and formatting

## Input Data

### Required Input Files in `resources/`

1. **Externally processed Enhancer-Gene pairs**
   - DC TAP-seq data: `resources/combine_val_data_and_format/DC_TAP_Seq_data.tsv`
     - This file was taken from the [DC TAP Seq processing pipeline](https://github.com/jamesgalante/DC_TAP_Paper/blob/main/results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv)
   - ENCODE-rE2G Training dataset: `resources/combine_val_data_and_format/EPCrisprBenchmark_ensemble_data_GRCh38.tsv`
   - Other test datasets: `resources/create_encode_output/ENCODE/EPCrisprBenchmark/`
   
2. **Raw data input created in Perturb_Seq_Test_Set_Preprocessing - see Data Processing section**
   - Klann et al. 2021: `resources/sceptre_setup/Klann/`
   - Morris et al. 2023: `resources/sceptre_setup/Morrisv1/` & `resources/sceptre_setup/Morrisv2/`
   - Xie et al. 2019: `resources/sceptre_setup/Xie/`

### Data Preprocessing

(Optional as these results are included on SYNAPSE) Before running the main pipeline, you must first preprocess the raw CRISPR screen data.

See the comprehensive guide in: **[Perturb_Seq_Test_Set_Preprocessing/README.md](Perturb_Seq_Test_Set_Preprocessing/README.md)**

This preprocessing step includes:
- Converting raw sequencing data to count matrices
- Creating guide assignment matrices
- Generating metadata files

## Running the Pipeline

### Quick Start

1. **Clone the repository**:
```bash
git clone <repository-url>
cd ENCODE_Test_Dataset_Analysis
```

2. **Complete data preprocessing** (see `Perturb_Seq_Test_Set_Preprocessing/README.md`)
  - Optionally, download the preprocessed data from [SYNAPSE](INCLUDE LINK TO SYNAPSE)

3. **Run the complete pipeline**:
```bash
# For HPC with SLURM
snakemake --profile slurm_profile all

# For local execution (not recommended due to size of matrices)
snakemake --use-conda --cores 8 all
```

### Configuration Profiles

#### SLURM Profile (Recommended for HPC)

Create a profile (e.g., `~/.config/snakemake/slurm_profile/config.yaml`):

```yaml
jobs: 500
cluster: slurm
use-conda: true
notemp: true
default-resources:
    - runtime="13h"
    - mem="32G"
# Add your specific SLURM configuration:
# - slurm_account=your_account
# - slurm_partition=your_partition
# - slurm_extra="--nice"
```

Then run:
```bash
snakemake --profile slurm_profile all
```

#### Local Execution

For local machines or non-SLURM clusters:
```bash
snakemake --use-conda --cores <number_of_cores> all
```

## Pipeline Stages

### 1. SCEPTRE Setup (`sceptre_setup.smk`)
- Downloads genome annotations
- Creates SCEPTRE input objects for each dataset

### 2. Power Analysis (`sceptre_power_analysis.smk`)
- Runs SCEPTRE differential expression analysis
- Performs power simulations at multiple effect sizes (2%, 3%, 5%, 10%, 15%, 20%, 25%, 50%)
- Estimates statistical power for detecting enhancer-gene interactions

### 3. ENCODE Formatting (`create_encode_output.smk`)
- Formats results for ENCODE standards
- Filters based on distance to TSS and statistical power
- Integrates Other test datasets from `resources/create_encode_output/ENCODE/EPCrisprBenchmark/`

### 4. Final Integration (`combine_val_data_and_format.smk`)
- Integrates DC TAP-seq data
- Removes training set overlaps
- Resolves duplicates between datasets
- Produces final validation dataset

## Outputs
- Final Validation Dataset: `results/combine_val_data_and_format/Final_Validation_Dataset.tsv.gz`
