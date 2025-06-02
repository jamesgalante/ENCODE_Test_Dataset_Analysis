# Perturb-Seq Test Set Preprocessing

This directory contains the raw data processing for three Perturb-seq datasets used in the held-out dataset. Each Perturb-seq experiment is processed through dedicated R Markdown workflows following the methods from the original published studies. Each converts raw sequencing data into formatted count matrices compatible with analysis with the SCEPTRE package.

NOTE: This is optional as the final results are included on SYNAPSE, detailed in the parent README.md.

## Overview

All datasets follow a standardized processing workflow:
1. Download raw data from GEO/ENCODE repositories
  - Raw data is either automatically downloaded or needs to be manually downloaded following R Markdown instructions, which are also mentioned below.
2. Process the data according to the original study's methods
2. Filter out guides that do not map to hg38 or map unambiguously
3. Generate gene expression, guide expression, guide annotation, and metadata files

## Directory Structure

```
Perturb_Seq_Test_Set_Preprocessing/
├── Morris/
│   ├── STINGv1/
│   └── STINGv2/
├── Xie/
├── Klann/
└── candidate_cre_data/
│   └── creating_sample1_candidate_cres_bed_file/
```

## Resources

Run time for each R Markdown is variable. Downloading the data takes up a majority of the runtime, but should be no longer than 10 minutes per dataset. Due to the variability in dataset size, RAM requirements also differ between each dataset. Morris et al. 2023 and Klann et al. 2021 require 64GB of RAM while Xie et al. 2019 requires 184GB of RAM to process.

## Dependencies

All R Markdowns were run using the following dependencies:

### Software Versions
- **R version:** 4.1.2
- **Python version:** 3.9 (Analysis of Xie et al. 2019 only)

### R Packages

| Package | Version |
|---------|---------|
| data.table | 1.16.0 |
| dplyr | 1.1.4 |
| easylift | 1.0.0 |
| GEOquery | 2.62.2 |
| GenomicRanges | 1.46.1 |
| ggplot2 | 3.5.1 |
| grid | 4.1.2 |
| gridExtra | 2.3 |
| hdf5r | 1.3.10 |
| httr | 1.4.7 |
| jsonlite | 1.8.9 |
| Matrix | 1.6.5 |
| org.Hs.eg.db | 3.14.0 |
| plyr | 1.8.9 |
| readr | 2.1.5 |
| readxl | 1.4.3 |
| rtracklayer | 1.54.0 |
| Seurat | 5.0.1 |
| tidyr | 1.3.1 |
| tidyverse | 2.0.0 |
---

## Morris et al. (2023)

### Shared Resources
- **`supplementary_tables/science.adh7699_table_s3.xlsx`**
    - Supplementary Table 3 from Morris et al. 2023 detailing guide sequence information for both the PILOT screen (STINGv1) and the full screen (STINGv2)

### STINGv1 Dataset
```
Morris/STINGv1/
├── STINGv1_analysis.Rmd    # Processing workflow
└── blat_results.rds          # BLAT alignment results
```
**Data Source**: GEO (automatic download in Rmd)

### STINGv2 Dataset  
```
Morris/STINGv2/
├── STINGv2_analysis.Rmd    # Processing workflow
└── out.psl                   # BLAT alignment results
```
**Data Source**: GEO (manual download required - see lines 92-94 in Rmd or below)
```
# Run this within the STINGv2 directory
wget -O data/GSE171452_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE171452&format=file"
tar -xf data/GSE171452_RAW.tar -C data
find data/ -type f \( -name "*BeeSTINGseq*" -o -name "*STINGseq-v1*" -o -name "GSE171452_RAW.tar" \) -delete
```

---

## Xie et al. (2019)
```
Xie/
├── Xie_analysis.Rmd                     # Processing workflow
├── Global-analysis-K562-enhancers/      # Code from the Xie paper used to process the raw counts data
└── out.psl                              # BLAT alignment results
```
**Data Source**: GEO (manual download required - see lines 35-37 in Rmd or below)
```
# Run this within the Xie directory
wget -O data/GSE129826_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129826&format=file"
tar -xf data/GSE129826_RAW.tar -C data
gunzip data/*.txt.gz
```

---

## Klann et al. (2021)
Two R Markdowns are used for Klann, first run the `count_data_formatting.Rmd` then run the `guide_info_creation.Rmd`.
```
Klann/
├── count_data_formatting.Rmd                           # Processing workflow
├── guide_info_creation.Rmd                             # Guide coordinate processing
├── gr_object_final_df.csv                              # Lifted guide coordinates (hg19 -> hg38)
└── supplementary_table_16_scCERES_library_gRNAs.xlsx   # Guide annotations from supplementary table 16 of Klann paper
```
**Data Source**: ENCODE (manual download required - see lines 31-35 in count_data_formatting.Rmd or below)
```
mkdir data/Gene_Count_Data/
wget -O data/Gene_Count_Data/ENCFF904ZDX_RAW.tar "https://www.encodeproject.org/files/ENCFF904ZDX/@@download/ENCFF904ZDX.tar.gz"
tar -xf data/Gene_Count_Data/ENCFF904ZDX_RAW.tar -C data/Gene_Count_Data
mv data/Gene_Count_Data/data/gersbachlab/tsk10/novaseq_data/klann_6066_191205/k562_scLib_novaseq_aggr/outs/filtered_feature_bc_matrix/* data/Gene_Count_Data/
rm -r data/Gene_Count_Data/data/
```

---

## Candidate CRE Data
```
candidate_cre_data/
└── creating_sample1_candidate_cres_bed_file/
│   ├── creating_candidate_cres_bed_file.sh/
│   ├── submit_dnase_processing.sh/
│   └── dnase_processing_env.yml/

```
Contains processing workflow for defining candidate cis-regulatory elements (cres) based on K562 DNase-Seq. This file is used for overlapping guides from the original studies with elements taken into account by the ENCODE model.
The following steps were taken in `creating_sample1_candidate_cres_bed_file/` to create the `sample1_candidate_cres.bed` used in each Rmd
```
# Build the dnase_processing_env conda environment
conda env create -f dnase_processing_env.yml

# In submit_dnase_process.sh, modify line 13 to point to your conda.sh in order to activate dnase_processing_env
# E.g. line 13: `source /path/to/your/conda/etc/profile.d/conda.sh`

# Submit the job to slurm which will run `creating_candidate_cres_bed_file.sh`
sbatch submit_dnase_processing.sh

# This will produce a file called `sample1_candidate_cres.bed` in the `candidate_cre_data/` directory
```

---

## Processing Workflow

### Steps to Process All Datasets

1. **Create the DNase bed file in candidate_cre_data/**
   - Follow instructions in the `Candidate CRE Data` section above
   - This creates `sample1_candidate_cres.bed` which is used by all dataset processing workflows

2. **Follow R Markdown workflows**:
   - Datasets can be processed in any order
   - Open the main processing `.Rmd` file
   - Follow data download instructions when prompted
   - (Optional as output files are provided) Run guide alignment when prompted
   - Execute code chunks sequentially

### Expected Outputs

Each dataset produces standardized outputs in `results/` directories:
- **Gene expression matrices**: `dge.txt.gz`
- **Guide perturbation matrices**: `perturb_status.txt.gz`
- **Guide target annotations**: `guide_targets.tsv`
- **Cell metadata (STINGv2 and Xie only)**: `metadata.tsv.gz`

---

## Output Validation

After processing, verify that each dataset produces the required files:

```bash
# Check required outputs for each dataset
ls Morris/STINGv1/results/
ls Morris/STINGv2/results/  
ls Xie/results/
ls Klann/results/
```

Expected files:
- `dge.txt.gz` (gene expression matrix)
- `perturb_status.txt.gz` (guide perturbation matrix)
- `guide_targets.tsv` (guide annotations)

---

## Troubleshooting

### Common Issues
1. **Missing raw data**: Follow download instructions in respective R Markdown files
3. **Python dependencies**: Install required packages for Xie dataset processing
4. **Memory issues**: Some datasets require substantial RAM for processing (Xie et al. 2021 contains a very large matrix)

---

## Integration with Main Pipeline

Once preprocessing is complete, run the following bash commands within the parent directory of `git clone` to copy the results files into the resources of the main pipeline. After this is complete, the snakemake pipeline can be run following the instructions in the parent README.md:

```bash
cp -r Perturb_Seq_Test_Set_Preprocessing/Morris/STINGv1/results/ resources/sceptre_setup/Morrisv1/
cp -r Perturb_Seq_Test_Set_Preprocessing/Morris/STINGv2/results/ resources/sceptre_setup/Morrisv2/
cp -r Perturb_Seq_Test_Set_Preprocessing/Xie/results/ resources/sceptre_setup/Xie/
cp -r Perturb_Seq_Test_Set_Preprocessing/Klann/results/ resources/sceptre_setup/Klann/
```

---

## BLAT Setup

All resulting guide alignment files are provided, so this step is optional. However, if you need to run BLAT alignment yourself, follow the instructions below.

### Installing and Running UCSC BLAT Locally

This tutorial covers downloading and running UCSC's BLAT on your local network. It involves setting up a local server, which you can then submit jobs to through a client.

#### Step 1: Download BLAT Executables

Download BLAT executables from the UCSC downloads page. Choose the appropriate architecture for your system:

```bash
# For Mac with ARM64 architecture
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.arm64/ ./

# For Linux x86_64
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./

# For other architectures, visit: http://hgdownload.soe.ucsc.edu/admin/exe/
```

#### Step 2: Set Up Executable Permissions

Navigate to the blat directory and modify executable permissions:

```bash
cd blat
chmod +x gfServer gfClient blat
```

#### Step 3: Download Reference Genome

Download the hg38 reference genome:

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit
```

#### Step 4: Start the BLAT Server

In your blat directory (which should now contain `blat`, `gfServer`, `gfClient`, and `hg38.2bit`):

```bash
./gfServer start 127.0.0.1 1234 -stepSize=5 hg38.2bit
```

#### Step 5: Submit Queries with BLAT Client

After initializing the server, you can queue jobs using `gfClient`.

**Input file format**: The `to_be_BLAT.fa` file is automatically generated in FASTA format during processing:

```
>query1
ATCGGATCGATACG
>query2
TACTATCTACTACT
```

**Run alignment**: Open a new terminal window, navigate to the `blat` directory, and run:

```bash
./gfClient -minScore=20 -minIdentity=0 127.0.0.1 1234 . to_be_BLAT.fa out.psl
```

Where:
- `to_be_BLAT.fa` is the FASTA file with guide sequences to align
- `out.psl` is the output file containing alignment results
- `-minScore=20` sets minimum alignment score threshold
- `-minIdentity=0` sets minimum identity percentage (0 = no minimum)

#### Step 6: Process Output Files

The output `.psl` file contains alignment coordinates and can be processed using the provided R scripts in each dataset's processing workflow.

### Dataset-Specific BLAT Methods

- **Morris STINGv1**: Uses BLAT API integration (detailed in the R Markdown workflow)
- **Morris STINGv2**: Generates `to_be_BLAT.fa` file for local BLAT alignment
- **Xie dataset**: Generates `to_be_BLAT.fa` file for local BLAT alignment
- **Klann dataset**: Uses easyLift for coordinate conversion (detailed in the R Markdown workflow)

### Stopping the Server

When finished, stop the BLAT server:

```bash
./gfServer stop 127.0.0.1 1234
```
---