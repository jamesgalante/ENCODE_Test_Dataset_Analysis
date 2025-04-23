# This rule file processes the expt_pred_merged_annot outputs for the subset_upsampling_analysis

# Rule order for processing jurkat bed file before downloading others
ruleorder: peak_calling > filter_top_peaks > download_bed_files

# Download and index hg38 for processing Jurkat specific files with BWA
rule download_and_index_hg38:
  output:
    # Store the uncompressed reference genome and the BWA index files
    genome = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa",
    bwt = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa.bwt",
    sa = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa.sa"
  params:
    rsync_path = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["download_and_index_hg38"]["hg38"]
  resources:
    mem = "32G",
    time = "4:00:00"
  shell:
    """
    # Load necessary modules
    ml load system
    ml load biology
    ml load bwa
    
    # Ensure the output directory exists
    mkdir -p resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/
    
    # Download the reference genome using rsync
    rsync -avz --progress {params.rsync_path} resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa.gz
    
    # Uncompress the downloaded genome
    gunzip -f resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa.gz
    
    # Index the genome with BWA
    bwa index resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa
    """

# Process CTCF ChIP seq to be used in the analysis
rule process_jurkat_ctcf_chip:
  input:
    index = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa"
  output:
    bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files/ctcf.bam"
  params:
    srr_ids = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["process_jurkat_ctcf_chip"]["srr_ids"]
  resources:
    mem = "72G",
    time = "6:00:00"
  threads: 8
  shell:
    """
    ml system
    ml biology
    ml samtools
    ml sra-tools/3.0.7
    ml bwa
    
    mkdir -p results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/
    
    # Loop through the SRR IDs, download the FASTQ files
    for srr_id in {params.srr_ids}; do
        # Step 1: Download paired-end FASTQ files from SRA with original headers
        fasterq-dump $srr_id -O results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/ --split-files -F --gzip
    done
    
    for srr_id in {params.srr_ids}; do
      # Align each single-end FASTQ file separately
      bwa mem -t {threads} {input.index} results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.fastq > \
      results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.sam
      
      # Convert each SAM file to BAM
      samtools view -bS results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.sam > \
      results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.bam
    done
    
    # Merge the BAM files into one BAM file
    samtools merge -f -@ {threads} {output.bam} $(for srr_id in {params.srr_ids}; do echo results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.bam; done) 
    """

# This rule is not required by the new pipeline, but is left in for reference on the h3k27ac jurkat data processing pre-ABC model
# # Process H3K27ac ChIP seq to be used in the analysis
# rule process_jurkat_h3k27ac_chip:
#   input:
#     index = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa"
#   output:
#     bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files/h3k27ac.bam"
#   params:
#     srr_ids = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["process_jurkat_h3k27ac_chip"]["srr_ids"]
#   resources:
#     mem = "72G",
#     time = "6:00:00"
#   threads: 8
#   shell:
#     """
#     ml system
#     ml biology
#     ml samtools
#     ml sra-tools/3.0.7
#     ml bwa
# 
#     mkdir -p results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/
# 
#     # Download FASTQ files for each SRR ID
#     for srr_id in {params.srr_ids}; do
#         fasterq-dump $srr_id -O results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/ --split-files -F --gzip
#     done
# 
#     for srr_id in {params.srr_ids}; do
#       # Align each single-end FASTQ file separately
#       bwa mem -t {threads} {input.index} results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.fastq > \
#       results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.sam
# 
#       # Convert each SAM file to BAM
#       samtools view -bS results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.sam > \
#       results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.bam
#     done
# 
#     # Merge the BAM files into one BAM file
#     samtools merge -f -@ {threads} results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_aligned.bam \
#     $(for srr_id in {params.srr_ids}; do echo results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.bam; done)
# 
#     # https://www.biostars.org/p/318974/ has a discussion on deduplicating ChIP-seq data
#     # https://www.biostars.org/p/390305/ has a discussion mentioning markdup being a sufficient replacement for picard
#     # https://www.htslib.org/doc/samtools-markdup.html documentation on markdup
#     # Run samtools fixmate (required for markdup)
#     samtools fixmate -m results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_aligned.bam \
#         results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_fixmate.bam
# 
#     # Sort the BAM file before removing duplicates
#     samtools sort -@ {threads} -o results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_sorted.bam \
#         results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_fixmate.bam
# 
#     # Remove PCR duplicates using samtools markdup
#     samtools markdup -r results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_sorted.bam {output.bam}
#     """
    
# Process H3K27ac ChIP seq to be used in the analysis
rule process_jurkat_h3k27me3_chip:
  input:
    index = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa"
  output:
    bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files/h3k27me3.bam"
  params:
    srr_id = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["process_jurkat_h3k27me3_chip"]["srr_id"]
  resources:
    mem = "72G",
    time = "6:00:00"
  threads: 8
  shell:
    """
    ml system
    ml biology
    ml samtools
    ml sra-tools/3.0.7
    ml bwa
    
    mkdir -p results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/

    # Download FASTQ files for the indicated SRR ID
    fasterq-dump {params.srr_id} -O results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/ --split-files -F --gzip

    # Align the single-end FASTQ file
    bwa mem -t {threads} {input.index} results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.fastq > \
    results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.sam

    # Convert the SAM file to BAM
    samtools view -bS results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.sam > \
    results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.bam

    # Run samtools fixmate (required for markdup)
    samtools fixmate -m results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.bam \
        results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_fixmate.bam

    # Sort the BAM file before removing duplicates
    samtools sort -@ {threads} -o results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_sorted.bam \
        results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_fixmate.bam

    # Remove PCR duplicates using samtools markdup
    samtools markdup -r results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_sorted.bam {output.bam}
    """

# Process the bam file - sort and filter
rule process_bam_files:
  input:
    bam_file = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files/{bam_type}.bam"
  output:
    bam_file_sorted = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files_sorted/{bam_type}_sorted.bam"
  resources:
    mem = "24G",
    time = "4:00:00"
  shell:
    """
    # Load the necessary modules
    ml load system
    ml load biology
    ml load samtools 
    
    # Filtering to quality q30 and sorting the filtered BAM file
    samtools view -q 30 -b {input.bam_file} | samtools sort -o {output.bam_file_sorted}
    """
    
# # Perform peak calling on Jurkat h3k27me3 and ctcf bam files
# rule peak_calling:
#   input:
#     sorted_bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files_sorted/{bam_type}_sorted.bam"
#   output:
#     bed = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat.bed",
#     summits = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_summits.bed"
#   resources:
#     mem = "64G",
#     time = "4:00:00"
#   shell:
#     """
#     # Load necessary modules
#     ml system
#     ml biology
#     ml samtools
#     ml bedtools
#     ml py-macs2/2.2.9.1_py39
#     
#     # Make sure the output directory exists
#     mkdir -p $(dirname {output.bed})
# 
#     # Perform peak calling on the BAM file using MACS2
#     macs2 callpeak -t {input.sorted_bam} \
#                     --gsize 2.7e9 \
#                     -n jurkat \
#                     -f AUTO \
#                     -g hs \
#                     -p 0.1 \
#                     --call-summits \
#                     --outdir $(dirname {output.bed})
#                     
#     # Rename narrowPeak file to .bed
#     mv $(dirname {output.bed})/jurkat_peaks.narrowPeak {output.bed}
#     
#     # Retain only the top 150,000 peaks based on read count
#     # Count reads in each peak and retain the top 150K peaks by read count
#     bedtools coverage -a {output.bed} -b {input.sorted_bam} | \
#     sort -k7,7nr | \
#     head -n 150000 | \
#     cut -f1-6 > {output.bed}
#     """
    
rule peak_calling:
  input:
    sorted_bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files_sorted/{bam_type}_sorted.bam"
  output:
    bed = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_raw.bed",
    summits = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_summits.bed"
  resources:
    mem = "32G",
    time = "4:00:00"
  shell:
    """
    # Load necessary modules
    ml system
    ml biology
    ml samtools
    ml py-macs2/2.2.9.1_py39
    
    # Make sure the output directory exists
    mkdir -p $(dirname {output.bed})

    # Perform peak calling on the BAM file using MACS2
    macs2 callpeak -t {input.sorted_bam} \
                    --gsize 2.7e9 \
                    -n jurkat \
                    -f AUTO \
                    -g hs \
                    -p 0.1 \
                    --call-summits \
                    --outdir $(dirname {output.bed})
                    
    # Rename narrowPeak file to .bed
    mv $(dirname {output.bed})/jurkat_peaks.narrowPeak {output.bed}
    """
    
# Filter peaks to get top 150k by read count
rule filter_top_peaks:
  input:
    bed = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_raw.bed",
    bam = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/jurkat_bam_files_sorted/{bam_type}_sorted.bam"
  output:
    bed = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat.bed",
    coverage = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_coverage.bed",
    sorted = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/jurkat_sorted.bed"
  resources:
    mem = "88G",
    time = "2:00:00"
  shell:
    """
    # Load necessary modules
    ml system
    ml biology
    ml bedtools
    
    # Calculate coverage
    bedtools coverage -a {input.bed} -b {input.bam} > {output.coverage}
    
    # Sort by read count
    sort -k7,7nr {output.coverage} > {output.sorted}
    
    # Extract top 150K peaks and relevant columns
    head -n 150000 {output.sorted} | cut -f1-6 > {output.bed}
    """
    

# Rule to download BED files
rule download_bed_files:
  output:
    "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}.bed"
  params:
    lambda wildcards: config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]['download_bed_files'][wildcards.bam_type][wildcards.cell_type]
  resources:
    mem = "24G",
    time = "4:00:00"
  shell:
    """
    wget -O {output}.gz {params}
    gunzip -f {output}.gz
    """
    
    
# An intermediate file used in the ABC pipeline is used in downstream analysis for overlapping h3k27ac and dnase values with enhancers
# This file already existed for wtc11, k562, and gm12878
  # K562: /oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ABC/dnase_h3k27ac_intactHiC/K562/Neighborhoods/EnhancerList.txt
  # GM12878: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0425_GM12878_Extended/GM12878_DNase_H3K27ac_intact/Neighborhoods/EnhancerList.txt
  # WTC11: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0612_WTC11_for_EJ/wtc11/Neighborhoods/EnhancerList.txt
# Maya and I created this file for hct116 and jurkat using the processed h3k27ac I processed above for jurkat along with dnase from ENCODE
  #  "https://www.encodeproject.org/files/ENCFF361WBN/@@download/ENCFF361WBN.bam"
# And for hct116 by using the following two files
  # h3k27ac: "https://www.encodeproject.org/files/ENCFF340TPS/@@download/ENCFF340TPS.bam"
  # dnase: "https://www.encodeproject.org/files/ENCFF304MDQ/@@download/ENCFF304MDQ.bam"
# Jurkat and HCT116 outputs can be found here
  # HTC116: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_1029_for_JAG/HCT116/Neighborhoods/EnhancerList.txt
  # Jurkat: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G/results/2024_1029_for_JAG/Jurkat/Neighborhoods/EnhancerList.txt
# I then copied these intermediate files from the ABC pipeline into the resources folder in resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/ABC_outputs


### EQTL Questions
# CD14+ monocyte:
  # H3K27me3 experiment: https://www.encodeproject.org/experiments/ENCSR000ASK/
  # CTCF experiment: https://www.encodeproject.org/experiments/ENCSR000ATN/
  # EnhancerList: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_v1.0.0/ENCODE_rE2G/results/2024_0110_monocyte_T/CD14_pos_monocyte/Neighborhoods/EnhancerList.txt
# CD4+ T-cell:
  # H3K27me3 experiment: https://www.encodeproject.org/experiments/ENCSR043SBG/
  # CTCF experiment: https://www.encodeproject.org/experiments/ENCSR470KCE/
  # EnhancerList: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_v1.0.0/ENCODE_rE2G/results/2024_0110_monocyte_T/CD4_pos_T/Neighborhoods/EnhancerList.txt


rule overlap_h3k27me3_and_ctcf_with_enhancers:
  input:
    training_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    validation_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    unfiltered_k562_dc_tap_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt.gz",
    peak_bed_file = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}.bed", bam_type = ["ctcf", "h3k27me3"], cell_type = ["k562", "wtc11", "hct116", "gm12878", "jurkat"])
  output:
    training_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/h3k27me3_ctcf_overlap/training_expt_pred_merged_annot.txt",
    validation_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/h3k27me3_ctcf_overlap/validation_expt_pred_merged_annot.txt",
    unfiltered_k562_dc_tap_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/h3k27me3_ctcf_overlap/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/overlap_h3k27me3_and_ctcf_with_enhancers.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/overlap_h3k27me3_and_ctcf_with_enhancers.R"


rule overlap_dnase_and_h3k27ac_with_enhancers:
  input:
    training_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    validation_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    unfiltered_k562_dc_tap_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt.gz",
    enhancer_list = expand("resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/ABC_outputs/{cell_type}_EnhancerList.txt", cell_type = ["k562", "wtc11", "hct116", "gm12878", "jurkat"])
  output:
    training_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/dnase_h3k27ac_overlap/training_expt_pred_merged_annot.txt",
    validation_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/dnase_h3k27ac_overlap/validation_expt_pred_merged_annot.txt",
    unfiltered_k562_dc_tap_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/dnase_h3k27ac_overlap/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/overlap_dnase_and_h3k27ac_with_enhancers.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/overlap_dnase_and_h3k27ac_with_enhancers.R"



bam_types = ["h3k27me3_ctcf", "dnase_h3k27ac"]
rule combine_all_tracks:
  input:
    validation_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/validation_expt_pred_merged_annot.txt", bam_type = bam_types),
    training_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/training_expt_pred_merged_annot.txt", bam_type = bam_types),
    unfiltered_k562_dc_tap_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt", bam_type = bam_types)
  output:
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  params:
    model_threshold = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["combine_all_tracks"]["model_threshold"]
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/combine_all_tracks.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/combine_all_tracks.R"  





