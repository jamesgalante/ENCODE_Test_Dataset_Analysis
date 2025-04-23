# This rule file processes the expt_pred_merged_annot outputs for the subset_upsampling_analysis

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
    bam = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/ctcf_bam_files/jurkat.bam"
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
    
    mkdir -p resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/
    
    # Loop through the SRR IDs, download the FASTQ files
    for srr_id in {params.srr_ids}; do
        # Step 1: Download paired-end FASTQ files from SRA with original headers
        fasterq-dump $srr_id -O resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/ --split-files -F --gzip
    done
    
    for srr_id in {params.srr_ids}; do
      # Align each single-end FASTQ file separately
      bwa mem -t {threads} {input.index} resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.fastq > \
      resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.sam
      
      # Convert each SAM file to BAM
      samtools view -bS resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.sam > \
      resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.bam
    done
    
    # Merge the BAM files into one BAM file
    samtools merge -f -@ {threads} {output.bam} $(for srr_id in {params.srr_ids}; do echo resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_ctcf_chip/${{srr_id}}.bam; done) 
    """

# Process H3K27ac ChIP seq to be used in the analysis
rule process_jurkat_h3k27ac_chip:
  input:
    index = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa"
  output:
    bam = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/h3k27ac_bam_files/jurkat.bam"
  params:
    srr_ids = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["process_jurkat_h3k27ac_chip"]["srr_ids"]
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
    
    mkdir -p resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/
    
    # Download FASTQ files for each SRR ID
    for srr_id in {params.srr_ids}; do
        fasterq-dump $srr_id -O resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/ --split-files -F --gzip
    done
    
    for srr_id in {params.srr_ids}; do
      # Align each single-end FASTQ file separately
      bwa mem -t {threads} {input.index} resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.fastq > \
      resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.sam
      
      # Convert each SAM file to BAM
      samtools view -bS resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.sam > \
      resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.bam
    done
    
    # Merge the BAM files into one BAM file
    samtools merge -f -@ {threads} resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_aligned.bam \
    $(for srr_id in {params.srr_ids}; do echo resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/${{srr_id}}.bam; done) 

    # https://www.biostars.org/p/318974/ has a discussion on deduplicating ChIP-seq data
    # https://www.biostars.org/p/390305/ has a discussion mentioning markdup being a sufficient replacement for picard
    # https://www.htslib.org/doc/samtools-markdup.html documentation on markdup
    # Run samtools fixmate (required for markdup)
    samtools fixmate -m resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_aligned.bam \
        resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_fixmate.bam
    
    # Sort the BAM file before removing duplicates
    samtools sort -@ {threads} -o resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_sorted.bam \
        resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_fixmate.bam
    
    # Remove PCR duplicates using samtools markdup
    samtools markdup -r resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27ac_chip/jurkat_h3k27ac_sorted.bam {output.bam}
    """
    
# Process H3K27ac ChIP seq to be used in the analysis
rule process_jurkat_h3k27me3_chip:
  input:
    index = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/hg38/hg38.fa"
  output:
    bam = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/h3k27me3_bam_files/jurkat.bam"
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
    
    mkdir -p resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/

    # Download FASTQ files for the indicated SRR ID
    fasterq-dump {params.srr_id} -O resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/ --split-files -F --gzip

    # Align the single-end FASTQ file
    bwa mem -t {threads} {input.index} resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.fastq > \
    resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.sam

    # Convert the SAM file to BAM
    samtools view -bS resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.sam > \
    resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.bam

    # Run samtools fixmate (required for markdup)
    samtools fixmate -m resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}.bam \
        resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_fixmate.bam

    # Sort the BAM file before removing duplicates
    samtools sort -@ {threads} -o resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_sorted.bam \
        resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_fixmate.bam

    # Remove PCR duplicates using samtools markdup
    samtools markdup -r resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/process_jurkat_h3k27me3_chip/{params.srr_id}_sorted.bam {output.bam}
    """

# This rule will take links from the config file for each of the bam files and download the bam file
# Rule to download BAM files
rule download_bam_files:
  output:
    "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_bam_files/{cell_type}.bam"
  params:
    lambda wildcards: config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]['download_bam_files'][wildcards.bam_type][wildcards.cell_type]
  resources:
    mem = "24G",
    time = "4:00:00"
  shell:
    """
    wget -O {output} {params}
    """

# Process the bam file - sort and filter
rule process_bam_files:
  input:
    bam_file = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_bam_files/{cell_type}.bam"
  output:
    bam_file_sorted = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_bam_files_sorted/{cell_type}_sorted.bam"
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
    
    
rule peak_calling:
  input:
    sorted_bam = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_bam_files_sorted/{cell_type}_sorted.bam"
  output:
    narrowpeak = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}_peaks.narrowPeak",
    summits = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}_summits.bed"
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
    mkdir -p $(dirname {output.narrowpeak})

    # Perform peak calling on the BAM file using MACS2
    macs2 callpeak -t {input.sorted_bam} \
                   --gsize 2.7e9 \
                   -n {wildcards.cell_type} \
                   --outdir $(dirname {output.narrowpeak})
    """
    
    
rule compute_peak_bam_overlap:
  input:
    peaks = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}_peaks.narrowPeak",
    bam = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_bam_files_sorted/{cell_type}_sorted.bam"
  output:
    peak_counts = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}_peak_counts.bed"
  resources:
    # High memory not necessary for most things, but wtc11 required high memory
    mem = "64G",
    time = "2:00:00"
  shell:
    """
    # Load necessary modules
    ml load system
    ml load biology
    ml load bedtools/2.30.0
    ml load samtools
    
    # Define a unique temporary file name based on bam_type and cell_type
    tmp_chr_order=tmp_chr_order_{wildcards.bam_type}_{wildcards.cell_type}.txt

    # Extract chromosome names and lengths from BAM file and create genome file
    samtools view -H {input.bam} | grep -P '^@SQ' | cut -f 2,3 | \
    awk 'BEGIN{{OFS="\t"}}{{split($1, a, ":"); split($2, b, ":"); print a[2], b[2]}}' > $tmp_chr_order

    # Sort the narrowPeak file using the extracted chromosome order and directly pipe it into bedtools coverage
    bedtools sort -faidx $tmp_chr_order -i {input.peaks} | \
    bedtools coverage -sorted -g $tmp_chr_order -a stdin -b {input.bam} > {output.peak_counts}

    # Clean up the temporary chromosome order file
    rm $tmp_chr_order
    """
    
    
# Rule to connect validation analyzing to all other parts of the pipeline
# rule copy_expt_pred_merged_annot 
  # to the analyze validation sets resources folder
# Currently, we're just copying the expt_pred_merged_annot from the crispr comparison pipeline manually
# Also have to create the unfiltered_k562_dc_tap expt_pred_merged_annot file from the main pipeline (or copy it over rather)

# Overlap the bam files with the enhancers - different for each assay (e.g. H3K27ac peaks get extended)
rule overlap_bam_with_enhancers:
  input:
    training_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    validation_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    unfiltered_k562_dc_tap_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt.gz",
    peak_counts_files = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{{bam_type}}_peaks/{cell_type}_peak_counts.bed", cell_type = ["k562", "wtc11", "hct116", "gm12878", "jurkat"])
  output:
    training_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/training_expt_pred_merged_annot.txt",
    validation_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/validation_expt_pred_merged_annot.txt",
    unfiltered_k562_dc_tap_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/overlap_bam_{bam_type}_with_enhancers.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/OLD_overlap_bam_with_enhancers.R"
    

# Combine the bam overlapped files into the final expt_pred_merged_annot_files
bam_types = ["ctcf", "h3k27ac", "h3k27me3", "dnase"]
rule combine_overlapped_bam_enhancer_files_and_format:
  input:
    validation_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/validation_expt_pred_merged_annot.txt", bam_type = bam_types),
    training_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/training_expt_pred_merged_annot.txt", bam_type = bam_types),
    unfiltered_k562_dc_tap_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/unfiltered_k562_dc_tap_expt_pred_merged_annot.txt", bam_type = bam_types)
  output:
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_unfiltered_k562_dc_tap_expt_pred_merged_annot.txt"
  params:
    #model_threshold = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["combine_overlapped_bam_enhancer_files_and_format"]["model_threshold"]
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/combine_overlapped_bam_enhancer_files_and_format.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/OLD_combine_overlapped_bam_enhancer_files_and_format.R"


# This is an endpoint in the pipeline =========================================================






