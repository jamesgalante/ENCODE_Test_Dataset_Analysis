#!/bin/bash
# Guide to processing DNase-seq data for CRE identification

# 1. Download DNase-seq data from ENCODE (hg38)
echo "Step 1: Downloading DNase-seq data..."
wget https://www.encodeproject.org/files/ENCFF205FNC/@@download/ENCFF205FNC.bam

# 2. Sort the BAM file using sambamba
echo "Step 2: Sorting BAM file..."
sambamba sort -t 9 -m 9GB ENCFF205FNC.bam
# This produces ENCFF205FNC.sorted.bam and ENCFF205FNC.sorted.bam.bai

# 3. Call peaks using MACS2
echo "Step 3: Calling peaks with MACS2..."
macs2 callpeak \
  -t ENCFF205FNC.sorted.bam \
  -n sample1 \
  -f BAM \
  -g hs \
  -p .1 \
  --call-summits \
  --outdir ./ \
  2> call_peaks.log
# This produces sample1_peaks.narrowPeak and sample1_summits.bed

# 4. Download chromosome sizes file
echo "Step 4: Downloading chromosome sizes file..."
wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv

# 5. Create a properly sorted chromosome sizes file
echo "Step 5: Creating sorted chromosome sizes file..."
# Extract chromosome names from BAM header
samtools view -H ENCFF205FNC.sorted.bam | grep '^@SQ' | awk -F'\t' '{for(i=1;i<=NF;i++) if($i ~ /^SN:/) print substr($i, 4)}' > bam_chroms.txt

# Sort chromosome sizes file to match BAM header
awk 'NR==FNR{a[$1]=NR; next} {if (a[$1]!="") print a[$1] "\t" $0}' bam_chroms.txt GRCh38_EBV.chrom.sizes.tsv | sort -k1,1n | cut -f2- > sorted_GRCh38_EBV.chrom.sizes.tsv

# 6. Quantify DNase-seq reads in peaks
echo "Step 6: Quantifying DNase-seq reads in peaks..."
bedtools sort -faidx sorted_GRCh38_EBV.chrom.sizes.tsv -i sample1_peaks.narrowPeak | bedtools coverage -sorted -a stdin -b ENCFF205FNC.sorted.bam -g sorted_GRCh38_EBV.chrom.sizes.tsv -counts > sample1_peaks.counts

# 7. Create candidate CREs
echo "Step 7: Creating candidate CREs..."
bedtools sort -i sample1_peaks.counts -faidx sorted_GRCh38_EBV.chrom.sizes.tsv | \
bedtools merge -i stdin -c 11 -o max | \
sort -nr -k 4 | awk '(NR <= 150000)' | \
bedtools intersect -b stdin -a sample1_peaks.counts -wa | \
awk 'BEGIN {OFS="\t"} {print $1, $2 + $10, $2 + $10}' | \
bedtools slop -i stdin -b 250 -g sorted_GRCh38_EBV.chrom.sizes.tsv | \
bedtools sort -i stdin -faidx sorted_GRCh38_EBV.chrom.sizes.tsv | \
bedtools merge -i stdin | \
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3, "0", "."}' > sample1_candidate_cres.bed

echo "Processing completed successfully!"
echo "Output file: sample1_candidate_cres.bed"

# 8. Copy the final file into the parent directory for easy reference by each Rmd script
cp sample1_candidate_cres.bed ../