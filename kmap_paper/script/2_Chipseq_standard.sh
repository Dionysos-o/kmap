#!/bin/bash

# Prompt for sample number
SAMPLE=$1

# Directories
BASE_DIR="./"
SRA_DIR="${BASE_DIR}sra_files_single/"
FASTQ_DIR="${BASE_DIR}fastq_files_single/${SAMPLE}/"
SAM_DIR="${BASE_DIR}sam_files_single/${SAMPLE}/"
BAM_DIR="${BASE_DIR}bam_files_single/${SAMPLE}/"
PEAK_DIR="${BASE_DIR}peak_files_single/${SAMPLE}/"
FASTA_DIR="${BASE_DIR}fasta_files_single/${SAMPLE}/"

# Create directories
mkdir -p $SRA_DIR $FASTQ_DIR $SAM_DIR $BAM_DIR $PEAK_DIR $FASTA_DIR

# Download SRA data
echo "Downloading SRA data for sample $SAMPLE..."
prefetch -O $SRA_DIR $SAMPLE

# Convert SRA to FASTQ
echo "Converting SRA to FASTQ..."
fasterq-dump --split-files -O $FASTQ_DIR $SRA_DIR$SAMPLE.sra

# Reference genome
HG19_REFERENCE="./hg19.fa"

# Align FASTQ to hg19
echo "Aligning FASTQ to hg19..."
bowtie2 -p 20 -x $HG19_REFERENCE -1 ${FASTQ_DIR}${SAMPLE}_1.fastq -2 ${FASTQ_DIR}${SAMPLE}_2.fastq | samtools view -Sb - > ${BAM_DIR}${SAMPLE}.bam

# Call peaks
echo "Calling peaks..."
macs2 callpeak -t ${BAM_DIR}${SAMPLE}.bam -f BAMPE -g hs -n ${SAMPLE} --outdir $PEAK_DIR

# Convert narrowPeak to FASTA
echo "Converting narrowPeak to FASTA..."
bedtools getfasta -fi $HG19_REFERENCE -bed ${PEAK_DIR}${SAMPLE}_peaks.narrowPeak -fo ${FASTA_DIR}${SAMPLE}.fasta

echo "Pipeline for $SAMPLE completed successfully!"
