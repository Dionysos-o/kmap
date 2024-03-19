#!/bin/bash

# ?~O~P示?~T??~H??~S?~E??| ??~\??~O?
SAMPLE=$1

# SRA?~V~G件路?~D
BAM_DIR="./bam_files/"
PEAK_DIR="./peak_files_ref/${SAMPLE}"
FASTA_DIR="./fasta_files_ref/${SAMPLE}_ref"

mkdir -p $SRA_FILE
mkdir -p $FASTQ_DIR
mkdir -p $SAM_DIR
mkdir -p $BAM_DIR
mkdir -p $PEAK_DIR
mkdir -p $FASTA_DIR


HG19_REFERENCE="./hg19.fa"



# ?~C?~T?peaks?~L?~N??~WnarrowPeak?~V~G件
echo "Calling peaks..."

macs2 callpeak -t $BAM_DIR/${SAMPLE}.bam -c $BAM_DIR/SRR15358917.bam -f BAMPE -g hs -n ${SAMPLE} --outdir $PEAK_DIR/

# narrowPeak?~V~G件路?~D
NARROWPEAK="peaks_peaks.narrowPeak"

# 使?~T?bedtools?~N??~O~Vfasta?~O?~H~W
echo "Converting narrowPeak to FASTA..."
bedtools getfasta -fi $HG19_REFERENCE -bed $PEAK_DIR/${SAMPLE}"_peaks.narrowPeak" -fo $FASTA_DIR/${SAMPLE}_ref.fasta

echo "Done!"
