library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Function to process a single sample
process_sample <- function(sample_id) {
  # Load txdb
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  # Read genome sequence
  genome <- readDNAStringSet("/scratch/cs/infantbiome/fuc1/chip_seq/data/ref_genome/fasta/hg19.fa")
  names(genome) <- sapply(strsplit(names(genome), " "), `[`, 1)

  # Read peak file for given sample_id
  peak_file <- paste0("/scratch/cs/infantbiome/fuc1/chip_seq/data/peak_files_ref/", sample_id, "/", sample_id, "_peaks.narrowPeak")
  peak <- readPeakFile(peak_file)

 
	# annotate peak
	anno_peak <- annotatePeak(peak, tssRegion=c(-2000,200), TxDb=txdb)
	anno_peak_df <- as.data.frame(anno_peak)
	
	
	promoter_peaks <- anno_peak_df[anno_peak_df$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)"), ]
	
	
	enhancer_peaks <- anno_peak_df[!anno_peak_df$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)"), ]
	
	
	promoter_gr <- GRanges(seqnames = promoter_peaks$seqnames,
	                       ranges = IRanges(start = promoter_peaks$start, end = promoter_peaks$end),
	                       strand = promoter_peaks$strand)
	
	enhancer_gr <- GRanges(seqnames = enhancer_peaks$seqnames,
	                       ranges = IRanges(start = enhancer_peaks$start, end = enhancer_peaks$end),
	                       strand = enhancer_peaks$strand)
  # Export promoter and enhancer GRanges objects for each sample
  output_promoter_file <- paste0("/scratch/cs/infantbiome/fuc1/chip_seq/data/promoter/", sample_id, "_pro.bed")
  output_enhancer_file <- paste0("/scratch/cs/infantbiome/fuc1/chip_seq/data/enhancer/", sample_id, "_enh.bed")

  export(promoter_gr, output_promoter_file)
  export(enhancer_gr, output_enhancer_file)
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
process_sample(sample_id)