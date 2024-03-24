# Kmer Manifold Approximation and Projection (KMAP)
kmap is a package for visualizing kmers in 2D space. 
![image](https://github.com/Dionysos-o/kmap/blob/main/kmap_paper/kmapshownew.gif)
## Quick start
1，For installation of kmap, you can refer to :
2, make sure you have installed the following package in your enviroment
requirements:
```bash
numpy==1.16.4
pandas==0.24.2
matplotlib==3.1.0
seaborn==0.9.0
scikit-learn==0.21.2
```

## KMAP experiment
### HT-selex anylysis
Please ensure that the MEME suite has been installed in your environment by following the tutorial [here](https://meme-suite.org/meme/doc/install.html). We recommend using Docker for the installation.

You can use script `0_ht_selex_downlaod.sh` to download the HT-SELEX data from European Nucleotide Archive (ENA) under the accession number ERP001824 [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB3289)
The file is fasta format and you can feed it to kmap directly without any preprocessing. 

```bash
sh 1_ht_selex_kmap.sh
sh 1_ht_selex_meme.sh
```
calculate the similarity score of MEME and KMAP result
```bash
sh 1_cal_similarity.sh
```
use the following command to run the experiment with DREME and STREME
```bash
sh 1_ht_selex_dreme_streme.sh
```

### EWS Chipseq analysis
In Ewing Sarcoma data analysis, we use the source data from [GSE47753](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753).

You can follow the steps below to conduct thh Chip-seq pipline and get the enrichment fasta file for KMAP。Here we use SRR935526 as an example:
```bash
sh 2_Chipseq_standard.sh SRR935526
```
#### Analysis of H3K27ac CHip-seq in promoter and enhancer regions

Get the annotation of the peak, you can download the annotation file from [here](https://www.bioconductor.org/packages/release/data/annotation/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html)
```bash
 txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
 anno_peak <- annotatePeak(peak, tssRegion=c(-2000,200), TxDb=txdb)
```
Divide the peak into promoter and enhancer(we define the promoter region as TSS:(-1kb, +2kb), the enhancer region is the rest))
```bash
anno_peak_df <- as.data.frame(anno_peak)
promoter_peaks <- anno_peak_df[anno_peak_df$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)"), ]
enhancer_peaks <- anno_peak_df[!anno_peak_df$annotation %in% c("Promoter (<=1kb)", "Promoter (1-2kb)"), ]
```
Export the peak to bed file 
```bash
promoter_gr <- GRanges(seqnames = promoter_peaks$seqnames,
	                       ranges = IRanges(start = promoter_peaks$start, end = promoter_peaks$end),
	                       strand = promoter_peaks$strand)
export(promoter_gr, output_promoter_file)
```
or simply use the 
```bash
sh 2_spilt_promoter_enhancer.sh
```
### Gene edtting data
use the following command to do the anylysis of gene editing data:
```bash
sh 3_ge_align.sh
python 3_kmap_ge.py
```

[comment]: <> (Release commands)
[comment]: <> (python -m build) 
[comment]: <> (python3 -m twine upload --repository testpypi dist/*)

