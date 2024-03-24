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
(make sure your envrioment has been installed with the MEME suit, please refer to the tutotial. we recommand you to use docker way to install it)
You can use script `ht_selex_downlaod.sh` to download the HT-SELEX data from European Nucleotide Archive (ENA) under the accession number ERP001824 [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753)
The file is fasta format and you can feed it to kmap directly without any preprocessing. 

```bash
sh 1_ht_selex_kmap.sh
sh 1_ht_selex_meme.sh
```
calculate the similarity score of MEME and KMAP result
```bash
1_cal_similarity.sh
```
use the following command to run the experiment with DREME and STREME
```bash
sh 1_ht_selex_dreme_streme.sh
```

### EWS Chipseq analysis
- **Step 0**: confirm kmap is successfully installed
```bash
kmap --help
```
- **Step 1**: create a directory `test` and save the test fasta file [test.fa](./tests/test.fa) in this directory


### Auxiliary functions

`

[comment]: <> (Release commands)
[comment]: <> (python -m build) 
[comment]: <> (python3 -m twine upload --repository testpypi dist/*)

