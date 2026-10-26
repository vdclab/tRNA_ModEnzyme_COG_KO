### Description  
This project contains commands and scripts that haven been generated and used for Queuosine gene/protein analyses in the manuscript:   
#### "Comparative Genomics Analyses of Bacterial tRNA Modification Protein Families Using Clusters of Orthologous Genes (COG) and KEGG Orthology (KO) Classification Systems"   
Authors: Colbie J. Reed, Yifeng Yuan, and Valérie de Crécy-Lagard   

### Dependencies   
python v3.13 https://www.python.org/   
MAFFT https://mafft.cbrc.jp/alignment/server/index.html   
conda packages: gtdbtk-2.4.0   

### Usage   
1. Download genomes

```bash
cd tree
download_genomes_datasets.sh data/1322accessions.txt

# move genome fasta files to a folder
dir_genome=${folder of the genomes}

# Runs CheckM2 on your genome folder and keeps the genomes with at least 90% completeness, at most 5% contamination, and at least 0.1 Mb (100,000 bp) total length.
checkm2_filter.sh -g ${dir_genome} -o checkm2_out -t 2
```

2. xxxx

### Help and Issues  
Please contact Yifeng Yuan at yuanyifeng@ufl.edu or Colbie J. Reed at colbie.j.reed@gmail.com   

###  Authors
Colbie J. Reed, Ph.D, Yifeng Yuan, Ph.D., and Valerie de Crecy-Lagard (Principal Investigator)

Version History v1.0.0 --2026/09/022 --creation of the project.
