This project contains commands and scripts that haven been generated and used for Queuosine gene/protein analyses in the manuscript:   
#### "Comparative Genomics Analyses of Bacterial tRNA Modification Protein Families Using Clusters of Orthologous Genes (COG) and KEGG Orthology (KO) Classification Systems"   
Authors: Colbie J. Reed, Yifeng Yuan, and Valérie de Crécy-Lagard   

### Dependencies   
python v3.13 https://www.python.org/   
NCBI dataset https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/   
checkM2 https://github.com/chklovski/CheckM2   
MAFFT https://mafft.cbrc.jp/alignment/server/index.html   
conda packages: 
    gtdbtk-2.4.0   
    matplotlib>=3.9


### Usage   
1. phylogenetic analysis

```bash
# 1.1 Download genomes
cd tree
download_genomes_datasets.sh data/1322accessions.txt

# move downloaded genome fasta files to a folder. Replace soft link with fna files.

# 1.2 Quality check
dir_genome=${folder of the genomes}

# Runs CheckM2 on your genome folder and keeps the genomes with at least 90% completeness, at most 5% contamination, and at least 0.1 Mb (100,000 bp) total length.
checkm2_filter.sh -g ${dir_genome} -o checkm2_out -t 2

# 1.3 build tree
gtdbtk_fasttree_wag.sh checkm2_out
#output is tree/gtdbtk_out/bac120.fasttree.unrooted.tree
```

The output tree (tree/bac120.fasttree.unrooted.tree as demo) can be visualized in iTOL v6.8 (https://itol.embl.de/) using the annotation files:   
tree/taxid_tgtbinary_TrmAFO_thick.txt (show the presence of COG and KO of TrmA, TrmFO )   
tree/taxid_tgtbinary_mnmEGA_thick.txt (show the presence of COG and KO of MnmE, MnmG, and MnmA)   
tree/treelabel_phylum.txt   

2. COG and KO concordance   
``` bash
cd concordance

# calculate the prevalence of COG or KO per phylum
# Cell value (%) = (number of genomes in that phylum with the key present) / (total number of genomes in the matrix) * 100
python presence_by_phylum.py COG_AP.txt KO_AP.txt genome_phylum.txt results

# calculate the concordance of COG and KO per genome. Build a genome x (COG|KO pair) agreement matrix.
# 1=jointly presence or absence; 0= disagreement between COG and KO presence.
# For pairs where only one key is in the matrices: the column holds that key's own presence (1) / absence (0); the missing counterpart is not counted.
python cog_ko_agreement.py COG_AP.txt KO_AP.txt COG_KO_pair.txt results/concordance.tsv

# Explore COG-KO concordance across genomes and phyla.
# phyla with fewer genomes are pooled as "Other" in plots.
python cog_ko_explore.py COG_AP.txt KO_AP.txt COG_KO_pair.txt genome_phylum.txt results

python pipeline_cog_ko_agreement.py

```

3. paralog
```
# 3.1 convert list of locus rags to matrix (genome x locus_tag)
cd locus_tag

python data/COG.txt COG_tag.txt
python data/KO.txt KO_tag.txt

# 3.2 paralog stats
# with COG_count.txt, KO_count.txt, genome_phylum.txt in the same directory
cd ../paralog

run_count_distribution.sh 

# output: results/_count_distribution.tsv, results/count_distribution_by_phylum.tsv, results/_phylum_groups.tsv, results/_count_distribution_heatmap.png, results/_count_distribution_by_phylum_heatmap.png
# results/COG_mean_copy_by_phylum.tsv, results/KO_mean_copy_by_phylum.tsv, results/COG_KO_mean_copy_heatmap.png

# calculate correlation between KO prevalence and copy-number-when-present across KOs:
python script_correlation_dosage_prevalence.py
# output: ko_copy_prevalence/output correlation_dosage_prevalence.tsv

# KO copy-number vs prevalence analysis across bacterial genomes.
python script_ko_prevalence.py
# output: ko_copy_prevalence/fig1_overall_vs_phylum_mean.png, output: ko_copy_prevalence/fig2_top_heterogeneous.png, ko_copy_prevalence/ko_prevalence_copynumber.tsv, ko_copy_prevalence/ko_prevalence_by_phylum.tsv

# generate figure 3. KO prevalence vs copy-number-when-present are uncorrelated.
python script_correlation_dosage_prevalence.py
# output: ko_copy_prevalence/fig3_dosage_vs_prevalence.svg/.png and ko_copy_prevalence/fig3_dosage_vs_prevalence_points.tsv

# OLS regression equation for fig3 (dosage-when-present ~ prevalence).
python script_ols_equation.py
# output: printed equation + ols_equation_dosage_prevalence.tsv in ko_copy_prevalence/
```

4. Genome size correlation

```bash
cd genome_size

# Correlation between tRNA-modification KO presence and genome size.
python script_ko_size_correlation.py
# output: in the ko_size_correlation folder: ko_size_correlation.tsv, sum_ko_size_correlation.tsv, correlation_linear_vs_log_long.tsv, correlation_linear_vs_log_wide.tsv, ko_size_correlation_by_phylum.tsv, fig1_sum_vs_size.png, fig2_per_ko_correlations.png, fig3_linear_vs_log_examples.png.
```

5. Traits ~ modification correlation

# Analysis scripts — tRNA modification presence/absence vs genome traits and between modifications

All scripts are self-contained R scripts (R 4.4.2; phylolm 2.6.5, ape 5.8.1, phytools 2.5.2, ggplot2 4.0.3).
Inputs are read from `/mnt/user-uploads/input_v4/` (`mod_AP.txt`, `genome_traits.txt`, `select_trait.txt`,
`bac120.fasttree.unrooted.tree`) and, for the Aquificota-rooted scripts, `/mnt/user-uploads/root.txt`.
Run each with `Rscript <script>.R`. Outputs are written to the `/mnt/results/` folders listed below.

Shared design elements across all scripts:
- QC gates stop the run on any mismatch (matrix dimensions, 0/1 values, genome-ID agreement, trait counts).
- Presence filter: exclude modifications present in < 5% or > 95% of the analysed genomes (fraction rule).
- Primary test: phyloglm, logistic MPLE, btol = 10; OR with 95% Wald CI; BH correction per trait/pair family.
- Effect size: phylolm (model = lambda), r_phylo = t / sqrt(t^2 + n - 2); non-gated, separate BH family.
- Refit protocol: flagged fits (non-convergence, alphaWarn != 0, warnings, |slope| > 5, missing coef/SE)
  are refitted from the default start plus 9 random starts (set.seed(123) immediately before the draws;
  coefficients ~ U(-2, 2) drawn first, then alpha uniform on the log scale of the search range); the fit
  with the highest penalized log-likelihood among converged starts is kept.
- Flags file lists every tested/untested unit with reasons and result flags.

## Scripts

| # | Script | Analysis | Rooting | Outputs |
|---|--------|----------|---------|---------|
| 1 | `01_mod_as_response_midpoint.R` | Binary trait ~ modification, modification as response (mod ~ trait) | midpoint (per pruned trait tree) | `trait_mod_correlation/` |
| 2 | `02_trait_as_response_midpoint.R` | Binary trait ~ modification, trait as response (trait ~ mod) | midpoint (per pruned trait tree) | `trait_mod_correlation_v2_trait_response/` |
| 3 | `03_trait_mod_both_orientations_aquificota_root.R` | Both orientations (a: mod ~ trait, b: trait ~ mod) in one run | Aquificota outgroup (per pruned trait tree) | `trait_mod_correlation_v2_aquificota_root/<orientation>/<trait>/` |
| 4 | `04_mod_mod_midpoint.R` | Modification ~ modification, all 630 unordered pairs (336 testable) | midpoint (one tree for all pairs) | `mod_mod_correlation/` |
| 5 | `05_mod_mod_aquificota_root.R` | Modification ~ modification, identical design | Aquificota outgroup (one tree) | `mod_mod_correlation_aquificota_root/` |

Notes:
- Scripts 1 and 2 differ only in the model orientation; running both implements the
  "switch predictor and response" robustness check. Script 3 reproduces both orientations
  under the alternative rooting in a single execution.
- Scripts 4 and 5 differ only in the rooting; `sensitivity_table.tsv` and
  `sensitivity_forest.png/svg` in `mod_mod_correlation/` compare their significant pairs
  (generated by `sensitivity_table_fig.R`).
- Scripts 1-3 prune the tree to the genomes annotated for each trait before rooting, so the
  rooting is performed once per trait tree; scripts 4-5 use all 1,322 genomes and root once.
- Each script writes a `flags.tsv` (all units, tested yes/no, untested reasons, result flags)
  and prints package versions at completion.


### Help and Issues  
Please contact Yifeng Yuan at yuanyifeng@ufl.edu or Colbie J. Reed at colbie.j.reed@gmail.com   

###  Authors
Colbie J. Reed, Ph.D, Yifeng Yuan, Ph.D., and Valerie de Crecy-Lagard (Principal Investigator)

Version History v1.0.0 --2026/09/022 --creation of the project.
