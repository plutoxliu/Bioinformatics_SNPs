## Overview
This repository contains detailed steps to analyse Genotyping-by-Sequencing (GBS) data.<br>
For usage please see detailed comments in each file <br>
I would appretiate acknowledgement or citation (Liu et al. 2026 https://doi.org/10.1111/jpy.70220) if using this repository
<br><br>
## Before analyses
- Genetic extractions and GBS library preparations steps are not described here, as they are subject to individual protocols and species. Deails for our preparations will be uploaded once the manuscript is accepted.
- Sequencing: I used Illumina NextSeq 2000 with 150 bp paired-end reads at the Otago Genomics Facility (University of Otago, New Zealand).
<br><br>
## Process of raw data
As I used customised adapter plates, raw reads need to be cleaned and demultiplexed. <br>
All codes were writen as batch jobs submitted to New Zealand eScience Infrastructure (NeSI), therefore users need to adapt them for their own use

### Prepare: 
- a codes.txt file that matches individual adapter codes to individual samples
- a sample.txt file that containes individual sample name
- a pop.txt file that matches individuals to their respective populations
- A reference genome (if available)

### Packages used:
- Stacks2: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775. Manual: https://catchenlab.life.illinois.edu/stacks/manual/
- Cutadapt: https://journal.embnet.org/index.php/embnetjournal/article/view/200/479. Manual: https://cutadapt.readthedocs.io/en/stable/
- BWA: https://academic.oup.com/bioinformatics/article/25/14/1754/225615. Manual: https://bio-bwa.sourceforge.net/bwa.shtml
- SAMrools: https://academic.oup.com/gigascience/article/10/2/giab008/6137722. Manual: https://www.htslib.org/doc/#manual-pages
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC: https://docs.seqera.io/multiqc

### Steps:
#### 01_cut.sh 
remove Illumina adapters and poly-G sequences
#### 02_demulti.sh 
Demultiplexing
#### 03_concat.py
concatenate forward and reverse reads into one file for each individual sample
#### 04_trimsize.sh
trim sequences to desired length, usually need to use fastqc and multiqc before to check qualities etc. before performing this step
#### 05_align.sh
algin sequences to a provided reference genome, after this step needs to clean up low quality samples, the easiest way is to check file size, but can also use Stacks
#### 06_refmap.sh
SNP calling
#### 07_SNP_filter
there are a lot of methods to filter SNPs including using Stacks or R packages. A reference code using Stacks is included.<br>
But I mostly used "SNPfiltR" package (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13618) as it's interactive and visualises different parameters. It has very detailed step-by-step guide on its website (https://devonderaad.github.io/SNPfiltR/). <br>
I only changes some parameters tailored to my SNP data. Therefore detailed R scripts not attached here. My filtering information can be found in the manuscipt once accepted.

## Phylogenomic analysis
This contains some example analysis could be done for phylogenomic research with SNP data, including 
### 01_data_QC_diversity
1.  Data Import & Population Assignment
2.  Sequencing Depth Visualisation
3.  Missingness per Sample
4.  GDS File Creation & Diversity Metrics Setup
5.  Fst Estimation (Global & Pairwise)
6.  Heterozygosity & Diversity Indices
7.  OutFLANK – Outlier Loci Detection
8.  Neighbour-Joining Tree   
### 02_PCA_AMOVA
1. Principal Components Analysis (PCA)
2. AMOVA (Analysis of Molecular Variance)
### 03_structure_DAPC
1.  Structure Analysis (LEA / sNMF)
2.  DAPC (Discriminant Analysis of Principal Components)
### 04_mantel_partial_mantel
1. Mantel Test (Isolation by Distance)
2. Partial Mantel Test (IBD controlling for sampling density)

### To prepare, despite respective R packages listed at the start of each file, you will need:
- The population file used in raw data processing (including different levels of population if applicable, such as regions, subregions, sites etc.)
- Sample location file that has coordinates in Decimal Degrees format
