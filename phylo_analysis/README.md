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
