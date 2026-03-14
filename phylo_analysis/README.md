## Phylogenomic analysis
This contains some example analysis could be done for phylogenomic research with SNP data, including 
- SNP missingness and depth information
- Basic diversity indices
- OUTFLANK (to detect outliers that may have adaptive potential)
- PCA
- Neighbor-joining tree
- Structure (admixture) analysis
- DAPC
- Mantel and partial Mantel test (to test isolation by distance accounting for spatial autocorrelation)

### To prepare, despite respective R packages listed at the start of each file, you will need:
- The population file used in raw data processing (including different levels of population if applicable, such as regions, subregions, sites etc.)
- Sample location file that has coorinates in Decimal Degrees format
