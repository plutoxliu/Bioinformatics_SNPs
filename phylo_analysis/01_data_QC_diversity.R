################################################################################
##                                                                            ##
##              POPULATION GENOMICS — DATA IMPORT, QC & DIVERSITY            ##
##                                                                            ##
##  Sections:                                                                 ##
##    1.  Data Import & Population Assignment                                 ##
##    2.  Sequencing Depth Visualisation                                      ##
##    3.  Missingness per Sample                                              ##
##    4.  GDS File Creation & Diversity Metrics Setup                         ##
##    5.  Fst Estimation (Global & Pairwise)                                  ##
##    6.  Heterozygosity & Diversity Indices                                  ##
##    7.  OutFLANK – Outlier Loci Detection                                   ##
##    8.  Neighbour-Joining Tree                                              ##
##                                                                            ##
##  Output objects passed to other scripts:                                   ##
##    vcf, GBS, samples, D.ind, D.pop, D.ind.dist, D.pop.dist                ##
##    (save with save() / load() or re-run this script first)                 ##
##                                                                            ##
################################################################################


# ==============================================================================
# LIBRARIES
# ==============================================================================

library(vcfR)
library(adegenet)
library(ape)
library(dartR)
library(mmod)
library(reshape2)
library(tidyverse)
library(hierfstat)
library(stringr)
library(ggplot2)
library(StAMPP)
library(SNPRelate)   # required for snpgdsOpen / snpgdsFst


################################################################################
# SECTION 1: DATA IMPORT & POPULATION ASSIGNMENT
################################################################################

# ------------------------------------------------------------------------------
# Read VCF and convert to genlight
# ------------------------------------------------------------------------------

vcf <- read.vcfR("RE50.vcf")

head(vcf)
head(getFIX(vcf))

# Convert to genlight format (required by dartR / adegenet functions)
GBS <- vcfR2genlight(vcf)
length(GBS@ind.names)   # check number of individuals loaded

# ------------------------------------------------------------------------------
# Build sample metadata table
# Sample names are structured as: Region_Site_IndividualNumber
# We extract Region and Site by splitting on "_"
# ------------------------------------------------------------------------------

samples <- data.frame(GBS@ind.names)
samples$Individual <- samples$GBS.ind.names
samples$Region     <- word(samples$Individual, 1, sep = "_")   # first field
samples$Site       <- word(samples$Individual, 2, sep = "_")   # second field

samples$Region <- as.factor(samples$Region)
samples$Site   <- as.factor(samples$Site)

# NOTE: This data frame contains:
#   GBS.ind.names – original sequence names from the VCF
#   Individual    – same as above (set during raw-data processing)
#   Region        – broad grouping variable
#   Site          – finer-scale grouping variable
# IMPORTANT: Never reorder the rows of this data frame — row order must
# match the individual order in the genlight object throughout.

# ------------------------------------------------------------------------------
# Assign population labels to the genlight object
# Region is used as the primary population variable
# ------------------------------------------------------------------------------

GBS@pop  <- samples$Region
pop(GBS) <- GBS@pop
pop(GBS)   # confirm assignments


################################################################################
# SECTION 2: SEQUENCING DEPTH VISUALISATION
################################################################################

# Extract per-sample, per-locus read depth from the DP format field
dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
dp[1:4, 1:6]   # quick sanity check

# Reshape to long format for ggplot
dpf <- melt(dp, varnames = c("Index", "Sample"),
            value.name = "Depth", na.rm = TRUE)

# Colour palette – one colour per sample, recycled as needed.
# NOTE: Change the recycling length (349) to match your actual sample count.
palette_depth <- rep_len(
  c("#FF0000", "#FF6E00", "#FFC300", "#FFFF00", "#AAD500",
    "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082",
    "#812BA6", "#B857CA", "#D03A87"),
  349
)

# Boxplot of depth per sample on a log2 scale
ggplot(dpf, aes(x = Sample, y = Depth)) +
  geom_boxplot(fill = palette_depth) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(
    trans        = scales::log2_trans(),
    expand       = c(0, 0),
    breaks       = c(1, 10, 100, 1000, 5000),
    minor_breaks = c(1:10, 2:10 * 10, 2:10 * 100, 2:5 * 1000)
  ) +
  theme(
    panel.grid.major.y = element_line(color = "#A9A9A9", linewidth = 0.6),
    panel.grid.minor.y = element_line(color = "#C0C0C0", linewidth = 0.2)
  )


################################################################################
# SECTION 3: MISSINGNESS PER SAMPLE
################################################################################

# Calculate the proportion of missing genotypes for each sample
myMiss           <- apply(dp, MARGIN = 2, function(x) sum(is.na(x)))
myMiss           <- myMiss / nrow(vcf)
myMiss           <- data.frame(levels(dpf$Sample), myMiss)
colnames(myMiss) <- c("Sample", "Missing")

# Bar chart of missingness per sample
# NOTE: palette_depth length should match sample count (see Section 2 note)
ggplot(myMiss, aes(x = Sample, y = Missing)) +
  geom_col(fill = palette_depth) +
  theme_bw() +
  labs(y = "Missingness (%)") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(angle = 60, hjust = 1)
  ) +
  scale_y_continuous(expand = c(0, 0))


################################################################################
# SECTION 4: GDS FILE CREATION & DIVERSITY METRICS SETUP
################################################################################

# dartR and SNPRelate functions operate on a GDS (Genomic Data Structure) file.
# Here we prepare the genlight metadata and export to GDS format.

# Build locus metrics table (chromosome and position)
metrics_df <- data.frame(
  snp_chr = as.character(GBS@chromosome),
  snp_pos = GBS@position,
  stringsAsFactors = FALSE
)
GBS@other$loc.metrics <- metrics_df

# Assign hierarchical strata (Region = broad, Site = fine-scale)
strata_df   <- data.frame(Pop = samples$Region, Subpop = samples$Site)
strata(GBS) <- strata_df
setPop(GBS) <- ~Pop

# Store population labels inside the ind.metrics slot (used by some dartR functions)
GBS@other$ind.metrics$population    <- strata_df$Pop
GBS@other$ind.metrics$subpopulation <- strata_df$Subpop

# SNPRelate requires numeric chromosome IDs, so convert scaffold names
scaffold_names   <- as.character(GBS@chromosome)
scaffold_numeric <- as.numeric(as.factor(scaffold_names))
GBS@other$loc.metrics$chrom_num <- scaffold_numeric
GBS@other$loc.metrics$pos       <- GBS@position

# Export genlight to GDS format
gl2gds(GBS,
       outfile = "ALL.gds",
       outpath = getwd(),
       snp_pos = "snp_pos",
       snp_chr = "chrom_num")


################################################################################
# SECTION 5: FST ESTIMATION (GLOBAL & PAIRWISE)
################################################################################

# ------------------------------------------------------------------------------
# 5a. Global Fst via SNPRelate (Weir & Cockerham 1984)
# ------------------------------------------------------------------------------

GDS         <- snpgdsOpen("ALL.gds", readonly = FALSE)
samp_folder <- addfolder.gdsn(GDS, "sample.annot")

# Write population codes into the GDS annotation folder
add.gdsn(samp_folder, "pop",    val = as.character(pop(GBS)))
add.gdsn(samp_folder, "subpop", val = as.character(GBS@strata$Subpop))

# Read back population codes to confirm
pop_code <- read.gdsn(index.gdsn(GDS, path = "sample.annot/pop"))
table(pop_code)

# Compute Fst (W&C84 method)
# sample.id = NULL includes all samples
fst_result <- snpgdsFst(GDS,
                        sample.id  = NULL,
                        population = as.factor(pop_code),
                        method     = "W&C84")

fst_result$Fst              # Weighted Fst across all SNPs
fst_result$MeanFst          # Mean Fst across all SNPs
summary(fst_result$FstSNP)  # Per-SNP Fst distribution

snpgdsClose(GDS)   # Always close the GDS connection when finished

# ------------------------------------------------------------------------------
# 5b. Pairwise Fst / Nei's distance among populations (StAMPP)
# ------------------------------------------------------------------------------

# Ensure Region is set as the active population
GBS@pop  <- samples$Region
pop(GBS) <- GBS@pop

# Nei's 1972 genetic distance between individuals and populations
D.ind <- stamppNeisD(GBS, pop = FALSE)   # individual-level distances
D.pop <- stamppNeisD(GBS, pop = TRUE)    # population-level distances

write.csv(D.pop, "FstRegion.csv")

# Ensure ploidy is set correctly (required by some StAMPP operations)
GBS@ploidy <- as.integer(ploidy(GBS))


################################################################################
# SECTION 6: HETEROZYGOSITY & DIVERSITY INDICES
################################################################################

# dartR compliance check populates required slots and validates the object
GBS <- gl.compliance.check(GBS)

# Recalculate locus metrics (fills @other$loc.metrics with Ho, He, etc.)
GBS <- gl.recalc.metrics(GBS)
GBS   # print summary

# Observed and expected heterozygosity by population and by individual
gl.report.heterozygosity(GBS, method = "pop")
gl.report.heterozygosity(GBS, method = "ind")

# General diversity indices (Hill numbers: D0, D1, D2)
# Run first at the Region level, then at the Site level
setPop(GBS)     <- ~Pop
diversity_stats <- gl.report.diversity(GBS)

diversity <- data.frame(
  pop   = rownames(diversity_stats$nlocpairpop),
  nloci = diversity_stats$nlocpop,
  Da0   = diversity_stats$zero_D_alpha,   # species richness equivalent
  Da1   = diversity_stats$one_D_alpha,    # Shannon-entropy equivalent
  Da2   = diversity_stats$two_D_alpha     # Simpson-diversity equivalent
)
write.csv(diversity, "diversity_region.csv", row.names = FALSE)

# Repeat at the sub-population (Site) level if needed
setPop(GBS) <- ~Subpop


################################################################################
# SECTION 7: OUTFLANK – OUTLIER LOCI (POTENTIAL SELECTION)
################################################################################

# Filter loci: remove those with high missingness and monomorphic loci
# threshold = 0.05 means loci present in at least 95% of individuals
GBS_clean <- gl.filter.callrate(GBS, method = "loc", threshold = 0.05)
GBS_clean <- gl.filter.monomorphs(GBS_clean)
GBS_clean
table(pop(GBS_clean))

# NOTE: OutFLANK requires at least 2 individuals per population.
# If any population has too few samples, merge them with a neighbouring group.
# Example: merge Tasmania ("TAS") into the broader Australia ("AUS") group.
pop(GBS_clean)[pop(GBS_clean) == "TAS"] <- "AUS"
table(pop(GBS_clean))

# Run OutFLANK
# LeftTrimFraction / RightTrimFraction: proportion of loci trimmed before
#   fitting the neutral Fst distribution (removes extreme low/high Fst loci).
# Hmin: minimum observed heterozygosity – loci below this are excluded.
# qthreshold: false discovery rate threshold for calling outliers.
outlier_results <- gl.outflank(
  gi                = GBS_clean,
  plot              = TRUE,
  LeftTrimFraction  = 0.03,
  RightTrimFraction = 0.05,
  Hmin              = 0.05,
  qthreshold        = 0.2
)

# Summarise results
outlier_results$outflank$FSTbar                  # Mean Fst of neutral loci
outlier_results$outflank$dfInferred              # Inferred d.f. for neutral model
outlier_results$outflank$numberHighFstOutliers   # Directional-selection candidates
outlier_results$outflank$numberLowFstOutliers    # Balancing-selection candidates

# Extract the outlier loci (flagged == TRUE, excluding NAs)
outliers <- outlier_results$outflank$results[
  !is.na(outlier_results$outflank$results$OutlierFlag) &
    outlier_results$outflank$results$OutlierFlag == TRUE, ]
outliers
outlier_names <- outliers$LocusName
outlier_names

write.csv(outlier_results$outflank$results, "outflank.csv")


################################################################################
# SECTION 8: NEIGHBOUR-JOINING TREE
################################################################################

# Reset plotting parameters
par(mar   = c(4, 3, 2, 3))
par(mfrow = c(1, 1))

# Build NJ tree from Euclidean distances on the allele-frequency matrix
tree <- njs(dist(as.matrix(GBS)))

# Colour palettes:
#   palette_tree: one colour per individual (recycled), for tip labels
#   palette_pop:  one colour per population, for tip point symbols
# NOTE: Change the recycling length (6) to match your number of populations.
palette_tree <- rep_len(
  c("#FF0000", "#FF6E00", "#FFC300", "#FFFF00", "#AAD500",
    "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082",
    "#812BA6", "#B857CA", "#D03A87"),
  6
)
palette_pop <- c("#FF0000", "#AAD500", "#0000FF",
                 "#812BA6", "#FFFF00", "#005555")

# Plot tree
plot(tree, "phylogram",
     cex             = 0.75,
     use.edge.length = FALSE,
     font            = 2,
     node.pos        = 1,
     edge.width      = 2,
     label.offset    = 0.5)
tiplabels(pch = 20, cex = 1, col = palette_pop[as.numeric(pop(GBS))])
axisPhylo()
title("Neighbour-joining tree of filtered data")

write.tree(tree, "ALLtree.txt")


################################################################################
# END OF SCRIPT
# Key objects for downstream scripts:
#   vcf, GBS, samples  →  used in 02_PCA_AMOVA.R, 03_structure_DAPC.R,
#                          and 04_mantel_partial_mantel.R
#   D.ind, D.pop       →  used in 02_PCA_AMOVA.R and 04_mantel_partial_mantel.R
################################################################################
