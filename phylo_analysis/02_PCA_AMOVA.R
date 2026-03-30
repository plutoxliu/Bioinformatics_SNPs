################################################################################
##                                                                            ##
##                POPULATION GENOMICS — PCA & AMOVA                           ##
##                                                                            ##
##  Sections:                                                                 ##
##    1.  Principal Components Analysis (PCA)                                 ##
##    2.  AMOVA (Analysis of Molecular Variance)                              ##
##                                                                            ##
##  Requires objects from 01_data_QC_diversity.R:                             ##
##    vcf, GBS, samples, D.ind, D.pop                                         ##
##                                                                            ##
################################################################################


# ==============================================================================
# LIBRARIES
# ==============================================================================

library(adegenet)
library(StAMPP)
library(ggplot2)
library(pegas)
library(gplots)    # for heatmap.2


################################################################################
# SECTION 1: PRINCIPAL COMPONENTS ANALYSIS (PCA)
################################################################################

# ------------------------------------------------------------------------------
# 9a. Eigenvalue exploration
# Run a first PCA with more PCs to decide how many to retain
# ------------------------------------------------------------------------------

data_pca <- glPca(GBS, nf = 3)

# NOTE: If "NAs detected" error appears, uncomment and run the block below:
# toRemove <- is.na(glMean(GBS, alleleAsUnit = FALSE))
# GBS2     <- GBS[, !toRemove]
# data_pca <- glPca(GBS2, nf = 4)

# Eigenvalues and percentage variance explained per PC
eig.val  <- data_pca$eig
eig.perc <- 100 * data_pca$eig / sum(data_pca$eig)
eigen    <- data.frame(eig.val, eig.perc)
head(eigen)
# write.csv(eigen, "eigen-summary.csv", row.names = TRUE, quote = FALSE)

# Scree plot – use this to decide how many PCs to retain below
dev.off()
plot.new()
par(mar = c(3, 3, 3, 3))
barplot(data_pca$eig, main = "Eigenvalues",
        col = heat.colors(length(data_pca$eig)))

# ------------------------------------------------------------------------------
# 9b. PCA scores and plotting
# ------------------------------------------------------------------------------

# Run final PCA retaining 3 PCs (adjust nf as needed based on scree plot)
pca1 <- glPca(GBS, nf = 3)

# Basic scatter plot of all individuals
scatter(pca1, ratio = 0.2)

# Colour palette for populations
palette_pop2 <- c("#e31a1c", "#33a02c", "#0072B2", "#fb9a99", "#F0E442","#6a3d9a","#b2df8a", "#56B4E9", "#ff7f00",  "#a6cee3")

# Population-level scatter with confidence ellipses
s.class(pca1$scores, fac = GBS$pop, col = palette_pop2)
text(pca1$scores, labels = indNames(GBS), cex = 0.5)
add.scatter.eig(pca1$eig, 3, 2, 1,
                pos = "bottomleft", inset = 0.01, ratio = 0.18)

# Extract scores into a data frame for ggplot
pca_scores     <- as.data.frame(pca1$scores)
pca_scores$pop <- pop(GBS)

# NOTE: palette length (146 here) must match the number of individuals.
palette_ind <- rep_len(
  c("#FF0000", "#FF6E00", "#FFC300", "#FFFF00", "#AAD500",
    "#008000", "#005555", "#0000FF", "#3200AC", "#4B0082",
    "#812BA6", "#B857CA", "#D03A87"),
  286
)

# PC1 vs PC2
ggplot(pca_scores, aes(PC1, PC2)) +
  stat_ellipse(aes(colour = pop)) +
  geom_point(aes(colour = pop), size = 3, alpha = 0.7, shape = 19) +
  scale_color_manual(values = palette_pop2) +
  theme_classic()

# PC2 vs PC3
ggplot(pca_scores, aes(PC2, PC3)) +
  geom_point(aes(colour = pop), size = 3, alpha = 0.7, shape = 19) +
  stat_ellipse(aes(colour = pop)) +
  scale_color_manual(values = palette_pop2) +
  theme_classic()

# PC1 vs PC3
ggplot(pca_scores, aes(PC1, PC3)) +
  geom_point(aes(colour = pop), size = 3, alpha = 0.7, shape = 19) +
  stat_ellipse(aes(colour = pop)) +
  scale_color_manual(values = palette_pop2) +
  theme_classic()


################################################################################
# SECTION 2: AMOVA (ANALYSIS OF MOLECULAR VARIANCE)
################################################################################

# Reset population to Region level
GBS@pop  <- samples$Region
pop(GBS) <- GBS@pop

# Visual check of genotype matrix (0/1/2; white = missing) — takes a moment
glPlot(GBS)

# Count missing SNPs per sample and export
# Row 7 of summary() contains the NA counts
x <- summary(t(as.matrix(GBS)))
write.table(x[7, ], file = "missing.persample.txt", sep = "\t")

# Recalculate Nei's distance matrices at Region level
D.ind      <- stamppNeisD(GBS, pop = FALSE)
D.pop      <- stamppNeisD(GBS, pop = TRUE)
GBS@ploidy <- as.integer(ploidy(GBS))

# ------------------------------------------------------------------------------
# Heatmaps of distance matrices
# ------------------------------------------------------------------------------

colnames(D.ind) <- rownames(D.ind)
pdf("Neis_dist_heatmap.pdf", width = 10, height = 10)
heatmap.2(D.ind, trace = "none", cexRow = 0.4, cexCol = 0.4)
dev.off()

colnames(D.pop) <- rownames(D.pop)
pdf("Neis_dist_heatmap_pop.pdf", width = 10, height = 10)
heatmap.2(D.pop, trace = "none", cexRow = 0.4, cexCol = 0.4)
dev.off()

# ------------------------------------------------------------------------------
# Convert to dist objects (required by pegas::amova)
# ------------------------------------------------------------------------------

colnames(D.ind) <- rownames(D.ind)
D.ind.dist      <- as.dist(D.ind, diag = TRUE)
attr(D.ind.dist, "Labels") <- rownames(D.ind)

colnames(D.pop) <- rownames(D.pop)
D.pop.dist      <- as.dist(D.pop, diag = TRUE)
attr(D.pop.dist, "Labels") <- rownames(D.pop)

# ------------------------------------------------------------------------------
# AMOVA
# Grouping factors must be factors derived from the samples table.
# NOTE: column names in 'samples' are case-sensitive — use "Site" not "site".
# ------------------------------------------------------------------------------

sites   <- as.factor(samples$Site)     # fine-scale grouping
regions <- as.factor(samples$Region)   # broad grouping

rownames(samples) <- samples$Individual

# One-level AMOVA: variation partitioned among regions
(res_one <- pegas::amova(D.ind.dist ~ regions, nperm = 10000))

# Hierarchical AMOVA: regions as top level, sites nested within regions
(res_hier <- pegas::amova(D.ind.dist ~ regions / sites, nperm = 10000))


################################################################################
# END OF SCRIPT
################################################################################
