################################################################################
##                                                                            ##
##          POPULATION GENOMICS — STRUCTURE ANALYSIS & DAPC                   ##
##                                                                            ##
##  Sections:                                                                 ##
##   1.  Structure Analysis (LEA / sNMF)                                      ##
##   2.  DAPC (Discriminant Analysis of Principal Components)                 ##
##                                                                            ##
##  Requires objects from 01_data_QC_diversity.R:                             ##
##    vcf, GBS, samples                                                       ##
##                                                                            ##
################################################################################


# ==============================================================================
# LIBRARIES
# ==============================================================================

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)
library(stringr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)


# Shared colour palette for populations — used across both sections
palette_pop2 <- c("#e31a1c", "#33a02c", "#0072B2", "#fb9a99", "#F0E442","#6a3d9a","#b2df8a", "#56B4E9", "#ff7f00",  "#a6cee3")


################################################################################
# SECTION 1: STRUCTURE ANALYSIS (LEA / sNMF)
################################################################################

# LEA's sNMF requires a .geno file where genotypes are encoded as:
#   0 = homozygous reference
#   1 = heterozygous
#   2 = homozygous alternate
#   9 = missing data

# ------------------------------------------------------------------------------
# 10a. Encode genotypes and write .geno file
# ------------------------------------------------------------------------------

tidy_geno <- vcf %>%
  extract_gt_tidy() %>%
  # Keep only the genotype column; drop all other FORMAT fields
  select(-gt_DP, -gt_GT_alleles, -gt_AD, -gt_HQ, -gt_GL, -gt_GQ, -gt_MY) %>%
  mutate(
    gt1 = str_split_fixed(gt_GT, "/", n = 2)[, 1],
    gt2 = str_split_fixed(gt_GT, "/", n = 2)[, 2],
    geno_code = case_when(
      gt1 == 0 & gt2 == 0      ~ 0,   # homozygous reference
      gt1 == 0 & gt2 == 1      ~ 1,   # heterozygous
      gt1 == 1 & gt2 == 0      ~ 1,
      gt1 == 1 & gt2 == 1      ~ 2,   # homozygous alternate
      gt1 == "" | gt2 == ""    ~ 9    # missing
    )
  ) %>%
  select(-gt_GT, -gt1, -gt2)

# Pivot so individuals are columns and SNPs are rows (LEA's expected format)
geno <- tidy_geno %>%
  pivot_wider(names_from = Indiv, values_from = geno_code) %>%
  select(-Key)

write.table(geno, "SARGgeno.geno",
            col.names = FALSE, row.names = FALSE, sep = "")

# ------------------------------------------------------------------------------
# 10b. Run sNMF for K = 1:10
# repetitions = 10 runs per K ensures results are stable across runs.
# alpha = 1000 is appropriate for smaller datasets (< 10,000 SNPs).
# ------------------------------------------------------------------------------

vcf_snmf <- snmf(
  input.file  = "SARGgeno.geno",
  K           = 1:10,
  entropy     = TRUE,
  repetitions = 20,
  project     = "new",
  alpha       = 1000
)

# Cross-entropy plot: the best K is at the elbow (minimum or clear inflection)
plot(vcf_snmf, cex = 1.2, col = "lightblue", pch = 19)

# Select the run with the lowest cross-entropy for the chosen K
# NOTE: Adjust chosen_K based on the cross-entropy plot above
chosen_K <- 3
ce       <- cross.entropy(vcf_snmf, K = chosen_K)
best_run <- which.min(ce)
best_run

# ------------------------------------------------------------------------------
# 10c. Extract Q-matrix and prepare for plotting
# Each row is an individual; each column is an ancestral population (cluster)
# ------------------------------------------------------------------------------
pop<-samples
q_mat           <- LEA::Q(vcf_snmf, K = chosen_K, run = best_run)
colnames(q_mat) <- paste0("P", seq_len(chosen_K))

q_df <- q_mat %>%
  as_tibble() %>%
  mutate(
    individual = pop$Individual,
    region     = pop$Region,
    order      = pop$Region   # controls plot ordering
  )

# Convert to long format for ggplot
q_df_long <- q_df %>%
  pivot_longer(cols      = starts_with("P"),
               names_to  = "pop",
               values_to = "q")

# Arrange individuals by region for a geographically ordered plot
q_df_prates <- q_df_long %>%
  arrange(order) %>%
  mutate(individual = forcats::fct_inorder(factor(individual)))

# ------------------------------------------------------------------------------
# 10d. Pie charts – mean admixture proportion per region
# ------------------------------------------------------------------------------

q_df_pie <- q_df_prates %>%
  group_by(region, pop) %>%
  summarise(mean = mean(q), .groups = "drop")

ggplot(q_df_pie,
       aes(x = " ", y = mean, group = region, color = pop, fill = pop)) +
  geom_bar(stat = "identity", color = NA) +
  coord_polar("y", start = 0) +
  facet_grid(. ~ region) +
  scale_fill_manual(values = palette_pop2) +
  theme_void()

ggsave("pies.pdf", device = "pdf", width = 12, height = 6, bg = "transparent")
# ------------------------------------------------------------------------------
# 10e. Admixture bar plot (STRUCTURE-style)
# ------------------------------------------------------------------------------

p_admix <- q_df_prates %>%
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  facet_grid(cols = vars(region), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_pop2) +
  labs(fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.spacing.x  = unit(0, "lines"),
    axis.line        = element_blank(),
    axis.text        = element_blank(),
    strip.background = element_rect(fill = "transparent", color = "black"),
    panel.background = element_blank(),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 90, hjust = 1, size = 8)
  )
p_admix


################################################################################
# SECTION 2: DAPC (DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS)
################################################################################

# ------------------------------------------------------------------------------
# 11a. K-means BIC plot to identify the optimal number of clusters
# 10 replicate runs account for stochastic variation in k-means initialisation
# ------------------------------------------------------------------------------

gl_data <- vcfR2genlight(vcf)   # fresh genlight (in case GBS has been modified)

maxK  <- 10
myMat <- matrix(nrow = 10, ncol = maxK)
colnames(myMat) <- 1:maxK

for (i in 1:nrow(myMat)) {
  grp        <- find.clusters(gl_data, n.pca = 20,
                               choose.n.clust = FALSE, max.n.clust = maxK)
  myMat[i, ] <- grp$Kstat   # BIC value for each K
}

bic_df          <- melt(myMat)
colnames(bic_df)[1:3] <- c("Group", "K", "BIC")
bic_df$K        <- as.factor(bic_df$K)

# BIC plot: the optimal K is typically at the elbow (minimum or clear inflection)
p1 <- ggplot(bic_df, aes(x = K, y = BIC)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Number of groups (K)")
p1

# ------------------------------------------------------------------------------
# 11b. Run DAPC for the selected range of K values
# n.pca = 140 retains many PCs; adjust based on your dataset size.
# set.seed() ensures reproducible cluster assignments across runs.
# ------------------------------------------------------------------------------

my_k <- c(2:3)   # range of K values to explore (minimum 2)

grp_l  <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for (i in seq_along(dapc_l)) {
  set.seed(9)
  grp_l[[i]]  <- find.clusters(gl_data, n.pca = 140, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(gl_data,
                      pop   = grp_l[[i]]$grp,
                      n.pca = 140,
                      n.da  = my_k[i])
}

# ------------------------------------------------------------------------------
# 11c. Scatter plot of LD1 vs LD2 for the highest K run
# ------------------------------------------------------------------------------

dapc_coords       <- as.data.frame(dapc_l[[length(dapc_l)]]$ind.coord)
dapc_coords$Group <- dapc_l[[length(dapc_l)]]$grp

my_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

p2 <- ggplot(dapc_coords, aes(x = LD1, y = LD2,
                               color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  scale_color_manual(values = my_pal) +
  scale_fill_manual(values  = paste0(my_pal, "66"))   # semi-transparent fill
p2

# ------------------------------------------------------------------------------
# 11d. Posterior probability bar plots across K values
# ------------------------------------------------------------------------------

# Collect posterior probabilities for all K values into one data frame
posterior_df <- NULL
for (i in seq_along(dapc_l)) {
  tmp             <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K           <- my_k[i]
  tmp$Isolate     <- rownames(tmp)
  tmp             <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region      <- pop$samples.site   # assign geographic label
  posterior_df    <- rbind(posterior_df, tmp)
}

grp.labs        <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(posterior_df, aes(x = Isolate, y = Posterior, fill = Group)) +
  geom_bar(stat = "identity") +
  facet_grid(K ~ Region,
             scales   = "free_x",
             space    = "free",
             labeller = labeller(K = grp.labs)) +
  theme_bw() +
  ylab("Posterior membership probability") +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 90, hjust = 1, size = 8)
  ) +
  scale_fill_manual(values = palette_pop2)
p3

# ------------------------------------------------------------------------------
# 11e. Combined multi-panel figure (BIC + scatter + posterior)
# ------------------------------------------------------------------------------

# tiff("dapc_k3_5_dapc.tiff", width = 6.5, height = 6.5,
#      units = "in", compression = "lzw", res = 300)
ggarrange(
  ggarrange(p1, p2, ncol = 2, labels = c("A", "B")),
  p3,
  nrow    = 2,
  labels  = c("", "C"),
  heights = c(1, 2)
)


################################################################################
# END OF SCRIPT
################################################################################
