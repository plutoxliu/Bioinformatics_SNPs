################################################################################
##                                                                            ##
##        POPULATION GENOMICS — MANTEL & PARTIAL MANTEL TESTS                 ##
##                                                                            ##
##  Sections:                                                                 ##
##   1.  Mantel Test (Isolation by Distance)                                  ##
##   2.  Partial Mantel Test (IBD controlling for sampling density)           ##
##                                                                            ##
##  Requires objects from 01_data_QC_diversity.R:                             ##
##    vcf, GBS, samples, D.ind, D.pop, D.ind.dist, D.pop.dist                 ##
##                                                                            ##
################################################################################


# ==============================================================================
# LIBRARIES
# ==============================================================================

library(vcfR)
library(adegenet)
library(StAMPP)
library(vegan)
library(tidyverse)
library(stringr)
library(ggplot2)
library(MASS)    # kde2d for density plots
library(ks)      # Hpi / kde for kernel density estimation
library(dartR)   # gi2gl


################################################################################
# SECTION 13: MANTEL TEST (ISOLATION BY DISTANCE)
################################################################################

# Convert VCF to genind format (required by ade4 / poppr functions)
GBS_genind <- vcfR2genind(vcf)

pops      <- samples$Region
subpops   <- samples$Site
my_strata <- data.frame(populations = pops, subpopulations = subpops)
strata(GBS_genind)  <- my_strata
setPop(GBS_genind)  <- ~populations
GBS_genind@pop

# ------------------------------------------------------------------------------
# 13a. Build geographic distance matrix
# ------------------------------------------------------------------------------

# Load site coordinates file
# Expected columns: Site, Region, Lat, Long
coor <- read.csv("coor.csv", header = TRUE)
coor <- coor %>%
  filter(!duplicated(Site)) %>%
  dplyr::select(Region, 
                Site, 
                Lat, 
                Lon)

# Join coordinates onto the sample metadata table
samples$Site<-toupper(samples$Site)
coor$Site<-toupper(coor$Site)

pop_geo  <- left_join(samples, coor, by = "Site")

gidtable <- data.frame(subpops = pop_geo$Site, pops = pop_geo$Region)
gidtable$subpops <- str_to_upper(gidtable$subpops)
gidtable <- left_join(gidtable, coor, by = c("subpops" = "Site"))

# Optional: use regional centroid coordinates instead of site-level coordinates
region_mean <- coor %>%
  group_by(Region) %>%
  summarise(x_mean = mean(Long), y_mean = mean(Lat)) %>%
  ungroup()

gidtable <- left_join(gidtable, region_mean, by = "Region")

# Attach individual-level coordinates to the genind object (one row per individual)
GBS_genind@other$xy           <- cbind(x=gidtable$Lon, y=gidtable$Lat)
rownames(GBS_genind@other$xy) <- gidtable$subpops

#Site-level coordinate table (one row per site, in the same order as D.pop)
rownames(D.pop) <- str_to_upper(rownames(D.pop))
site_coords <- coor %>%
  mutate(Site = str_to_upper(Site)) %>%
  filter(Site %in% rownames(D.pop)) %>%          # keep only sites in your genetic matrix
  arrange(match(Site, rownames(D.pop))) %>%       # enforce same order as D.pop
  dplyr::select(Site, Lon, Lat)

# Confirm order matches genetic distance matrix
stopifnot(all(site_coords$Site == rownames(D.pop)))

# Geographic distance matrix (site × site)
xy_mat <- as.matrix(dplyr::select(site_coords, Lon, Lat))
rownames(xy_mat) <- site_coords$Site

Dgeo <- dist(xy_mat)

# ------------------------------------------------------------------------------
# 13b. Mantel test
# ------------------------------------------------------------------------------

IBD <- mantel.randtest(Dgeo, D.pop.dist, nrepet = 10000)
IBD

# Basic IBD scatter plot
plot(Dgeo, D.pop.dist, pch = 20, cex = 0.5)
abline(lm(D.pop.dist ~ Dgeo))

# ------------------------------------------------------------------------------
# 13c. KDE-density IBD plot
# Overlaying a 2D kernel density estimate highlights areas of high pair density
# ------------------------------------------------------------------------------

dens  <- kde2d(Dgeo, D.pop.dist, n = 300, lims = c(-1, 16, 0, 0.08))
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))

plot(Dgeo, D.pop.dist, pch = 20, cex = 0.5)
image(dens, col = transp(myPal(300), 0.7), add = TRUE)
abline(lm(D.pop.dist ~ Dgeo))
title("Correlation of Genetic and Geographical Distances")


################################################################################
# SECTION 14: PARTIAL MANTEL TEST (IBD CONTROLLING FOR SAMPLING DENSITY)
################################################################################

# This analysis tests whether the genetic–geographic distance correlation
# holds after accounting for local sampling density (estimated via KDE).
# Denser sampling areas could create spurious IBD signals if unaccounted for.

# ------------------------------------------------------------------------------
# 14a. Attach individual-level coordinates to genind and convert to genlight
# ------------------------------------------------------------------------------

GBS_genind@other$xy           <- cbind(x=pop_geo$Lon, y=pop_geo$Lat)
rownames(GBS_genind@other$xy) <- pop_geo$Site
pop_geo$Region<-as.factor(pop_geo$Region)
pop_geo$Site<-as.factor(pop_geo$Site)
GBS@pop  <- pop_geo$Site
pop(GBS) <- GBS@pop
pop(GBS)

GBS_gl <- gi2gl(GBS_genind)   # convert genind → genlight

# ------------------------------------------------------------------------------
# 14b. Build population-level coordinate table
# Average individual coordinates within each population
# ------------------------------------------------------------------------------

gidtable2 <- data.frame(GBS_gl@other$xy,
                         pop = row.names(GBS_gl@other$xy))
gidtable2 <- gidtable2 %>%
  summarise(x = mean(x), y = mean(y), .by = pop)

# Coordinate matrices for individuals and populations
ind_den <- as.matrix(GBS_genind@other$xy)
pop_den <- as.matrix(gidtable2[, c("x", "y")])
rownames(pop_den) <- gidtable2$pop

# ------------------------------------------------------------------------------
# 14c. Kernel density estimation (KDE)
# Hpi()     – selects the optimal bandwidth matrix via plug-in estimation
# kde()     – computes the bivariate density surface over the coordinate space
# predict() – evaluates the density at each sampled location
# ------------------------------------------------------------------------------

H_ind <- Hpi(x = ind_den)
H_pop <- Hpi(x = pop_den)

fhat_ind <- kde(x = ind_den, H = H_ind)
fhat_pop <- kde(x = pop_den, H = H_pop)

individual_densities <- predict(fhat_ind, x = ind_den)
pop_densities        <- predict(fhat_pop, x = pop_den)

names(individual_densities) <- rownames(ind_den)
names(pop_densities)        <- rownames(pop_den)

# Convert density values to dissimilarity distance matrices
Dkde_ind <- dist(individual_densities)
Dkde_pop <- dist(pop_densities)
Dgeo_ind <- dist(ind_den)
Dgeo_pop <- dist(pop_den)

# ------------------------------------------------------------------------------
# 14d. Recalculate Nei's distances at the site level for partial Mantel
# (Re-run here because population assignment was changed to Site above)
# ------------------------------------------------------------------------------

D.ind <- stamppNeisD(GBS, pop = FALSE)
D.pop <- stamppNeisD(GBS, pop = TRUE)

colnames(D.ind) <- rownames(D.ind)
D.ind.dist      <- as.dist(D.ind, diag = TRUE)
attr(D.ind.dist, "Labels") <- rownames(D.ind)

colnames(D.pop) <- rownames(D.pop)
D.pop.dist      <- as.dist(D.pop, diag = TRUE)
attr(D.pop.dist, "Labels") <- rownames(D.pop)

# Quick alignment check before running tests
head(D.pop)
head(pop_den)
head(pop_densities)

# ------------------------------------------------------------------------------
# 14e. Partial Mantel tests
# Tests: Genetic distance ~ Geographic distance | KDE density
# Increase permutations to 9999 for publication-quality p-values
# ------------------------------------------------------------------------------

# Population level
mantel.partial(D.pop.dist, Dgeo_pop, Dkde_pop, permutations = 1000)

# Individual level
mantel.partial(D.ind.dist, Dgeo_ind, Dkde_ind, permutations = 1000)

# ------------------------------------------------------------------------------
# 14f. IBD scatter plot with KDE density contours
# NOTE: Update the subtitle with your actual Mantel r and p-value after running
# ------------------------------------------------------------------------------

df_ibd <- data.frame(
  GenDist = as.vector(D.pop.dist),
  GeoDist = as.vector(Dgeo_pop)
)

p_ibd <- ggplot(df_ibd, aes(x = GeoDist, y = GenDist)) +
  # KDE density fill – shows where most population pairs fall
  stat_density_2d(aes(fill = after_stat(level)),
                  geom = "polygon", alpha = 0.3) +
  # Individual pair points
  geom_point(alpha = 0.2, size = 0.5, color = "slategray") +
  # Mantel regression line
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(
    title    = "IBD Scatterplot with KDE Density Contours (by site)",
    subtitle = "Partial Mantel r = 0.249 (p = 0.170)",   # update after running
    x        = "Geographic Distance",
    y        = "Genetic Distance",
    fill     = "Pair Density"
  )

p_ibd
ggsave(plot = p_ibd, "IBD_KDE_site.png",
       width = 8, height = 6, dpi = 300)


################################################################################
# END OF SCRIPT
################################################################################
