# 10/1/2024
# Aaron Tayal
# I would like to import the carabid count data from Powdermill 2022 
# properly

library(ggplot2)
library(dplyr)

# Import the carabid counts into the dataframe "carab". I would like to keep
# all column classes as character, because in the csv, blank cells in a species 
# column indicate true zeros, whereas "NA" indicates the trap contents were
# lost or destroyed by raccoons.
carab <- read.csv("Aaron_PNR_formatted_csvs/PNR2022_carabid_counts.csv",
                  na.strings = "NA", colClasses = "character")

sum(is.na(carab)) # Looks like we have 336 NA values. These are due to 
# 7 rows * 48 columns where the NA values occur.
carab[is.na(carab$Agonum_fidele), 1:3] # This lists the samples that are missing
# for various reasons (raccoons, mislabelling, misplacing vial)

carab_species <- colnames(carab)[which(colnames(carab)=="Agonum_fidele"):dim(carab)[2]] 
# This creates a vector containing the names of the carabid species.

# Now I want to replace blank cells with zeros for the species count columns.
# To do this, iterate over species and over rows, and only write a zero in 
# the column if there is currently a blank value "" and not if there is NA
for (species in carab_species) {
  for (row in 1:nrow(carab)){
    if (carab[row, species] == "" && !is.na(carab[row, species])){
      carab[row, species] <- "0"
    }
  }
}

# Now change the column classes to integers:
for (species in carab_species) {
  carab[,species] <- as.integer(carab[,species])
}
carab$Total_Carabidae_during_pinning <- as.integer(carab$Total_Carabidae_during_pinning)
carab$Total_Carabidae_from_sum_of_species_counts <- as.integer(carab$Total_Carabidae_from_sum_of_species_counts)

# Change the columns that indicate collection interval and trap location into 
# factors
carab$Collection_interval <- as.factor(carab$Collection_interval)
carab$PNR_Code <- as.factor(carab$PNR_Code)

# Graph the number of carabids captured, with trap location in the x-axis:
ggplot(data=carab, mapping=aes(x=PNR_Code, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()
# Now with collection interval on the x-axis:
ggplot(data=carab, mapping=aes(x=Collection_interval, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()

# Import information about the trap locations, such as transect and area:
trap_locations <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("integer", "integer", "factor", "factor", "factor", "numeric", "numeric"))

#Remove rows with missing carabid count#########################################
# Note: 7 rows have missing carabid count
carab_no_missing <- carab %>% filter(!is.na(carab[,carab_species[1]]))

#taxonomic diversity metrics###################################################


str(carab)
carab$Treatment <- as.factor(carab$Treatment)

library(hillR)

## total abundance
carab$abund <- rowSums(carab[,15:60])
str(carab)
hist(carab$abund)
boxplot(abund ~ Treatment, data = carab,
        xlab = "", ylab = "Ground beetle Abundance")
stripchart(abund ~ Treatment, data = carab, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)


## Hill Numbers
# q = 0 (default) to get species richness, q = 1 to get shannon entropy,
# q = 2 will give inverse Simpson.
# margin 1 is the default, indicates that sites are rows

# species richness
carab$rich <- hill_taxa(carab[,15:60], q = 0, MARGIN = 1)
str(carab)
hist(carab$rich)
boxplot(rich ~ Treatment, data = carab,
        xlab = "", ylab = "Ground beetle Richness")
stripchart(rich ~ Treatment, data = carab, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# Shannon diversity
carab$sha.div <- hill_taxa(carab[,15:60], q = 1, MARGIN = 1)
str(carab)
hist(carab$sha.div)
boxplot(sha.div ~ Treatment, data = carab,
        xlab = "", ylab = "Ground beetle Shannon Diversity")
stripchart(sha.div ~ Treatment, data = carab, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# Simpson diversity
carab$sim.div <- hill_taxa(carab[,15:60], q = 2, MARGIN = 1)
str(carab)
hist(carab$sim.div)
boxplot(sim.div ~ Treatment, data = carab,
        xlab = "", ylab = "Ground beetle Simpson Diversity")
stripchart(sim.div ~ Treatment, data = carab, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

# evenness
library(chemodiv)
#  For q = 2, more weight is put on compounds with high proportions
carab$eve <- calcDiv(carab[,15:60], type = "HillEven", q = 2)
str(carab)
hist(carab$eve$HillEven)
boxplot(eve$HillEven ~ Treatment, data = carab,
        xlab = "", ylab = "Ground beetle Evenness")
stripchart(eve$HillEven ~ Treatment, data = carab, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

#Estimate species richness with accumulation curves#####################################################################
## 
# individual-based rarefaction by treatment, jackknife estimates by treatment
# we are separating out each treatment so we can get a species accumulation curve for each, and then we will graph it

# Note (Aaron): I substituted the dataframe where missing carabid rows are removed:

levels(carab_no_missing$Treatment)
TF <- carab_no_missing[which(carab_no_missing$Treatment == "F"),]
TS <- carab_no_missing[which(carab_no_missing$Treatment == "S"),]
TW <- carab_no_missing[which(carab_no_missing$Treatment == "W"),]

# the rarefaction method standardizes the sample sizes so that we are comparing species richness
# at equivalent abundances
library(vegan)

sp.TF <- specaccum(TF[,15:60], method = "rarefaction", permutations = 100, gamma = "jack2") 
sp.TS <- specaccum(TS[,15:60], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.TW <- specaccum(TW[,15:60], method = "rarefaction", permutations = 100, gamma = "jack2")


# make the plot

png("Rarefaction.png", width = 1200, height = 1000, pointsize = 30)

par(mfrow=c(1,1))
par(mar=c(5,6,4,2))

plot(sp.TS, pch = 19, col = "brown4", xvar = c("individuals"), lty = 4, lwd = 5,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 300), ylim = c(0, 40))
plot(sp.TF, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 5, col = "palegreen4")
plot(sp.TW, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 5, col = "goldenrod2")

legend("bottomright", legend = c("Salvaged", "Forest", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1.2, bty = "n", lwd = 5,
       col = c("brown4", "palegreen3", "goldenrod2"))

dev.off()


## calculates the estimated species richness of a population using first- and second-order jackknife estimators
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons
library(fossil)

specnumber(carab[,15:60], groups = carab$Treatment)

# forest
jack1(TF[,15:60], taxa.row = FALSE, abund = TRUE)
jack2(TF[,15:60], taxa.row = FALSE, abund = TRUE)

# salvaged
jack1(TS[,15:60], taxa.row = FALSE, abund = TRUE)
jack2(TS[,15:60], taxa.row = FALSE, abund = TRUE)

# windthrow
jack1(TW[,15:60], taxa.row = FALSE, abund = TRUE)
jack2(TW[,15:60], taxa.row = FALSE, abund = TRUE)


#Nonmetric multidimensional scaling (NMDS)#####################################
# compute a dissimilarity matrix
# the method option let's you indicate which dissimilarity metric to calculate
# will calculate the a bray-curtis dissimilarity matrix for abundance-based data
dis.matrix <- vegdist(carab[,15:60], method = "bray")
dis.matrix

# run the nonmetric multidimensional scaling model
nmds.carabid <- metaMDS(dis.matrix, trymax = 500, autotransform = TRUE, k = 2)
nmds.carabid # stress is quality of fit
stressplot(nmds.carabid)
plot(nmds.carabid) # basic plot with no treatment distinctions

# plot the NMDS model by treatment
ordiplot(nmds.carabid, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-1, 1))
points(nmds.carabid, dis = "sites", select = which(carab$Treatment=="W"), pch = 15, cex = 2, col = "goldenrod2")
points(nmds.carabid, dis = "sites", select = which(carab$Treatment=="F"), pch = 16, cex = 2, col = "palegreen4")
points(nmds.carabid, dis = "sites", select = which(carab$Treatment=="S"), pch = 16, cex = 2, col = "brown4")

ordiellipse(nmds.carabid, carab$Treatment, draw = "lines", col = c("goldenrod2", "palegreen4", "brown4"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("bottomright", legend = c("Windthrow", "Forest", Salvaged),
       pch = c(15, 16), cex = 1.5, bty = "n", col = c("goldenrod2", "palegreen4", "brown4"))

## Test for differences in ant composition among treatments
# PERMANOVA tests whether the group centroid of communities differs among groups
# in multivariate space (e.g. different community composition)
adonis2(dis.matrix ~ carab$Treatment, permutations = 999)

# BETADISPER tests whether the dispersion of a group from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogeneity of variances
c.beta <- betadisper(dis.matrix, carab$Treatment, type = c("median"))
anova(c.beta)
plot(c.beta)
boxplot(c.beta, ylab = "Distance to median")
TukeyHSD(c.beta, which = "group", conf.level = 0.95)

