# 10/1/2024
# Aaron Tayal
# Carabid beetles from Powdermill Nature Reserve 2022, in wind-disturbed,
# salvage-logged, and undisturbed forest plots.

library(ggplot2)
library(dplyr)

# Import data #################################################################

# Import the carabid counts into the dataframe "carab". I would like to keep
# all column classes as character, because in the csv, blank cells in a species 
# column indicate true zeros, whereas "NA" indicates the trap contents were
# lost or destroyed by raccoons.
carab0 <- read.csv("Aaron_PNR_formatted_data/PNR2022_carabid_counts.csv",
                  na.strings = "NA", colClasses = "character")

sum(is.na(carab0)) # Looks like we have 336 NA values. These are due to 
# 7 rows * 48 columns where the NA values occur.
carab0[is.na(carab0$Agonoleptus_thoracicus), 3:4] # This lists the samples that are missing
# for various reasons (raccoons, mislabelling, misplacing vial)

carab_species <- colnames(carab0)[which(colnames(carab0)=="Agonoleptus_thoracicus"):dim(carab0)[2]] 
# This creates a vector containing the names of the carabid species.
# Note: In the datasheet, Agonoleptus thoracicus MUST BE THE FIRST species column

# Now I want to replace blank cells with zeros for the species count columns.
# To do this, iterate over species and over rows, and only write a zero in 
# the column if there is currently a blank value "" and not if there is NA
for (species in carab_species) {
  for (row in 1:nrow(carab0)){
    if (carab0[row, species] == "" && !is.na(carab0[row, species])){
      carab0[row, species] <- "0"
    }
  }
}

# Now change the column classes to integers:
for (species in carab_species) {
  carab0[,species] <- as.integer(carab0[,species])
}
carab0$Total_Carabidae_during_pinning <- as.integer(carab0$Total_Carabidae_during_pinning)
carab0$Total_Carabidae_from_sum_of_species_counts <- as.integer(carab0$Total_Carabidae_from_sum_of_species_counts)

# Change the columns that indicate collection interval and trap location into 
# factors
carab0$Collection_interval <- as.factor(carab0$Collection_interval)
carab0$Plot <- as.factor(carab0$Plot)
carab0$Set_date <- as.factor(carab0$Set_date)
carab0$Collection_date <- as.factor(carab0$Collection_date)

carab0 <- carab0 %>% select(-Treatment, -PNR_Code) # get rid of unnecessary columns

# Graph the number of carabids captured, with trap location in the x-axis:
ggplot(data=carab0, mapping=aes(x=Plot, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()
# Now with collection interval on the x-axis:
ggplot(data=carab0, mapping=aes(x=Collection_interval, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_jitter(height=0, width=0.1, alpha=0.5) +
  theme_classic()
# Graph a histogram of the number of carabids captured:
ggplot(data=carab0, mapping=aes(x=Total_Carabidae_from_sum_of_species_counts)) +
  geom_histogram(breaks=seq(-0.5,23.5,1)) + theme_classic()

# Import information about the trap locations, such as transect and area:
trap_locations <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv")
library(forcats)
trap_locations <- mutate(trap_locations, Treatment = fct_recode(Treatment, 
                            "Forest" = "F", "Windthrow" = "W", "Salvaged" = "S"))
trap_locations$Plot <- as.factor(trap_locations$Plot)

# Merge the carab dataframe with the trap_locations dataframe:
carab1 <- right_join(trap_locations, carab0, by="Plot") %>% arrange(Interval)
table(carab1$Plot)
table(carab1$Interval)

# Look at which treatments have missing data:
carab1[is.na(carab1$Agonum_fidele), c("Interval", "Plot", "Treatment")] 
# The missing rows are all from the undisturbed forest treatment

# Deal with missing data and dumped traps ######################################

# Adjust carabid counts based on dumped traps. Each pitfall has two collection 
# cups, and sometimes one or both of those cups was dumped. 2/2 cups were dumped
# four times (resulting in missing data). 1/2 cups were dumped 12 times (1 of
# those 12 is also missing data due to losing the vial):
carab1[carab1$Traps_Dumped_copied_from_PNR_ENV_2022 == 1, c("Plot", 
                                                          "Collection_interval",
                                                          "Total_Carabidae_from_sum_of_species_counts",
                                                          "Treatment")]
# As you can see, more traps got dumped in the earlier collection intervals.
# To try and adjust the carabid counts for the dumped cups, I will multiply
# by 2 the counts of those rows where 1 out of 2 cups were dumped.
carab2 <- carab1
carab2[carab2$Traps_Dumped_copied_from_PNR_ENV_2022 == 1, 
           c("Total_Carabidae_from_sum_of_species_counts", carab_species)] <- 
  carab1[carab1$Traps_Dumped_copied_from_PNR_ENV_2022 == 1, 
            c("Total_Carabidae_from_sum_of_species_counts", carab_species)] * 2
# Now, the dataframe carab2 has adjusted carabid counts for those samples where
# 1 of the 2 cups was dumped:
carab2$Total_Carabidae_from_sum_of_species_counts - carab1$Total_Carabidae_from_sum_of_species_counts

ggplot(data=carab2, mapping=aes(x=Collection_interval, y=Total_Carabidae_from_sum_of_species_counts,
                                color=Traps_Dumped_copied_from_PNR_ENV_2022)) +
  geom_jitter(height=0, width=0.1, alpha=0.5) + theme_classic() +
  theme(legend.title=element_blank()) # Color dots if the abundance was corrected

# Remove rows with missing carabid count (Some analyses require these rows to be 
# removed). Note: 7 rows have missing carabid count
carab2_no_missing <- carab2 %>% filter(!is.na(carab2[,carab_species[1]]))

# Pool across sampling interval ################################################

# The trap intervals are:
levels(carab2$Set_Date)
levels(carab2$Collection_Date)
# Each trap was collected every 14-15 days.
# Look at missing data:
#View(carab0[is.na(carab0$Agonum_fidele), 3:4]) 
# Collection interval 0 has 13 or 14 days
# Collection interval 1 has 14 days
# Collection interval 4 has 15 days

# When were traps initially set up:
#View(carab2[carab2$Collection_interval==0, c("Plot", "Set_Date")]) 
# When were they collected:
#View(carab2[carab2$Collection_interval==7, c("Plot", "Collection_Date")])
# Sept 20th - June 1st = 30+31+31+20 = 112
# Sept 20th - June 2nd = 111
# Write down how long each trap was out:
# THIS SEEMS TO BE INCORRECT!
days_active <- data.frame(Plot=c( 41,  42,  43,  44,  45,  46,  47,     48,  49,  50,  51,     52,  53,     54,  55,  56,  57,  58,  59,  60,  61,        62,  63,  64), 
                   days_active=c(112, 111, 111, 112, 111, 111, 111, 111-13, 111, 111, 111, 111-13, 112, 112-14, 112, 112, 112, 112, 112, 112, 112, 111-13-14, 112, 112))
days_active$Plot <- as.factor(days_active$Plot)

# Use dplyr::summarize to create a new dataframe with a row for every Plot. The
# columns will be the total number of each species captured, across all 
# collection intervals.

carab2_by_plot_init <- carab2 %>% group_by(Plot) %>% 
  summarise(num_nonmissing_intervals = n() - sum(is.na(Agonoleptus_thoracicus)),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

carab2_by_plot_adding_plot_info <- right_join(trap_locations, 
                                              carab2_by_plot_init, by="Plot")

carab2_by_plot_adding_days_active <- right_join(days_active, carab2_by_plot_adding_plot_info, 
                             by="Plot")

carab2_by_plot <- carab2_by_plot_adding_days_active

days_trap_operational <- carab2_by_plot[,c("Plot", "Treatment", "days_active")]
#write.csv(days_trap_operational)

# The below code divides abundances of each species at each trap by the days
# active. 
carab2_by_plot_beetles_per_day <- carab2_by_plot # make a new data table for the 
# standardized values
for (species in carab_species) {
  carab2_by_plot_beetles_per_day[,species] <- 
    carab2_by_plot_adding_days_active[,species] / 
    carab2_by_plot_adding_days_active[,"days_active"]
}

# Taxonomic diversity metrics ##################################################

library(hillR) # A package for computing diversity indices

# Total abundance (***incorporating the adjustment for when 1 cup lost***)
carab2_by_plot$abundance <- rowSums(carab2_by_plot[,carab_species])
hist(carab2_by_plot$abundance, breaks=seq(0,100,5)) # The total number of 
# carabids caught at each Plot across the summer of 2022
carab2_by_plot$abundance_per_day <- carab2_by_plot$abundance / 
  carab2_by_plot$days_active
ggplot(data=carab2_by_plot, aes(x=Treatment, y=abundance_per_day)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() + ylim(0,0.8)+
  ylab("Ground beetles caught per day") + xlab("Forest disturbance")

# What is the coefficient of variation among the beetles/day at 
# each plot?
# (the ratio of the standard deviation to the mean)
sd(carab2_by_plot$abundance_per_day) / 
  mean(carab2_by_plot$abundance_per_day) 

# Hill Numbers
# q = 0 (default) to get species richness, 
# q = 1 to get shannon entropy,
# q = 2 will give inverse Simpson.
# MARGIN = 1 is the default, indicates that sites are rows

# species richness (the number of species found at each plot over the season)
carab2_by_plot$richness_test <- rowSums(carab2_by_plot[,carab_species] > 0)
carab2_by_plot$richness <- hill_taxa(carab2_by_plot[,carab_species], q = 0, 
                                     MARGIN = 1)
hist(carab2_by_plot$richness, breaks=seq(0,20,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=richness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() + 
  ylab("Carabid species richness") + scale_y_continuous(breaks = seq(5,18,1))

# Shannon diversity: exp(H')
carab2_by_plot$shannon_diversity <- hill_taxa(carab2_by_plot[,carab_species], 
                                              q = 1, MARGIN = 1)
hist(carab2_by_plot$shannon_diversity, breaks=seq(0,15,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() +
  ylab("Shannon diversity: exp(H)")

# Simpson diversity: inverse Simpson index (1/D)
carab2_by_plot$simpson_diversity <- hill_taxa(carab2_by_plot[,carab_species], 
                                              q = 2, MARGIN = 1)
hist(carab2_by_plot$simpson_diversity, breaks=seq(0,15,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=simpson_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() + 
  ylab("Inverse Simpson index")

# evenness
library(chemodiv)
#  Equal to the Hill diversity at a value of q, divided by the species richness.
# Min value is 1/species richness and max value is 1.
carab2_by_plot$evenness_test <- carab2_by_plot$simpson_diversity / 
                                                    carab2_by_plot$richness
carab2_by_plot$evenness <- calcDiv(carab2_by_plot[,carab_species], type = "HillEven", 
                                   q = 2)$HillEven
hist(carab2_by_plot$evenness)
ggplot(data=carab2_by_plot, aes(x=Treatment, y=evenness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() +
  ylab("Simpson evenness")

# Assess sampling effort using accumulation curves #############################

# Individual-based rarefaction by treatment, jackknife estimates by treatment.
# We are separating out each treatment so we can get a species accumulation curve 
# for each, and then we will graph it

table(carab2_no_missing$Treatment)
TF <- carab2_no_missing[which(carab2_no_missing$Treatment == "F"),]
TS <- carab2_no_missing[which(carab2_no_missing$Treatment == "S"),]
TW <- carab2_no_missing[which(carab2_no_missing$Treatment == "W"),]

# The rarefaction method standardizes the sample sizes so that we are comparing 
# species richness at equivalent abundances
library(vegan)

sp.TF <- specaccum(TF[,carab_species], method = "rarefaction", 
                   permutations = 100, gamma = "jack2") 
sp.TS <- specaccum(TS[,carab_species], method = "rarefaction", 
                   permutations = 100, gamma = "jack2")
sp.TW <- specaccum(TW[,carab_species], method = "rarefaction", 
                   permutations = 100, gamma = "jack2")


# Make the plot using number of individuals on the x-axis:

par(mfrow=c(1,1))
par(mar=c(5,6,4,2))

plot(sp.TS, pch = 19, col = "brown4", xvar = c("individuals"), lty = 4, lwd = 2,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 400), ylim = c(0, 40))
plot(sp.TF, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 2, col = "palegreen4")
plot(sp.TW, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 2, col = "goldenrod2")

legend("bottomright", legend = c("Salvaged", "Forest", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 0.8, bty = "n", lwd = 5,
       col = c("brown4", "palegreen3", "goldenrod2"))

par(mfrow=c(1,1))
par(mar=c(5,6,4,2))

plot(sp.TS, pch = 19, col = "brown4", xvar = c("sites"), lty = 4, lwd = 2,
     ylab = "Species Richness", xlab = "Number of pitfalls", xlim = c(0, 100), ylim = c(0, 40))
plot(sp.TF, add = TRUE, pch = 15, xvar = c("sites"), lty = 1, lwd = 2, col = "palegreen4")
plot(sp.TW, add = TRUE, pch = 4, xvar = c("sites"), lty = 2, lwd = 2, col = "goldenrod2")

legend("bottomright", legend = c("Salvaged", "Forest", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 0.8, bty = "n", lwd = 5,
       col = c("brown4", "palegreen3", "goldenrod2"))

# Calculates the estimated species richness of a population using first- and 
# second-order jackknife estimators
# first-order jackknife estimates are based on the number of singletons
# second-order jackknife estimates are based on the number of singletons and doubletons

specnumber(carab2_no_missing[,carab_species], groups = carab2_no_missing$Treatment)

library(fossil)
# forest
jack1(TF[,carab_species], taxa.row = FALSE, abund = TRUE)
jack2(TF[,carab_species], taxa.row = FALSE, abund = TRUE)

# salvaged
jack1(TS[,carab_species], taxa.row = FALSE, abund = TRUE)
jack2(TS[,carab_species], taxa.row = FALSE, abund = TRUE)

# windthrow
jack1(TW[,carab_species], taxa.row = FALSE, abund = TRUE)
jack2(TW[,carab_species], taxa.row = FALSE, abund = TRUE)

# Chao1 estimator ##############################################################
# Below, calculate the Chao1 estimator of species richness, using the SpadeR
# package:

carab2_by_treatment <- carab2_by_plot %>% group_by(Treatment) %>%
  summarize(across(all_of(carab_species), ~sum(.)))

forest_counts <- carab2_by_treatment[carab2_by_treatment$Treatment=="Forest", 
                                    carab_species]
windthrow_counts <- carab2_by_treatment[carab2_by_treatment$Treatment=="Windthrow", 
                                        carab_species]
salvaged_counts <- carab2_by_treatment[carab2_by_treatment$Treatment=="Salvaged", 
                                       carab_species]

library(SpadeR)
# Undisturbed forest:
ChaoSpecies(forest_counts, datatype = "abundance", k = 10, conf=0.95)

# Windthrow:
ChaoSpecies(windthrow_counts, datatype = "abundance", k = 10, conf=0.95)

# Salvaged:
ChaoSpecies(salvaged_counts, datatype = "abundance", k = 10, conf=0.95)

# Beta-diversity ##############################################################
# Compute a dissimilarity matrix for abundance-based data

# Note: The dataset I'm using here are counts per day of the different species.
dis.matrix <- vegdist(carab2_by_plot_beetles_per_day[,carab_species], 
                      method = "bray") # Bray-Curtis dissimilarity
dis.matrix 
max(dis.matrix)
min(dis.matrix) # The values of the dissimilarity matrix are between 
# 0.32 and 0.97

# Run the Permutation-based Multivariate Analysis of Variance (test for
# treatment differences between the centroids of each treatment group in 
# species space):
adonis2(dis.matrix ~ carab2_by_plot_beetles_per_day$Treatment, permutations = 999)
# There are not significant differences in the centroids of each treatment group

# Run the Analysis of Multivariate Homogeneity of Group Dispersions
# Finds the dispersion of each treatment from its spatial median. Then tests
# if treatments differ in their dispersions.
beta_dispersion <- betadisper(d = dis.matrix, 
                              group = carab2_by_plot_beetles_per_day$Treatment, 
                              type = c("median"))
beta_dispersion
plot(beta_dispersion)
boxplot(beta_dispersion, ylab = "Distance to median", xlab="")
anova(beta_dispersion) # The three treatment groups differ in their dispersions

# Note: it would not be surprising if the forest had more dispersion than the 
# salvage or windthrow, because there were more plots in the forest.

TukeyHSD(beta_dispersion, which = "group", conf.level = 0.95)
# The windthrow and salvaged treatments do not SIGNIFICANTLY 
# differ in their dispersion.

# Run the nonmetric multidimensional scaling visualization:
nmds.carabid <- metaMDS(dis.matrix, trymax = 500, k = 2)
nmds.carabid # stress is quality of fit
stressplot(nmds.carabid)
plot(nmds.carabid) # basic plot with no treatment distinctions

# plot the NMDS model by treatment
ordiplot(nmds.carabid, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Treatment=="Forest"), pch = 15, cex = 1, col = "palegreen4")
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Treatment=="Salvaged"), pch = 16, cex = 1, col = "brown4")
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Treatment=="Windthrow"), pch = 17, cex = 1, col = "goldenrod2")

ordiellipse(nmds.carabid, carab2_by_plot$Treatment, draw = "lines", col = c("palegreen4", "brown4", "goldenrod2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("topleft", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("palegreen4", "brown4", "goldenrod2"))


# plot the NMDS model by area:
ordiplot(nmds.carabid, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Area=="northeast"), pch = 15, cex = 1, col = "blue")
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Area=="southwest"), pch = 16, cex = 1, col = "purple")

ordiellipse(nmds.carabid, carab2_by_plot$Area, draw = "lines", col = c("blue", "purple"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("topleft", legend = c("northeast", "southwest"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("blue", "purple"))

# plot the NMDS model with plot 49 highlighted. Exploration of the environmental
# dataset indicated that plot 49 has a very open canopy, with high ground 
# vegetation cover and high ground vegetation height
ordiplot(nmds.carabid, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Plot!=49), pch = 15, cex = 1, col = "black")
points(nmds.carabid, dis = "sites", select = which(carab2_by_plot$Plot==49), pch = 16, cex = 1, col = "red")
# A quick look at the carabid data shows that plot 49 had the highest captures of 
# Pterostichus stygicus, Agonum fidele, Anisodactylus nigerrimus, melanopus, and 
# harrisi, Chlaenius emarginatus, and maybe others (need to investigate more)

# In the paper Perry et al. 2018, Chlaenius emarginatus, Cyclotrachelus sigillatus,
# and Pterostichus stygicus are found in higher abundance in canopy gaps. This matches
# with what I'm finding at Plot 49. Additionally, Agonum species are typically
# indicative of water-saturated soils (source: Larochelle and Lariviere as well
# as the paper on impacts of EAB on ground beetles)

# Investigating the most common species ##########################################

# There were 46 species found in 2022

# What were the most commonly caught carabids?

colSums(carab2_by_plot[, carab_species])[order(colSums(carab2_by_plot[, carab_species]), decreasing=T)]
# Pterostichus adoxus (120), Sphaeroderus stenostomus (115), 
# Pterostichus stygicus (90), Platynus angustatus (77), Pterostichus tristis (74),
# Pterostichus lachrymosus (68), 
# Cyclotrachelus_unknown_likely_sigillatus_or_convivus_or_mixed (54),
# P. moestus (50), Carabus goryi (43)

# How many species were common (>= 30 individuals collected) vs rare (< 30 
# individuals)
sum(colSums(carab2_by_plot[, carab_species]) >= 30) # 9 species
sum(colSums(carab2_by_plot[, carab_species]) < 30) # 37 species

# Graph activity-abundance of Pterostichus adoxus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Pterostichus_adoxus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Sphaeroderus stenostomus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Sphaeroderus_stenostomus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Pterostichus stygicus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Pterostichus_stygicus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Platynus angustatus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Platynus_angustatus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Pterostichus tristis:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Pterostichus_tristis)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Pterostichus lachrymosus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Pterostichus_lachrymosus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Cyclotrachelus_unknown_likely_sigillatus_or_convivus_or_mixed:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Cyclotrachelus_unknown_likely_sigillatus_or_convivus_or_mixed)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Pterostichus moestus:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Pterostichus_moestus)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph activity-abundance of Carabus goryi:
ggplot(data=carab2_by_plot, aes(x=Treatment, y=Carabus_goryi)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Export data table ###########################################################

# Here, I'll write a new csv file with the 24 plots as rows

write.csv(carab2_by_plot, "Aaron_PNR_formatted_csvs/PNR2022_carabids_by_plot.csv",
          row.names=FALSE)


