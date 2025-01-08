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
carab0 <- read.csv("Aaron_PNR_formatted_csvs/PNR2022_carabid_counts.csv",
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
carab0$Set_Date <- as.factor(carab0$Set_Date)
carab0$Collection_Date <- as.factor(carab0$Collection_Date)

carab0 <- carab0 %>% select(-Treatment, -PNR_Code) # get rid of unnecessary columns

# Graph the number of carabids captured, with trap location in the x-axis:
ggplot(data=carab0, mapping=aes(x=Plot, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()
# Now with collection interval on the x-axis:
ggplot(data=carab0, mapping=aes(x=Collection_interval, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_jitter(height=0, width=0.1, alpha=0.5) +
  theme_classic()

# Import information about the trap locations, such as transect and area:
trap_locations <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("factor", "integer", "factor", "factor", "factor", "numeric", "numeric"))
library(forcats)
trap_locations <- mutate(trap_locations, Treatment = fct_recode(Treatment, 
                            "Forest" = "F", "Windthrow" = "W", "Salvaged" = "S"))

# Merge the carab dataframe with the trap_locations dataframe:
carab1 <- right_join(trap_locations, carab0, by="Plot") %>% arrange(Collection_interval)
table(carab1$Plot)
table(carab1$Collection_interval)

# Look at which treatments have missing data:
carab1[is.na(carab1$Agonum_fidele), c("Collection_interval", "Plot", "Treatment")] 
# The missing rows are all from the undisturbed forest treatment

# Deal with missing data and dumped traps ######################################

# Adjust carabid counts based on dumped traps. Each pitfall has two collection 
# cups, and sometimes one or both of those cups was dumped. 2/2 cups were dumped
# four times (resulting in missing data). 1/2 cups were dumped 12 times (1 of
# those 12 is also missing data due to losing the vial):
carab1[carab1$Traps_Dumped_copied_from_PNR_ENV_2022 == 1, c("Plot", 
                                                          "Collection_interval",
                                                          "Total_Carabidae_from_sum_of_species_counts")]
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
days_active <- data.frame(Plot=c( 41,  42,  43,  44,     45,  46,  47,     48,  49,    50,  51,     52,  53,     54,  55,  56,  57,  58,  59,  60,  61,        62,  63,  64), 
                   days_active=c(112, 111, 111, 112, 111-15, 111, 111, 111-13, 111,111-15, 111, 111-13, 112, 112-14, 112, 112, 112, 112, 112, 112, 112, 111-13-14, 112, 112))
days_active$Plot <- as.factor(days_active$Plot)

# Use dplyr::summarize to create a new dataframe with a row for every Plot. The
# columns will be the total number of each species captured, across all 
# collection intervals.

carab2_by_plot_init <- carab2 %>% group_by(Plot) %>% 
  summarise(across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

carab2_by_plot_adding_plot_info <- right_join(trap_locations, 
                                              carab2_by_plot_init, by="Plot")

carab2_by_plot_adding_days_active <- right_join(days_active, carab2_by_plot_adding_plot_info, 
                             by="Plot")

carab2_by_plot <- carab2_by_plot_adding_days_active
# The below code divides abundances of each species at each trap by the days
# active. However, I commented it out because I'm not sure I want to convert 
# the count data to fractions.
#for (species in carab_species) {
#  carab2_by_plot[,species] <- carab2_by_plot_adding_days_active[,species] / 
#                              carab2_by_plot_adding_days_active[,"days_active"]
#}

# Taxonomic diversity metrics ##################################################

library(hillR)

# Total abundance:
carab2_by_plot$abundance <- rowSums(carab2_by_plot[,carab_species])
hist(carab2_by_plot$abundance, breaks=seq(0,100,5)) # The total number of 
# carabids caught at each Plot across the summer of 2022
ggplot(data=carab2_by_plot, aes(x=Treatment, y=abundance/days_active)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic() + ylim(0,0.8)+
  ylab("Ground beetles caught per day") + xlab("Forest disturbance")

# Hill Numbers
# q = 0 (default) to get species richness, 
# q = 1 to get shannon entropy,
# q = 2 will give inverse Simpson.
# MARGIN = 1 is the default, indicates that sites are rows

# species richness
carab2_by_plot$richness <- hill_taxa(carab2_by_plot[,carab_species], q = 0, 
                                     MARGIN = 1)
hist(carab2_by_plot$richness, breaks=seq(0,20,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=richness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic()

# Shannon diversity: exp(H')
carab2_by_plot$shannon_diversity <- hill_taxa(carab2_by_plot[,carab_species], 
                                              q = 1, MARGIN = 1)
hist(carab2_by_plot$shannon_diversity, breaks=seq(0,15,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic()

# Simpson diversity: inverse Simpson index (1/D)
carab2_by_plot$simpson_diversity <- hill_taxa(carab2_by_plot[,carab_species], 
                                              q = 2, MARGIN = 1)
hist(carab2_by_plot$simpson_diversity, breaks=seq(0,15,1))
ggplot(data=carab2_by_plot, aes(x=Treatment, y=simpson_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic()

# evenness
library(chemodiv)
#  For q = 2, more weight is put on compounds with high proportions (I don't get this)
carab2_by_plot$evenness <- calcDiv(carab2_by_plot[,carab_species], type = "HillEven", 
                                   q = 2)$HillEven
hist(carab2_by_plot$evenness)
ggplot(data=carab2_by_plot, aes(x=Treatment, y=evenness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + theme_classic()

#Estimate species richness with accumulation curves#####################################################################

# Individual-based rarefaction by treatment, jackknife estimates by treatment.
# We are separating out each treatment so we can get a species accumulation curve for each, and then we will graph it

# Note: I substituted the dataframe where missing carabid rows are removed:

table(carab2_no_missing$Treatment)
TF <- carab2_no_missing[which(carab2_no_missing$Treatment == "Forest"),]
TS <- carab2_no_missing[which(carab2_no_missing$Treatment == "Salvaged"),]
TW <- carab2_no_missing[which(carab2_no_missing$Treatment == "Windthrow"),]

# the rarefaction method standardizes the sample sizes so that we are comparing 
# species richness at equivalent abundances
library(vegan)

sp.TF <- specaccum(TF[,carab_species], method = "rarefaction", permutations = 100, gamma = "jack2") 
sp.TS <- specaccum(TS[,carab_species], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.TW <- specaccum(TW[,carab_species], method = "rarefaction", permutations = 100, gamma = "jack2")


# make the plot

#png("Rarefaction.png", width = 1200, height = 1000, pointsize = 30)

par(mfrow=c(1,1))
par(mar=c(5,6,4,2))

plot(sp.TS, pch = 19, col = "brown4", xvar = c("individuals"), lty = 4, lwd = 2,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 400), ylim = c(0, 40))
plot(sp.TF, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 2, col = "palegreen4")
plot(sp.TW, add = TRUE, pch = 4, xvar = c("individuals"), lty = 2, lwd = 2, col = "goldenrod2")

legend("bottomright", legend = c("Salvaged", "Forest", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1.2, bty = "n", lwd = 5,
       col = c("brown4", "palegreen3", "goldenrod2"))

#dev.off()


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


# Nonmetric multidimensional scaling (NMDS) ###################################
# Compute a dissimilarity matrix
# The method option let's you indicate which dissimilarity metric to calculate
# We will calculate the bray-curtis dissimilarity matrix for abundance-based data
dis.matrix <- vegdist(carab2_by_plot[,carab_species], method = "bray")
dis.matrix

# Run the nonmetric multidimensional scaling model
nmds.carabid <- metaMDS(dis.matrix, trymax = 500, autotransform = TRUE, k = 2)
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


# Test for differences in carabid composition among treatments
# PERMANOVA tests whether the group centroid of communities differs among groups
# in multivariate space (e.g. different community composition)
# adonis2(dis.matrix ~ carab$Treatment, permutations = 999)

# BETADISPER tests whether the dispersion of a group from its spatial median is different
# between groups (i.e. species redundancy across space)
# analysis of multivariate homogeneity of group dispersions (variances)
# multivariate analogue of Levene's test for homogeneity of variances
#c.beta <- betadisper(dis.matrix, carab$Treatment, type = c("median"))
#anova(c.beta)
#plot(c.beta)
#boxplot(c.beta, ylab = "Distance to median")
#TukeyHSD(c.beta, which = "group", conf.level = 0.95)

# Investigating the most common species ##########################################

# What were the most commonly caught carabids?

colSums(carab2_by_plot[, carab_species])[order(colSums(carab2_by_plot[, carab_species]), decreasing=T)]





