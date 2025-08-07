# Aaron Tayal
# 3/30/2025
# Powdermill ground beetles stats analysis for Aaron's MS Thesis

library(ggplot2)
theme_set(theme_classic())
library(dplyr)
library(lubridate)
library(vegan) # For functions specaccum
library(hillR) # for Shannon diversity
library(SpadeR) # for Chao estimator
library(picante) # used for finding mean pairwise distance in trait space,
# and Rao's pairwise distance
library(geometry) # for dot product
library(FD) # for community-weighted mean
#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

carab_0 <- read.csv("Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts.csv")

# Fixing column classes ########################################################

carab_0$Year <- as.factor(carab_0$Year)
carab_0$Set_date <- mdy(carab_0$Set_date)
carab_0$Collection_date <- mdy(carab_0$Collection_date)
carab_0$Plot <- as.integer(carab_0$Plot)
carab_0$Treatment <- as.factor(carab_0$Treatment)
carab_0$Interval <- as.integer(carab_0$Interval)

# Make a list of the species found (using the column names)
carab_species <- colnames(carab_0)[which(colnames(carab_0) == "Agonoleptus_thoracicus"):dim(carab_0)[2]] # 55 species

# Import plot information ######################################################

plot_locations <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv")
library(forcats)
plot_locations <- mutate(plot_locations, Treatment = 
                           fct_recode(Treatment, "Forest" = "F", 
                                      "Windthrow" = "W", "Salvaged" = "S"))
plot_locations$Plot <- as.integer(plot_locations$Plot)
plot_info_columns <- colnames(plot_locations)
# Note: there are 25 rows in the plot_locations data table, because plot 65 
# was used in 2015 but plot 63 was used instead of 65 in 2022. This is because
# there was a reason why the plot needed to be moved by ~27 m

# IMPORTANT NOTE: Right now, I changed the plot number for the 2015 carabid
# counts dataset to plot 63. This facilitates comparisons between years.
# But in reality, the trap was set out at plot 65

# Investigate accumulation of ground beetle species over the season ############

carab_by_interval_2015 <- carab_0 %>% filter(Year==2015) %>% group_by(Interval) %>%
  summarise(num_nonmissing_plots = n() - sum(is.na(Agonoleptus_thoracicus)),
            Set_date = first(Set_date), Collection_date = first(Collection_date),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

carab_by_interval_2022 <- carab_0 %>% filter(Year==2022) %>% group_by(Interval) %>%
  summarise(num_nonmissing_plots = n() - sum(is.na(Agonoleptus_thoracicus)),
            Set_date = first(Set_date), Collection_date = first(Collection_date),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

# The specaccum function creates a species accumulation curve. The "collector"
# method means the species are accumulated in the order of the trapping
# intervals: in other words, over the summer season.
plot(specaccum(comm=carab_by_interval_2015[,carab_species], method="collector"),
     xlab="Interval") + title("2015 accumulation of species over the season") 

plot(specaccum(comm=carab_by_interval_2022[,carab_species], method="collector"),
     xlab="Interval") + title("2022 accumulation of species over the season") 

# Additional species were caught in 2022 in September:
September_only_species <- 
  names(which((colSums(carab_by_interval_2022[7:8,carab_species]) > 0) & 
                (colSums(carab_by_interval_2022[1:6,carab_species]) == 0)))

#View(carab_by_interval_2022[, September_only_species])
#View(carab_by_interval_2015[, September_only_species])
# Of the nine species found only after Aug 11 in the 2022 sampling, most (8/9) 
# were not found at all in 2015. These include "Amerizus_unknown","Myas_coracinus",
# "Patrobus_longicornis"  ,"Platynus_hypolithos","Pterostichus_atratus", "Scaphinotus_andrewsii"  
# "Scaphinotus_ridingsii" and "Synuchus_impunctatus"

last_2_intervals_2015 <- 
  names(which((colSums(carab_by_interval_2015[5:6,carab_species]) > 0) & 
                (colSums(carab_by_interval_2015[1:4,carab_species]) == 0)))

# Make a graph of the captures of Platynus angustatus over the season:
ggplot(data=carab_by_interval_2015, aes(x=Collection_date, y=Platynus_angustatus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Collection_date, y=Platynus_angustatus)) +
  geom_line() + geom_point() + 
  ylab(expression(paste("Number of  ", italic("Platynus angustatus")))) + 
  xlab("Collection date")+
  ggtitle(expression(paste("Captures of  ", italic("Platynus angustatus"), " over the 2022 season")))

ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Pterostichus_tristis)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Pterostichus_tristis)) +
  geom_line()

ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Pterostichus_moestus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Pterostichus_moestus)) +
  geom_line()

ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Dicaelus_politus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Dicaelus_politus)) +
  geom_line()

ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Chlaenius_emarginatus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Chlaenius_emarginatus)) +
  geom_line()

ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Sphaeroderus_stenostomus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Sphaeroderus_stenostomus)) +
  geom_line()

# Investigate carabid abundance over the season ################################

carab_by_interval_2015$total_count_stdz <- 
  rowSums(carab_by_interval_2015[,carab_species]) / 
  carab_by_interval_2015$num_nonmissing_plots
# Gives the number of carabids caught at an average plot during the ~14 day interval

plot(carab_by_interval_2015$Interval, carab_by_interval_2015$total_count_stdz)

carab_by_interval_2022$total_count_stdz <- 
  rowSums(carab_by_interval_2022[,carab_species]) / 
  carab_by_interval_2022$num_nonmissing_plots

plot(carab_by_interval_2022$Interval, carab_by_interval_2022$total_count_stdz)

# Data standardization - by plot ###############################################

# First, calculate totals across all sampling intervals for each plot in each year

# In 2015, it seems the interval 7/22-8/5 had a lot of traps lost.
# 2015 had 6 total intervals, and the traps were out from 5/27 to 8/17/2015, 
# which is 82 days: 
yday("2015-8-17") - yday("2015-5-27")
82/6 # **About** 14 days per interval and *about* 84 days per season.
# Actually it was the last interval (8/5-8/17/2015) that only has 12 days rather 
# than 14.

yday("2022-9-20") - yday("2022-6-1") # 111 days from June 1st to Sept 20th
yday("2022-8-23") - yday("2022-6-1") # 83 days from June 1st to Aug 23

# When creating the plot-level dataframe, I'm choosing to exclude the final
# two trap intervals of the 2022 data (the ones collected in September)
g_0 <- carab_0 %>% filter(month(Collection_date) != 9) %>% group_by(Year, Plot) %>%
  summarise(num_nonmissing_intervals = n() - sum(is.na(Agonoleptus_thoracicus)),
            approx_trap_days = 14 * num_nonmissing_intervals,
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

# Now join the information about the plot locations:
g <- right_join(plot_locations, g_0, by="Plot") %>% arrange(Year)

# Now divide the species counts by the approximate number of days active out 
# of 84 days (to account for some trap contents being lost for a few plots)
g_stdz <- g
g_stdz[,carab_species] <- g[,carab_species] * 84 / g$approx_trap_days

# Investigate total counts of ground beetles ##################################
# How many ground beetles were collected?

g$total_count <- rowSums(g[,carab_species])
sum(g %>% filter(Year==2015) %>% select(total_count)) # 934 ground beetles in 2015
sum(g %>% filter(Year==2022) %>% select(total_count)) # 603 ground beetles in 
# the first six intervals of 2022.

g_stdz$total_count_stdz <- rowSums(g_stdz[,carab_species]) # totals for the standardized counts
sum(g_stdz %>% filter(Year==2015) %>% select(total_count_stdz)) 
sum(g_stdz %>% filter(Year==2022) %>% select(total_count_stdz)) 


# Investigate number of species found #########################################

# Make a list of species found in each year (including all intervals)
species_2015 <- names(which((colSums(carab_by_interval_2015[,carab_species]) > 0)))
species_2022 <- names(which((colSums(carab_by_interval_2022[,carab_species]) > 0)))
species_2022_excluding_last_2_int <- names(which((colSums(carab_by_interval_2022[1:6,carab_species]) > 0)))
species_2022_last_2_int <- names(which((colSums(carab_by_interval_2022[7:8,carab_species]) > 0)))

# How many species were found in both years?
base::intersect(species_2015, species_2022)
# Overall, 28 species were found in both years

# How many species were found only in 2022?
base::setdiff(species_2022, species_2015)
# 18 species were found in 2022 but not in 2015

# How many spp were found only in 2015?
base::setdiff(species_2015, species_2022)
# 9 species were found in 2015 but not 2022

# How many species were found in either 2015 and/or in 2022 in the first 6 intervals?
species_2015_2022_1st_6 <- union(species_2015, species_2022_excluding_last_2_int)

# Make a table of the counts for each species and each year:
count_table <- data.frame(species = names(colSums(carab_by_interval_2015[,carab_species])),
                          counts_2015 = unname(colSums(carab_by_interval_2015[,carab_species])),
                          counts_2022 = unname(colSums(carab_by_interval_2022[,carab_species])),
                          counts_2022_excluding_last_2_int = unname(colSums(carab_by_interval_2022[1:6,carab_species])),
                          counts_2022_last_2_int = unname(colSums(carab_by_interval_2022[7:8,carab_species])))

#write.csv(count_table, "Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts_summary.csv")


# Species accumulation curves #################################################
# Using the raw, unstandardized counts of carabids
# Exclude the final two sample intervals from 2022.

# Create subsets of rows for each treatment:
forest_2015 <- g %>% filter(Treatment == "Forest", Year==2015)
salvaged_2015 <- g %>% filter(Treatment == "Salvaged", Year==2015)
windthrow_2015 <- g %>% filter(Treatment == "Windthrow", Year==2015)

forest_2022 <- g %>% filter(Treatment == "Forest", Year==2022)
salvaged_2022 <- g %>% filter(Treatment == "Salvaged", Year==2022)
windthrow_2022 <- g %>% filter(Treatment == "Windthrow", Year==2022)

# Now create accumulation curves:
forest_2015_accum <- specaccum(comm=forest_2015[carab_species], method="random", permutations=100)
salvaged_2015_accum <- specaccum(comm=salvaged_2015[carab_species], method="random", permutations=100)
windthrow_2015_accum <- specaccum(comm=windthrow_2015[carab_species], method="random", permutations=100)

# Make a graph for 2015:
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
plot(forest_2015_accum, pch = 19, col = "palegreen3", xvar = c("sites"), lty = 4, 
     lwd = 2, ylab = "Species Richness", xlab = "Number of plots", 
     xlim = c(0, 12), ylim = c(0, 35))
plot(salvaged_2015_accum, add = TRUE, pch = 15, xvar = c("sites"), lty = 1, 
     lwd = 2, col = "goldenrod2")
plot(windthrow_2015_accum, add = TRUE, pch = 4, xvar = c("sites"), lty = 2, 
     lwd = 2, col = "brown4")

legend("bottomright", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 0.6, bty = "n", lwd = 5,
       col = c("palegreen3", "goldenrod2", "brown4"))
title("2015 ground beetle species accumulation")

# Now for 2022:

forest_2022_accum <- specaccum(comm=forest_2022[carab_species], method="random", permutations=100)
salvaged_2022_accum <- specaccum(comm=salvaged_2022[carab_species], method="random", permutations=100)
windthrow_2022_accum <- specaccum(comm=windthrow_2022[carab_species], method="random", permutations=100)

# Make a graph for 2022:
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
plot(forest_2022_accum, pch = 19, col = "palegreen3", xvar = c("sites"), lty = 4, 
     lwd = 2, ylab = "Species Richness", xlab = "Number of plots", 
     xlim = c(0, 12), ylim = c(0, 35))
plot(salvaged_2022_accum, add = TRUE, pch = 15, xvar = c("sites"), lty = 1, 
     lwd = 2, col = "goldenrod2")
plot(windthrow_2022_accum, add = TRUE, pch = 4, xvar = c("sites"), lty = 2, 
     lwd = 2, col = "brown4")

legend("bottomright", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 0.6, bty = "n", lwd = 5,
       col = c("palegreen3", "goldenrod2", "brown4"))
title("2022 ground beetle species accumulation")


# Taxonomic alpha-diversity: Chao estimators ###################################
# still omitting the last two intervals of the 2022 data

total_counts_2015 <- colSums(g %>% filter(Year==2015) %>% select(all_of(carab_species)))
total_counts_2022 <- colSums(g %>% filter(Year==2022) %>% select(all_of(carab_species)))

# Total counts 2015:
ChaoSpecies(total_counts_2015, datatype = "abundance", k = 10, conf=0.95)
# 37 observed species
# Chao1 estimates 47.114 species. 95% conf. int: 39.262 -> 82.233

# Total counts 2022:
ChaoSpecies(total_counts_2022, datatype = "abundance", k = 10, conf=0.95)
# 37 observed species
# Chao1 estimates 79.180 species. 95% conf. int: 46.183 -> 230.745


# Taxonomic alpha-diversity: plot level species richness and Shannon diversity ####

g_stdz$sp_rich <- rowSums(g_stdz[,carab_species] > 0)

# Plot species richness as a function of treatment:
ggplot(data=g_stdz %>% filter(Year==2015), aes(x=Treatment, y=sp_rich)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Number of species") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle Species Richness")

ggplot(data=g_stdz %>% filter(Year==2022), aes(x=Treatment, y=sp_rich)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Number of species") + xlab("Forest disturbance") + 
  ggtitle("2022 Ground Beetle Species Richness")

# Hill Numbers
# q = 0 (default) to get species richness, 
# q = 1 to get shannon entropy,
# q = 2 will give inverse Simpson.
# MARGIN = 1 is the default, indicates that sites are rows

# Shannon diversity: exp(-Σp_i*log(p_i))
g_stdz$shannon_diversity <- 
  hill_taxa(g_stdz[,carab_species], q = 1, MARGIN = 1)

# Plot Shannon diversity as a function of treatment:
ggplot(data=g_stdz %>% filter(Year==2015), aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Shannon diversity") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle Shannon diversity")

ggplot(data=g_stdz %>% filter(Year==2022), aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + ylim(0,15) +
  ylab("Shannon diversity") + xlab("Forest disturbance") +
  ggtitle("2022 Ground Beetle Shannon diversity")


# Functional alpha diversity ###################################################
# I will calculate functional alpha diversity using the weighted mean pairwise 
# distance between species found at a plot. The distance is the distance in 
# trait space (larger distance = traits are more distict). Distances will
# be weighted by the product of the relative abundances of the two species.

# Important note: as the species richness increases, the mean pairwise distance
# tends to become less variable (source: Swenson)

dist_0 <- read.csv("Aaron_PNR_formatted_data/carabid_dist_in_trait_space.csv")

# Check if the species column names match between the distance matrix and the
# species abundance data:
all.equal(colnames(dist_0), carab_species) # True

dist <- as.matrix(dist_0)
rownames(dist) <- colnames(dist)

g_stdz_matrix <- as.matrix(g_stdz[,carab_species])

# Compute the weighted mean_pairwise_distance for each plot
g_stdz$mean_pairwise_distance <-
  picante::mpd(samp = g_stdz_matrix,
               dis = dist,
               abundance.weighted = TRUE)

# How does mean pairwise distance relate to spp richness?
ggplot(data=g_stdz) + geom_point(aes(x=sp_rich, y=mean_pairwise_distance))
# The mean_pairwise_distance seems to be positively correlated with spp richness

# How does mean pairwise distance relate to activity-abundance?
ggplot(data=g_stdz) + geom_point(aes(x=total_count_stdz, y=mean_pairwise_distance))
# Does not seem to be super related to activity-abundance

# Plot functional alpha-diversity as a function of treatment:
ggplot(data=g_stdz %>% filter(Year==2015), aes(x=Treatment, y=mean_pairwise_distance)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,0.3)+
  ylab("Mean pairwise distance") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle functional alpha-diversity")

ggplot(data=g_stdz %>% filter(Year==2022), aes(x=Treatment, y=mean_pairwise_distance)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,0.3)+
  ylab("Mean pairwise distance") + xlab("Forest disturbance") + 
  ggtitle("2022 Ground Beetle functional alpha-diversity")

# Taxonomic beta-diversity ###################################################

# First, I need to convert my data table to relative abundance data:
g_stdz_rel <- vegan::decostand(g_stdz[carab_species], method = "total")

# Now compute distance matrix between all 48 plots:
dist_spp_space <- vegdist(g_stdz_rel, method = "bray")
summary(as.vector(dist_spp_space))
hist(as.vector(dist_spp_space))
# A distance of 1 indicates complete dissimilarity.
# I guess if there is no commonality in the captured species, then the dissimilarity
# is 1.

# Now create nonmetric multidimensional scaling ordinations:
nmds <- metaMDS(dist_spp_space, trymax = 500, k = 2)
nmds # stress is quality of fit
stressplot(nmds)
g_stdz$NMDS1 <- nmds$points[,1]
g_stdz$NMDS2 <- nmds$points[,2]

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Year)) + 
  geom_point(size=2)

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment)) + 
  geom_point(size=2) + scale_color_manual(values = c("Forest" = "palegreen3", 
                          "Salvaged" = "goldenrod2", "Windthrow" = "brown4"))

# Now run the Permutational Multivariate Analysis of Variance:
adonis2(dist_spp_space ~ g_stdz$Treatment + g_stdz$Year + g_stdz$Treatment * g_stdz$Year, 
        permutations = 999)

adonis2(dist_spp_space ~ g_stdz$Year + g_stdz$Treatment + g_stdz$Treatment * g_stdz$Year, 
        permutations = 999) # Changing the ordering of the explanatory variables 
# changes the p-values

# Pairwise test to determine if any of the treatment groups differ:
pairwise.adonis(dist_spp_space, g_stdz$Treatment)

# Run the Analysis of Multivariate Homogeneity of Group Dispersions:

beta_dispersion_treatment <- betadisper(d = dist_spp_space, 
                              group = g_stdz$Treatment, 
                              type = c("median"))
beta_dispersion_treatment
plot(beta_dispersion_treatment)
boxplot(beta_dispersion_treatment, ylab = "Distance to median", xlab="")
anova(beta_dispersion_treatment)

beta_dispersion_year <- betadisper(d = dist_spp_space, 
                                        group = g_stdz$Year, 
                                        type = c("median"))
beta_dispersion_year
plot(beta_dispersion_year)
boxplot(beta_dispersion_year, ylab = "Distance to median", xlab="")
anova(beta_dispersion_year)

# Functional beta-diversity ####################################################

# First I need to calculate a dissimilarity matrix between all pairwise combinations
# of plots. Plots should be more similar if their species share similar traits,
# and more different if their species have different traits. I'll use the comdist
# function within the picante package:
dist_functional_beta <- comdist(comm = g_stdz_matrix, # the activity-abundance data for each plot
                                     dis = dist, # the distance between each carabid species in trait space
                                     abundance.weighted = T) # mean-pairwise distances will be weighted by species abundances

# Now create an NMDS ordination:
nmds_functional <- metaMDS(dist_functional_beta, trymax = 500, k = 2)
nmds_functional # stress is quality of fit
stressplot(nmds_functional)

g_stdz$NMDS_functional1 <- nmds_functional$points[,1]
g_stdz$NMDS_functional2 <- nmds_functional$points[,2]

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Year)) + 
  geom_point(size=2)

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Treatment)) + 
  geom_point(size=2) + scale_color_manual(values = c("Forest" = "palegreen3", 
                          "Salvaged" = "goldenrod2", "Windthrow" = "brown4"))


# Now run a PERMANOVA to test the null hypothesis that the centroids of each 
# treatment group are identical and their dispersions are identical:
adonis2(dist_functional_beta ~ g_stdz$Treatment + g_stdz$Year + 
          g_stdz$Treatment * g_stdz$Year, permutations = 999)

adonis2(dist_functional_beta ~ g_stdz$Year + g_stdz$Treatment + 
          g_stdz$Treatment * g_stdz$Year, permutations = 999)

# Do a pairwise test to see if any of the treatment groups are different:
pairwise.adonis(dist_functional_beta, g_stdz$Treatment)

# Now run an Analysis of Multivariate Homogeneity of Group Dispersions:
beta_dispersion_functional_year <- betadisper(d = dist_functional_beta, 
                                              group = g_stdz$Year, 
                                              type = c("median"))
beta_dispersion_functional_year
plot(beta_dispersion_functional_year)
boxplot(beta_dispersion_functional_year, ylab = "Distance to median", xlab="")
anova(beta_dispersion_functional_year)

beta_dispersion_functional_treatment <- betadisper(d = dist_functional_beta, 
                                              group = g_stdz$Treatment, 
                                              type = c("median"))
beta_dispersion_functional_treatment
plot(beta_dispersion_functional_treatment)
boxplot(beta_dispersion_functional_treatment, ylab = "Distance to median", xlab="")
anova(beta_dispersion_functional_treatment)

# Investigate activity-abundance of open-habitat, eurytopic, and ###############
# forest-specialist carabids 



# Import the trait data:
traits <- read.csv("Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv")
rownames(traits) <- traits$Species

# Of the 47 species caught in the first 6 collecion intervals, how many were 
# open-habitat and eurytopic?
table(traits %>% filter(Species %in% species_2015_2022_1st_6) %>% select(Forest_affinity))

# make lists of the eurytopic, open-habitat, and forest specialist species
open_habitat_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "open habitat", "Species"] # 2 spp
eurytopic_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "eurytopic", "Species"] # 28 spp
forest_specialist_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "forest specialist", "Species"] # 23 spp
unknown_forest_affinity_spp <- traits[is.na(traits$Forest_affinity), "Species"] 
# (2 species have unknown forest affinity)

# Create a column for activity-abundance of open habitat species
g_stdz$open_habitat_spp_stdz <- rowSums(g_stdz[, open_habitat_spp])

# And for eurytopic species:
g_stdz$eurytopic_spp_stdz <- rowSums(g_stdz[, eurytopic_spp])

# And for forest specialists:
g_stdz$forest_specialist_spp_stdz <- rowSums(g_stdz[, forest_specialist_spp])

# Community-weighted mean trait values #########################################

# I want to calculate community-weighted mean trait values for each
# plot
CWM <- function(community, trait_name){
  geometry::dot(c(t(community[carab_species] /
                      sum(community[carab_species]))),
                traits[,trait_name])
}

geometry::dot(c(t(g_stdz[1, carab_species] / 
                    g_stdz$total_count_stdz[1])), 
              traits$PC1) # weighted mean is simply a dot product of a vector 
# of the relative abundances of each species, with a vector of the trait value 
# for each species

# Test out the function to make sure it works:
CWM(community = g_stdz[1,], trait_name = "PC1")

# I decided to use the function "functcomp" in the package FD to calculate 
# CWMs. It looks like the results match what I got with my function. Here,
# I'll calculate the CWM for each of the first three PC axes, as well
# as the eight numerical traits and Water_affinity and Flight_capability
# and standardized antenna length
numeric_traits_modified <- c("body_length", "eye_length_standard",
                             "eye_protrusion_ratio", "pronotum_width_standard",
                             "abdomen_width_standard", "rear_leg_length_standard", 
                             "antenna_rear_leg_ratio",
                             "rear_trochanter_length_standard")
CWMs <- FD::functcomp(traits[,c("PC1", "PC2", "PC3", 
                                     numeric_traits_modified,
                                     "antenna_length_standard",
                                     "eye_protrusion_standard",
                                     "Water_affinity", "Flight_capability")], 
          as.matrix(g_stdz[, carab_species]), CWM.type = "all")

# Join to the data table:
g_stdz_with_CWMs <- cbind(g_stdz, CWMs)

# Export data tables ###########################################################
# I will run the statistical models in another R script, so I need to export 
# the plot-level data tables

#write.csv(g_stdz_with_CWMs, file="Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts_by_plot_standardized.csv", row.names = F)

# Observation: the capture rate of ground beetles in 2015 varied widely between
# plots. One plot caught a ton of Chlaenius emarginatus and Pterostichus
# moestus


