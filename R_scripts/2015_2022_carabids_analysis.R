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
library(picante) # used for finding mean pairwise distance in trait space
library(geometry) # for dot product
library(FD) # for community-weighted mean

carab_0 <- read.csv("Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts.csv")

# Fixing column classes ########################################################

carab_0$Year <- as.factor(carab_0$Year)
carab_0$Set_date <- mdy(carab_0$Set_date)
carab_0$Collection_date <- mdy(carab_0$Collection_date)
carab_0$Plot <- as.integer(carab_0$Plot)
carab_0$Treatment <- as.factor(carab_0$Treatment)
carab_0$Interval <- as.integer(carab_0$Interval)

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

# Data standardization #########################################################

# First, calculate totals across all sampling intervals, for each year:

# Make a list of the species found (using the column names)
carab_species <- colnames(carab_0)[which(colnames(carab_0) == "Agonoleptus_thoracicus"):dim(carab_0)[2]]

# For 2015:
carab_by_plot_2015_0 <- carab_0 %>% filter(Year==2015) %>% group_by(Plot) %>%
  summarise(num_nonmissing_intervals = n() - sum(is.na(Agonoleptus_thoracicus)),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))
# The column num_nonmissing intervals is the number of intervals, for a given
# plot, which have carabid data (the trap contents were not lost)

# Now join the 2015 carabid count data to the plot information dataframe:
carab_by_plot_2015_1 <- right_join(plot_locations, carab_by_plot_2015_0, by="Plot")

# In 2015, it seems the interval 7/22-8/5 had a lot of traps lost.
# 2015 had 6 total intervals, and the traps were out from 5/27 to 8/17/2015, 
# which is 82 days: 
yday(mdy("8-17-2015")) - yday(mdy("5-27-2015"))
# Meanwhile, 6*14 = 84. Actually it was the last interval (8/5-8/17/2015) that
# only has 12 days rather than 14. Thus, whenever all 6 intervals are not 
# missing, then it was 82 days. If 5 intervals are not missing, then it was
# 68 days.

carab_by_plot_2015_1$days_active <- 0
for (row in 1:24){
  if (carab_by_plot_2015_1$num_nonmissing_intervals[row] == 6){
    carab_by_plot_2015_1$days_active[row] <- 82
  } else {
    carab_by_plot_2015_1$days_active[row] <- 68
  }
} # Now we have a column with the number of days each trap was active for

# Rearrange columns:
carab_by_plot_2015 <- tibble(carab_by_plot_2015_1[,c(plot_info_columns, "days_active", 
                                              "num_nonmissing_intervals", carab_species)])

# Now for 2022:
carab_by_plot_2022_0 <- carab_0 %>% filter(Year==2022) %>% group_by(Plot) %>%
  summarise(num_nonmissing_intervals = n() - sum(is.na(Agonoleptus_thoracicus)),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

# The number of days each plot had an active trap:
days_active_2022 <- data.frame(Plot=c(   41,  42,  43,  44,  45,  46,  47,     48,  49,  50,  51,     52,  53,     54,  55,  56,  57,  58,  59,  60,  61,        62,  63,  64), 
                          days_active=c(112, 111, 111, 112, 111, 111, 111, 111-13, 111, 111, 111, 111-13, 112, 112-14, 112, 112, 112, 112, 112, 112, 112, 111-13-14, 112, 112))

carab_by_plot_2022_1 <- tibble(full_join(days_active_2022, carab_by_plot_2022_0,
                                by="Plot"))

# Join to the plot locations information:
carab_by_plot_2022 <- tibble(right_join(plot_locations, carab_by_plot_2022_1, 
                                       by="Plot")) 

# Now divide the species counts by the number of days active to get the 
# standardized counts:
# For 2015:
days_active_2015_rep <- as.data.frame(rep(carab_by_plot_2015[,"days_active"], 
                                          length(carab_species)))

carab_by_plot_2015_stdz <- carab_by_plot_2015

carab_by_plot_2015_stdz[,carab_species] <- 
  carab_by_plot_2015[,carab_species] / days_active_2015_rep

# For 2022:
days_active_2022_rep <- as.data.frame(rep(carab_by_plot_2022[,"days_active"], 
                            length(carab_species)))

carab_by_plot_2022_stdz <- carab_by_plot_2022

carab_by_plot_2022_stdz[,carab_species] <- 
  carab_by_plot_2022[,carab_species] / days_active_2022_rep


# Investigate total counts of ground beetles ##################################
# How many ground beetles were collected?

carab_by_plot_2015$total_count <- rowSums(carab_by_plot_2015[,carab_species])
sum(carab_by_plot_2015$total_count) # 934 ground beetles in 2015

carab_by_plot_2022$total_count <- rowSums(carab_by_plot_2022[,carab_species])
sum(carab_by_plot_2022$total_count) # 857 ground beetles in 2022

carab_by_plot_2015_stdz$total_count_stdz <- rowSums(carab_by_plot_2015_stdz[,carab_species])
mean(carab_by_plot_2015_stdz$total_count_stdz) # Each trap, on average, caught
# about 0.5 ground beetles per day in 2015
var(carab_by_plot_2015_stdz$total_count_stdz)

carab_by_plot_2022_stdz$total_count_stdz <- rowSums(carab_by_plot_2022_stdz[,carab_species])
mean(carab_by_plot_2022_stdz$total_count_stdz) # Each trap, on average, caught
# about 0.33 ground beetles per day in 2022


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
# Of the nine species found only after Aug 23 in the 2022 sampling, most (8/9) 
# were not found at all in 2015. These include "Amerizus_unknown","Myas_coracinus",
# "Patrobus_longicornis"  ,"Platynus_hypolithos","Pterostichus_atratus", "Scaphinotus_andrewsii"  
# "Scaphinotus_ridingsii" and "Synuchus_impunctatus"

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

# Investigate number of species found #########################################

# Make a list of species found in each year:
species_2015 <- names(which((colSums(carab_by_interval_2015[,carab_species]) > 0)))
species_2022 <- names(which((colSums(carab_by_interval_2022[,carab_species]) > 0)))

# How many species were found in both years?
base::intersect(species_2015, species_2022)
# Overall, 28 species were found in both years

# How many species were found only in 2022?
base::setdiff(species_2022, species_2015)
# 18 species were found in 2022 but not in 2015

# How many spp were found only in 2015?
base::setdiff(species_2015, species_2022)
# 9 species were found in 2015 but not 2022

# Species accumulation curves #################################################
# Using the raw, unstandardized counts of carabids

# Create subsets of rows for each treatment:
forest_2015 <- carab_by_plot_2015 %>% filter(Treatment == "Forest")
salvaged_2015 <- carab_by_plot_2015 %>% filter(Treatment == "Salvaged")
windthrow_2015 <- carab_by_plot_2015 %>% filter(Treatment == "Windthrow")

forest_2022 <- carab_by_plot_2022 %>% filter(Treatment == "Forest")
salvaged_2022 <- carab_by_plot_2022 %>% filter(Treatment == "Salvaged")
windthrow_2022 <- carab_by_plot_2022 %>% filter(Treatment == "Windthrow")

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

# Taxonomic alpha-diversity: Chao estimator ###################################

carab_by_treatment_2015 <- carab_by_plot_2015 %>% group_by(Treatment) %>%
  summarize(across(all_of(carab_species), ~sum(.)))

carab_by_treatment_2022 <- carab_by_plot_2022 %>% group_by(Treatment) %>%
  summarize(across(all_of(carab_species), ~sum(.)))

forest_counts_2015 <- t(carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Forest", 
                                     carab_species])
windthrow_counts_2015 <- t(carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Windthrow", 
                                        carab_species])
salvaged_counts_2015 <- t(carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Salvaged", 
                                       carab_species])

forest_counts_2022 <- t(carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Forest", 
                                              carab_species])
windthrow_counts_2022 <- t(carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Windthrow", 
                                                 carab_species])
salvaged_counts_2022 <- t(carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Salvaged", 
                                                carab_species])

# Undisturbed forest 2015:
ChaoSpecies(forest_counts_2015, datatype = "abundance", k = 10, conf=0.95)
# 20 observed species
# Chao1 estimates 20.997 species. 95% conf. int: 20.090 -> 31.035

# Windthrow 2015:
ChaoSpecies(windthrow_counts_2015, datatype = "abundance", k = 10, conf=0.95)
# not working - potentially because there aren't any doubletons

# Salvaged 2015:
ChaoSpecies(salvaged_counts_2015, datatype = "abundance", k = 10, conf=0.95)
# 30 observed species
# Chao1 estimates 32.566 species, 95% conf. int: 30.450 -> 44.619

# Undisturbed forest 2022:
ChaoSpecies(forest_counts_2022, datatype = "abundance", k = 10, conf=0.95)
# 30 observed species
# Chao1 estimates 50.195 species, 95% conf. int: 33.964 -> 132.887

# Windthrow 2022:
ChaoSpecies(windthrow_counts_2022, datatype = "abundance", k = 10, conf=0.95)
# 26 observed species
# Chao1 estimates 46.077 species, 95% conf. int: 30.645 -> 112.779

# Salvaged 2022:
ChaoSpecies(salvaged_counts_2022, datatype = "abundance", k = 10, conf=0.95)
# 33 observed species
# Chao1 estimates 65.544 species, 95% conf. int: 41.132 -> 163.240

# Taxonomic alpha-diversity: plot level species richness and Shannon diversity ####

carab_by_plot_2015_stdz$species_richness <- 
  rowSums(carab_by_plot_2015_stdz[,carab_species] > 0)

carab_by_plot_2022_stdz$species_richness <- 
  rowSums(carab_by_plot_2022_stdz[,carab_species] > 0)

# Plot species richness as a function of treatment:
ggplot(data=carab_by_plot_2015_stdz, aes(x=Treatment, y=species_richness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Number of species") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle Species Richness")

ggplot(data=carab_by_plot_2022_stdz, aes(x=Treatment, y=species_richness)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Number of species") + xlab("Forest disturbance") + 
  ggtitle("2022 Ground Beetle Species Richness")

# Hill Numbers
# q = 0 (default) to get species richness, 
# q = 1 to get shannon entropy,
# q = 2 will give inverse Simpson.
# MARGIN = 1 is the default, indicates that sites are rows

# Shannon diversity: exp(-Σp_i*log(p_i))
carab_by_plot_2015_stdz$shannon_diversity <- 
  hill_taxa(carab_by_plot_2015_stdz[,carab_species], q = 1, MARGIN = 1)

carab_by_plot_2022_stdz$shannon_diversity <- 
  hill_taxa(carab_by_plot_2022_stdz[,carab_species], q = 1, MARGIN = 1)

carab_by_plot_2015

# Plot Shannon diversity as a function of treatment:
ggplot(data=carab_by_plot_2015_stdz, aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Shannon diversity") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle Shannon diversity")

ggplot(data=carab_by_plot_2022_stdz, aes(x=Treatment, y=shannon_diversity)) +
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

carab_by_plot_2015_stdz_matrix <- as.matrix(carab_by_plot_2015_stdz[,carab_species])
carab_by_plot_2022_stdz_matrix <- as.matrix(carab_by_plot_2022_stdz[,carab_species])

# Compute the weighted mean_pairwise_distance for each plot
carab_by_plot_2015_stdz$mean_pairwise_distance <-
  picante::mpd(samp = carab_by_plot_2015_stdz_matrix,
               dis = dist,
               abundance.weighted = TRUE)

carab_by_plot_2022_stdz$mean_pairwise_distance <-
  picante::mpd(samp = carab_by_plot_2022_stdz_matrix,
               dis = dist,
               abundance.weighted = TRUE)

# How does mean pairwise distance relate to spp richness?
ggplot(data=carab_by_plot_2015_stdz) + geom_point(aes(x=species_richness,
                                                      y=mean_pairwise_distance))
ggplot(data=carab_by_plot_2022_stdz) + geom_point(aes(x=species_richness,
                                                      y=mean_pairwise_distance))
# The mean_pairwise_distance seems to be positively correlated with spp richness
# for 2015, but not as much in 2022

# How does mean pairwise distance relate to activity-abundance?
ggplot(data=carab_by_plot_2015_stdz) + geom_point(aes(x=total_count_stdz,
                                                      y=mean_pairwise_distance))
ggplot(data=carab_by_plot_2022_stdz) + geom_point(aes(x=total_count_stdz,
                                                      y=mean_pairwise_distance))

# Plot functional alpha-diversity as a function of treatment:
ggplot(data=carab_by_plot_2015_stdz, aes(x=Treatment, y=mean_pairwise_distance)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,0.3)+
  ylab("Mean pairwise distance") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle functional alpha-diversity")

ggplot(data=carab_by_plot_2022_stdz, aes(x=Treatment, y=mean_pairwise_distance)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) + ylim(0,0.3) +
  ylab("Mean pairwise distance") + xlab("Forest disturbance") +
  ggtitle("2022 Ground Beetle functional alpha-diversity")

# Taxonomic beta-diversity ###################################################

# First, I need to convert my data tables to relative abundance data:
carab_by_plot_2015_stdz_rel <- vegan::decostand(carab_by_plot_2015_stdz[carab_species], 
                                                method = "total")

carab_by_plot_2022_stdz_rel <- vegan::decostand(carab_by_plot_2022_stdz[carab_species], 
                                                method = "total")

# Now compute distance matrix between all 24 plots:
dist_spp_space_2015 <- vegdist(carab_by_plot_2015_stdz_rel, method = "bray")
summary(as.vector(dist_spp_space_2015))
# I'm noticing a lot of 1s in this matrix (19 1s), which indicates complete dissimilarity.
# I guess if there is no commonality in the captured species, then the dissimilarity
# is 1.

dist_spp_space_2022 <- vegdist(carab_by_plot_2022_stdz_rel, method = "bray")
summary(as.vector(dist_spp_space_2022))

# Now run the Permutational Multivariate Analysis of Variance:
adonis2(dist_spp_space_2015 ~ carab_by_plot_2015_stdz$Treatment, permutations = 999)

adonis2(dist_spp_space_2022 ~ carab_by_plot_2015_stdz$Treatment, permutations = 999)

# Run the Analysis of Multivariate Homogeneity of Group Dispersions:

beta_dispersion_2015 <- betadisper(d = dist_spp_space_2015, 
                              group = carab_by_plot_2015_stdz$Treatment, 
                              type = c("median"))
beta_dispersion_2015
plot(beta_dispersion_2015)
boxplot(beta_dispersion_2015, ylab = "Distance to median", xlab="")
anova(beta_dispersion_2015)

beta_dispersion_2022 <- betadisper(d = dist_spp_space_2022, 
                                   group = carab_by_plot_2022_stdz$Treatment, 
                                   type = c("median"))
beta_dispersion_2022
plot(beta_dispersion_2022)
boxplot(beta_dispersion_2022, ylab = "Distance to median", xlab="")
anova(beta_dispersion_2022)

# Now create nonmetric multidimensional scaling ordinations (first for 2015):
nmds_2015 <- metaMDS(dist_spp_space_2015, trymax = 500, k = 2)
nmds_2015 # stress is quality of fit
stressplot(nmds_2015)
plot(nmds_2015) 

# plot the 2015 NMDS model by treatment:
ordiplot(nmds_2015, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds_2015, dis = "sites", select = which(carab_by_plot_2015_stdz$Treatment=="Forest"), pch = 15, cex = 1, col = "palegreen4")
points(nmds_2015, dis = "sites", select = which(carab_by_plot_2015_stdz$Treatment=="Salvaged"), pch = 16, cex = 1, col = "brown4")
points(nmds_2015, dis = "sites", select = which(carab_by_plot_2015_stdz$Treatment=="Windthrow"), pch = 17, cex = 1, col = "goldenrod2")
ordiellipse(nmds_2015, carab_by_plot_2022_stdz$Treatment, draw = "lines", col = c("palegreen4", "brown4", "goldenrod2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)
legend("topleft", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("palegreen4", "brown4", "goldenrod2"))

# plot the 2015 NMDS model by area:
ordiplot(nmds_2015, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds_2015, dis = "sites", select = which(carab_by_plot_2015_stdz$Area=="northeast"), pch = 15, cex = 1, col = "blue")
points(nmds_2015, dis = "sites", select = which(carab_by_plot_2015_stdz$Area=="southwest"), pch = 16, cex = 1, col = "purple")
ordiellipse(nmds_2015, carab_by_plot_2022_stdz$Area, draw = "lines", col = c("blue", "purple"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)
legend("topleft", legend = c("northeast", "southwest"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("blue", "purple"))

# Now create a NMDS ordination for 2022:
nmds_2022 <- metaMDS(dist_spp_space_2022, trymax = 500, k = 2)
nmds_2022 # stress is quality of fit
stressplot(nmds_2022)
plot(nmds_2022) 

# plot the 2022 NMDS model by treatment:
ordiplot(nmds_2022, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds_2022, dis = "sites", select = which(carab_by_plot_2022_stdz$Treatment=="Forest"), pch = 15, cex = 1, col = "palegreen4")
points(nmds_2022, dis = "sites", select = which(carab_by_plot_2022_stdz$Treatment=="Salvaged"), pch = 16, cex = 1, col = "brown4")
points(nmds_2022, dis = "sites", select = which(carab_by_plot_2022_stdz$Treatment=="Windthrow"), pch = 17, cex = 1, col = "goldenrod2")
ordiellipse(nmds_2022, carab_by_plot_2022_stdz$Treatment, draw = "lines", col = c("palegreen4", "brown4", "goldenrod2"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)
legend("topleft", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("palegreen4", "brown4", "goldenrod2"))

# plot the 2022 NMDS model by area:
ordiplot(nmds_2022, disp = "sites", type = "n", xlim = c(-1.5, 1.5), ylim = c(-2, 2))
points(nmds_2022, dis = "sites", select = which(carab_by_plot_2022_stdz$Area=="northeast"), pch = 15, cex = 1, col = "blue")
points(nmds_2022, dis = "sites", select = which(carab_by_plot_2022_stdz$Area=="southwest"), pch = 16, cex = 1, col = "purple")
ordiellipse(nmds_2022, carab_by_plot_2022_stdz$Area, draw = "lines", col = c("blue", "purple"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)
legend("topleft", legend = c("northeast", "southwest"),
       pch = c(15, 16, 17), cex = 1, bty = "n", col = c("blue", "purple"))

# Investigate activity-abundance of open-habitat, eurytopic, and ###############
# forest-specialist carabids 

# Import the trait data:
traits <- read.csv("Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv")
rownames(traits) <- traits$Species

open_habitat_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "open habitat", "Species"] # 2 spp
eurytopic_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "eurytopic", "Species"] # 28 spp
forest_specialist_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "forest specialist", "Species"] # 23 spp
unknown_forest_affinity_spp <- traits[is.na(traits$Forest_affinity), "Species"] 
# (2 species have unknown forest affinity)

# Create a column for activity-abundance of open habitat species
carab_by_plot_2015_stdz$open_habitat_spp_stdz <- 
  rowSums(carab_by_plot_2015_stdz[, open_habitat_spp])

# And for eurytopic species:
carab_by_plot_2015_stdz$eurytopic_spp_stdz <- 
  rowSums(carab_by_plot_2015_stdz[, eurytopic_spp])

# And for forest specialists:
carab_by_plot_2015_stdz$forest_specialist_spp_stdz <- 
  rowSums(carab_by_plot_2015_stdz[, forest_specialist_spp])

# Do the same for the 2022 data:
carab_by_plot_2022_stdz$open_habitat_spp_stdz <- 
  rowSums(carab_by_plot_2022_stdz[, open_habitat_spp])

carab_by_plot_2022_stdz$eurytopic_spp_stdz <- 
  rowSums(carab_by_plot_2022_stdz[, eurytopic_spp])

carab_by_plot_2022_stdz$forest_specialist_spp_stdz <- 
  rowSums(carab_by_plot_2022_stdz[, forest_specialist_spp])

# Community-weighted mean trait values #########################################

# I want to calculate community-weighted mean trait values for each
# plot
CWM <- function(community, trait_name){
  geometry::dot(c(t(community[carab_species] /
                      sum(community[carab_species]))),
                traits[,trait_name])
}

geometry::dot(c(t(carab_by_plot_2015_stdz[1, carab_species] / 
                    carab_by_plot_2015_stdz$total_count_stdz[1])), 
              traits$PC1) # weighted mean is simply a dot product of a vector 
# of the relative abundances of each species, with a vector of the trait value 
# for each species

# Test out the function to make sure it works:
CWM(community = carab_by_plot_2015_stdz[1,], 
    trait_name = "PC1")

# I decided to use the function "functcomp" in the package FD to calculate 
# CWMs. It looks like the results match what I got with my function. Here,
# I'll calculate the CWM for each of the first three PC axes, as well
# as Water_affinity and Flight_capability
CWMs_2015 <- FD::functcomp(traits[,c("PC1", "PC2", "PC3", "Water_affinity", "Flight_capability")], 
          as.matrix(carab_by_plot_2015_stdz[, carab_species]), CWM.type = "all")

# Join to the data table:
carab_by_plot_2015_stdz_with_CWMs <- cbind(carab_by_plot_2015_stdz, CWMs_2015)

# And for 2022:
CWMs_2022 <- functcomp(traits[,c("PC1", "PC2", "PC3", "Water_affinity", "Flight_capability")], 
                       as.matrix(carab_by_plot_2022_stdz[, carab_species]), CWM.type = "all")

carab_by_plot_2022_stdz_with_CWMs <- cbind(carab_by_plot_2022_stdz, CWMs_2015)

# Export data tables ###########################################################
# I will run the statistical models in another R script, so I need to export 
# the plot-level data tables

#write.csv(carab_by_plot_2015_stdz_with_CWMs, file="PNR2015_carabid_counts_by_plot_standardized.csv", row.names = F)

# Observation: the capture rate of ground beetles in 2015 varied widely between
# plots. One plot caught a ton of Chlaenius emarginatus and Pterostichus
# moestus

#write.csv(carab_by_plot_2022_stdz_with_CWMs, file="PNR2022_carabid_counts_by_plot_standardized.csv", row.names = F)


