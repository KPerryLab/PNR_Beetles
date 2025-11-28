# Aaron Tayal
# 3/30/2025
# Powdermill ground beetles stats analysis for Aaron's MS Thesis

library(ggplot2)
theme_set(theme_classic(base_size = 14))
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
library(rgl) # For 3 dimensional plots
library(ggpubr) # for the ggarrange function to put multiple graphs in one figure
library(gridGraphics) # for use with ggpubr
library(patchwork) # for lining up figures


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

# Additional species were caught in 2022 in September, and not any other interval
# in 2015 or 2022:
September_only_species <- 
  names(which((colSums(carab_by_interval_2022[7:8,carab_species]) > 0) & 
                (colSums(carab_by_interval_2015[1:6,carab_species]) == 0) &
                (colSums(carab_by_interval_2022[1:6,carab_species]) == 0)))

carab_species_1st_6 <- setdiff(carab_species, September_only_species)

# Make a graph of the captures of Platynus angustatus over the season:
ggplot(data=carab_by_interval_2015, aes(x=Collection_date, y=Platynus_angustatus)) +
  geom_line()
ggplot(data=carab_by_interval_2022, aes(x=Collection_date, y=Platynus_angustatus)) +
  geom_line() + geom_point() + 
  ylab(expression(paste("Number of  ", italic("Platynus angustatus")))) + 
  xlab("Collection date")+
  ggtitle(expression(paste("Captures of  ", italic("Platynus angustatus"), " over the 2022 season")))

# Pterostichus tristis:
ggplot(data=carab_by_interval_2015, aes(x=Set_date, y=Pterostichus_tristis)) +
  geom_line()+ geom_point()
ggplot(data=carab_by_interval_2022, aes(x=Set_date, y=Pterostichus_tristis)) +
  geom_line()+ geom_point()

# Investigate carabid abundance over the season ################################

carab_by_interval_2015$total_count_stdz <- 
  rowSums(carab_by_interval_2015[,carab_species]) / 
  carab_by_interval_2015$num_nonmissing_plots
# Gives the number of carabids caught at an average plot during the ~14 day interval

carab_by_interval_2022$total_count_stdz <- 
  rowSums(carab_by_interval_2022[,carab_species]) / 
  carab_by_interval_2022$num_nonmissing_plots

# The following graph shows how the abundance of ground beetles changed over
# the season for each year
ggplot(data=carab_by_interval_2015, aes(x=Interval, y=total_count_stdz)) +
  geom_line(color="red") + geom_point(color="red") + 
  geom_point(data=carab_by_interval_2022, aes(x=Interval, y=total_count_stdz), color="blue") +
  geom_line(data=carab_by_interval_2022, aes(x=Interval, y=total_count_stdz), color="blue") + ylim(0,12)

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

# How many species were found in both 2015 AND 2022 in the first 6 intervals?
base::intersect(species_2015, species_2022_excluding_last_2_int)


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

# Now remove any columns for species that were only collected in September:
g <- g %>% select(-all_of(September_only_species))

# Now divide the species counts by the approximate number of days active out 
# of 84 days (to account for some trap contents being lost for a few plots)
g_stdz <- g
g_stdz[,carab_species_1st_6] <- g[,carab_species_1st_6] * 84 / g$approx_trap_days

# Investigate total counts of ground beetles ##################################
# How many ground beetles were collected?

g$total_count <- rowSums(g[,carab_species_1st_6])
sum(g %>% filter(Year==2015) %>% select(total_count)) # 934 ground beetles in 2015
sum(g %>% filter(Year==2022) %>% select(total_count)) # 603 ground beetles in 
# the first six intervals of 2022.

g_stdz$total_count_stdz <- rowSums(g_stdz[,carab_species_1st_6]) # totals for the standardized counts
sum(g_stdz %>% filter(Year==2015) %>% select(total_count_stdz)) 
sum(g_stdz %>% filter(Year==2022) %>% select(total_count_stdz)) 


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
forest_2015_accum <- specaccum(comm=forest_2015[carab_species_1st_6], method="random", permutations=100)
salvaged_2015_accum <- specaccum(comm=salvaged_2015[carab_species_1st_6], method="random", permutations=100)
windthrow_2015_accum <- specaccum(comm=windthrow_2015[carab_species_1st_6], method="random", permutations=100)

# Make a graph for 2015:
par(mfrow=c(2,1))
par(mar=c(4,4,2,1))
plot(forest_2015_accum, col = "palegreen3", xvar = c("sites"), lty = 4, 
     lwd = 2, ylab = "Species Richness", xlab = "Number of plots", 
     xlim = c(0, 12), ylim = c(0, 35)) 
plot(salvaged_2015_accum, add=T, xvar = c("sites"), lty = 1, 
     lwd = 2, col = "goldenrod2") 
plot(windthrow_2015_accum, add=T, xvar = c("sites"), lty = 2, 
     lwd = 2, col = "brown4") 
legend("bottomright", legend = c("Forest", "Salvaged", "Windthrow"), 
       lty = c(4,1,2), cex = 0.7, bty = "n", lwd = 2,
       col = c("palegreen3", "goldenrod2", "brown4"))
mtext("A", side = 3, adj = 0, line = 0.5, cex = 1, font = 2)

# Now for 2022:
forest_2022_accum <- specaccum(comm=forest_2022[carab_species_1st_6], method="random", permutations=100)
salvaged_2022_accum <- specaccum(comm=salvaged_2022[carab_species_1st_6], method="random", permutations=100)
windthrow_2022_accum <- specaccum(comm=windthrow_2022[carab_species_1st_6], method="random", permutations=100)

# Make a graph for 2022:
par(mar=c(4,4,2,1))
plot(forest_2022_accum, col = "palegreen3", xvar = c("sites"), lty = 4, 
     lwd = 2, ylab = "Species Richness", xlab = "Number of plots", 
     xlim = c(0, 12), ylim = c(0, 35)) 
plot(salvaged_2022_accum, add=T, xvar = c("sites"), lty = 1, 
     lwd = 2, col = "goldenrod2") 
plot(windthrow_2022_accum, add=T, xvar = c("sites"), lty = 2, 
     lwd = 2, col = "brown4") 
legend("bottomright", legend = c("Forest", "Salvaged", "Windthrow"), 
       lty = c(4,1,2), cex = 0.7, bty = "n", lwd = 2,
       col = c("palegreen3", "goldenrod2", "brown4"))
mtext("B", side = 3, adj = 0, line = 0.5, cex = 1, font = 2)

# Taxonomic alpha-diversity: Chao estimators ###################################
# still omitting the last two intervals of the 2022 data

total_counts_2015 <- colSums(g %>% filter(Year==2015) %>% select(all_of(carab_species_1st_6)))
total_counts_2022 <- colSums(g %>% filter(Year==2022) %>% select(all_of(carab_species_1st_6)))

# Total counts 2015:
ChaoSpecies(total_counts_2015, datatype = "abundance", k = 10, conf=0.95)
# 37 observed species
# Chao1 estimates 47.114 species. 95% conf. int: 39.262 -> 82.233

# Total counts 2022:
ChaoSpecies(total_counts_2022, datatype = "abundance", k = 10, conf=0.95)
# 37 observed species
# Chao1 estimates 79.180 species. 95% conf. int: 46.183 -> 230.745


# Taxonomic alpha-diversity: plot level species richness and Shannon diversity ####

g_stdz$sp_rich <- rowSums(g_stdz[,carab_species_1st_6] > 0)

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
  hill_taxa(g_stdz[,carab_species_1st_6], q = 1, MARGIN = 1)

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

dist_0 <- read.csv("Aaron_PNR_formatted_data/carabid_dist_in_trait_space.csv") # Gower distances

# Check if the species column names match between the distance matrix and the
# species abundance data:
all.equal(colnames(dist_0), carab_species_1st_6) # True

dist <- as.matrix(dist_0)
rownames(dist) <- colnames(dist)

g_stdz_matrix <- as.matrix(g_stdz[,carab_species_1st_6])

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

treatment_colors = c("Forest" = "palegreen3", "Salvaged" = "goldenrod2", "Windthrow" = "brown4")
g_stdz$Treatment <- factor(g_stdz$Treatment, levels = c("Windthrow", "Salvaged", "Forest"))

# First, I need to convert my data table to relative abundance data:
g_stdz_rel <- vegan::decostand(g_stdz[carab_species_1st_6], method = "total")

# Now compute distance matrix between all 48 plots:
dist_spp_space <- vegdist(g_stdz_rel, method = "bray")
summary(as.vector(dist_spp_space))
hist(as.vector(dist_spp_space))
# A distance of 1 indicates complete dissimilarity.
# I guess if there is no commonality in the captured species between two plots, 
# then the dissimilarity is 1.

# Now create nonmetric multidimensional scaling ordinations: # TRYING 3 dimensions:
nmds3d <- metaMDS(dist_spp_space, trymax = 500, k = 3)
nmds3d # stress is quality of fit
stressplot(nmds3d)
g_stdz$NMDS3d1 <- nmds3d$points[,1]
g_stdz$NMDS3d2 <- nmds3d$points[,2]
g_stdz$NMDS3d3 <- nmds3d$points[,3]
#plot3d(g_stdz$NMDS3d1, g_stdz$NMDS3d2, g_stdz$NMDS3d3, col = treatment_colors)
#plot3d(g_stdz$NMDS3d1, g_stdz$NMDS3d2, g_stdz$NMDS3d3, col = g_stdz$Year)
area_colors = c("southwest"="green", "northeast"="blue")
#plot3d(g_stdz$NMDS3d1, g_stdz$NMDS3d2, g_stdz$NMDS3d3, col = area_colors)

# Now create nonmetric multidimensional scaling ordination with 2 dimensions:
nmds <- metaMDS(dist_spp_space, trymax = 500, k = 2)
nmds # stress is quality of fit
stressplot(nmds)
g_stdz$NMDS1 <- nmds$points[,1]
g_stdz$NMDS2 <- nmds$points[,2]

#Try a plot using the method suggested by ChatGPT (hopefully it adds ellipses)
plot(nmds, type = "n")
points(nmds, display = "sites", col = g_stdz$Treatment, pch = 19)
ordiellipse(nmds, group=g_stdz$Treatment, kind = "se", conf = 0.95, label = TRUE)

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Year)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Area)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Area, shape = Year)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Transect)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment)) + 
  geom_point(size=2) + scale_color_manual(values = treatment_colors)+ coord_fixed() +
  stat_ellipse(level = 0.95)

ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment, shape=Year)) + 
  geom_point(size=2) + scale_color_manual(values = treatment_colors) + coord_fixed() +
  stat_ellipse(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment, group = Treatment), level = 0.95)

taxonomic_NMDS <- ggplot(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment, shape=Year, group=Plot)) + 
  geom_point(size=2) + scale_color_manual(values = treatment_colors) +
  coord_fixed() + scale_x_continuous(limits = c(-2,2)) + scale_y_continuous(limits = c(-1.8,1.8)) +
  theme(axis.text = element_blank(),
        legend.background = element_rect(color = "black", linewidth = 0.2),
        legend.box.margin = margin(5, 5, 5, 20),
        plot.margin = margin(10,10,10,10)) +
  stat_ellipse(data=g_stdz, aes(x=NMDS1, y=NMDS2, color=Treatment, group = Treatment), level = 0.95)

taxonomic_NMDS

# Now run the Permutational Multivariate Analysis of Variance (PERMANOVA):

adonis2(dist_spp_space ~ Treatment * Year, data=g_stdz,
        permutations = 999, by="terms") # May need to change this to by="margin"

# Because the p-value of the adonis2 function is changing each time I run
# the line of code, I'll need to increase the number of permutations so that 
# I get a more consistent p-value:
#adonis2(dist_spp_space ~ Treatment * Year, data=g_stdz,
#        permutations = 99999, by="terms")

# I was initially only going to test for differences between salvaged and windthrow:
# OUTDATED: Make species abundance matrix for only the windthrow and salvaged plots:
#g_stdz_ws <- g_stdz %>% filter(Treatment != "Forest") # ws stands for windthrow salvaged
#g_stdz_rel_ws <- vegan::decostand(g_stdz_ws[carab_species], method = "total")
# Make a distance matrix between plots using this subset of plots:
#dist_spp_space_ws <- vegdist(g_stdz_rel_ws, method = "bray")

# JUST TO TRY IT OUT: Also, try a PERMANOVA with area (northeast blowdown or southwest blowdown):
adonis2(dist_spp_space ~ Area, data=g_stdz, permutations = 999, by="terms")

# Run the pairwise adonis test:
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
                                     dis = dist, # the Gower distance between each carabid species in trait space
                                     abundance.weighted = T) # mean-pairwise distances will be weighted by species abundances
hist(as.vector(dist_functional_beta))

# Now create nonmetric multidimensional scaling ordinations: # TRYING 3 dimensions:
nmds_functional_3d <- metaMDS(dist_functional_beta, trymax = 500, k = 3)
nmds_functional_3d # stress is quality of fit
stressplot(nmds_functional_3d)
g_stdz$NMDS_functional_3d1 <- nmds_functional_3d$points[,1]
g_stdz$NMDS_functional_3d2 <- nmds_functional_3d$points[,2]
g_stdz$NMDS_functional_3d3 <- nmds_functional_3d$points[,3]
#plot3d(g_stdz$NMDS_functional_3d1, g_stdz$NMDS_functional_3d2, g_stdz$NMDS_functional_3d3, 
#       col = treatment_colors)
#plot3d(g_stdz$NMDS_functional_3d1, g_stdz$NMDS_functional_3d2, g_stdz$NMDS_functional_3d3, 
#       col = g_stdz$Year)
#plot3d(g_stdz$NMDS_functional_3d1, g_stdz$NMDS_functional_3d2, g_stdz$NMDS_functional_3d3, 
#       col = area_colors)

# Now create an NMDS ordination with 2 dimensions:
nmds_functional <- metaMDS(dist_functional_beta, trymax = 500, k = 2)
nmds_functional # stress is quality of fit
stressplot(nmds_functional)

g_stdz$NMDS_functional1 <- nmds_functional$points[,1]
g_stdz$NMDS_functional2 <- nmds_functional$points[,2]

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Year)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Area)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Area, shape=Year)) + 
  geom_point(size=2)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Treatment)) + 
  geom_point(size=2) + scale_color_manual(values = treatment_colors)+ coord_fixed()

ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Treatment,
                        shape=Year, linetype = Year)) + coord_fixed() + theme(legend.background = element_rect(color = "black",
                                                                                              linewidth = 0.2),
                                                             legend.box.margin = margin(5, 5, 5, 20)) +
  geom_point(size=2) + scale_color_manual(values = treatment_colors) +
  stat_ellipse(level = 0.95)

functional_NMDS <- ggplot(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, color=Treatment,
                        shape=Year, group=Plot)) + 
  geom_point(size=2) + scale_color_manual(values = treatment_colors) +
  xlab("NMDS1") + ylab("NMDS2") + coord_fixed() +
  scale_x_continuous(limits = c(-0.18,0.18)) + scale_y_continuous(limits = c(-0.18,0.18)) +
  theme(axis.text = element_blank(),
        legend.background = element_rect(color = "black", linewidth = 0.2),
        legend.box.margin = margin(5, 5, 5, 20),
        plot.margin = margin(10,10,10,10)) +
  stat_ellipse(data=g_stdz, aes(x=NMDS_functional1, y=NMDS_functional2, 
                                color=Treatment, group = Treatment), level = 0.95)

functional_NMDS

ggarrange(taxonomic_NMDS, functional_NMDS,
          labels = c("A", "B"), ncol=1, nrow=2)

(taxonomic_NMDS / functional_NMDS) + plot_layout(heights = c(1, 1))

# Now run a PERMANOVA to test the null hypothesis that the centroids of each 
# treatment group are identical and their dispersions are identical:

# Make a distance matrix between plots for just the salvaged and windthrow plots:
#g_stdz_matrix_ws <- as.matrix(g_stdz_ws[,carab_species])
#dist_functional_beta_ws <- comdist(comm = g_stdz_matrix_ws, dis = dist,
#                                abundance.weighted = T)

#adonis2(dist_functional_beta ~ g_stdz$Treatment * g_stdz$Year,
#         permutations = 999, by="terms") 

adonis2(dist_functional_beta ~ g_stdz$Treatment * g_stdz$Year,
        permutations = 99999, by="terms") 

# Functional beta-diversity for each individual year: ##########################
# Because there was a significant interaction term, I need to interpret the 
# interaction, so I need to run separate PERMANOVAs for each year:

g_stdz_matrix_2015 <- as.matrix(g_stdz[g_stdz$Year == 2015 ,carab_species_1st_6])
dist_functional_beta_2015 <- comdist(comm = g_stdz_matrix_2015, # the activity-abundance data for each plot
                                dis = dist, # the Gower distance between each carabid species in trait space
                                abundance.weighted = T)
adonis2(dist_functional_beta_2015 ~ Treatment, data = g_stdz[g_stdz$Year == "2015",])
pairwise.adonis(dist_functional_beta_2015, g_stdz[g_stdz$Year == "2015","Treatment"])

g_stdz_matrix_2022 <- as.matrix(g_stdz[g_stdz$Year == 2022 ,carab_species_1st_6])
dist_functional_beta_2022 <- comdist(comm = g_stdz_matrix_2022, # the activity-abundance data for each plot
                                     dis = dist, # the Gower distance between each carabid species in trait space
                                     abundance.weighted = T)
adonis2(dist_functional_beta_2022 ~ Treatment, data = g_stdz[g_stdz$Year == "2022",])

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

# Mantel tests ################################################################

# To investigate geographic distance vs. distance in spp space 

# A Mantel test looks at correlations between distance matrices.

# Make a taxonomic distance matrix just for the 2015 data:
g_stdz_2015 <- g_stdz %>% filter(Year==2015)
g_stdz_rel_2015 <- vegan::decostand(g_stdz_2015[carab_species_1st_6], method = "total")
dist_spp_space_2015 <- vegdist(g_stdz_rel_2015, method = "bray")

# Make a taxonomic distance matrix just for the 2022 data:
g_stdz_2022 <- g_stdz %>% filter(Year==2022)
g_stdz_rel_2022 <- vegan::decostand(g_stdz_2022[carab_species_1st_6], method = "total")
dist_spp_space_2022 <- vegdist(g_stdz_rel_2022, method = "bray")

# Import the distance between plots in meters:
dist_geographic <- as.dist(read.csv("newer_maps/Powdermill_dist_matrix.csv") %>% select(-Plot))

# Examine the relationship between distance in geographic space and 
# distance in species space:
mantel(dist_spp_space_2015, dist_geographic)
mantel(dist_spp_space_2022, dist_geographic)

plot(as.vector(dist_spp_space_2015), as.vector(dist_geographic))
plot(as.vector(dist_spp_space_2022), as.vector(dist_geographic))


# Investigate activity-abundance of open-habitat, eurytopic, and ###############
# forest-specialist carabids 

# Import the trait data:
traits <- read.csv("Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv")
rownames(traits) <- traits$Species

# Of the 47 species caught in the first 6 collecion intervals, how many were 
# open-habitat and eurytopic?
table(traits %>% select(Forest_affinity))

# make lists of the eurytopic, open-habitat, and forest specialist species
open_habitat_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "open habitat", "Species"] # 2 spp
eurytopic_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "eurytopic", "Species"] # 24 spp
forest_specialist_spp <- traits[!is.na(traits$Forest_affinity) &
  traits$Forest_affinity == "forest specialist", "Species"] # 20 spp
unknown_forest_affinity_spp <- traits[is.na(traits$Forest_affinity), "Species"] 
# (1 species has unknown forest affinity)

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
  geometry::dot(c(t(community[carab_species_1st_6] /
                      sum(community[carab_species_1st_6]))),
                traits[,trait_name])
}

geometry::dot(c(t(g_stdz[1, carab_species_1st_6] / 
                    g_stdz$total_count_stdz[1])), 
              traits$PC1) # weighted mean is simply a dot product of a vector 
# of the relative abundances of each species, with a vector of the trait value 
# for each species

# Test out the function to make sure it works:
CWM(community = g_stdz[1,], trait_name = "PC1")

# I decided to use the function "functcomp" in the package FD to calculate 
# CWMs (more repeatable). # It looks like the results match what I got with my function. Here,
# I'll calculate the CWM for each of the first three PC axes, as well
# as the eight numerical traits and Water_affinity and Flight_capability
# and standardized antenna length
trait_list_standard <- c("body_length", "antenna_length_standard", "eye_protrusion_standard",
                         "eye_length_standard", "pronotum_width_standard",
                         "abdomen_width_standard", "rear_leg_length_standard", 
                         "rear_trochanter_length_standard", "Water_affinity", 
                         "Flight_capability")
CWMs <- FD::functcomp(traits[,c("PC1", "PC2", "PC3", trait_list_standard)], 
          as.matrix(g_stdz[, carab_species_1st_6]), CWM.type = "all")

# Join to the data table:
g_stdz_with_CWMs <- cbind(g_stdz, CWMs)

# Counts of ground beetles data table ##########################################
# I want to make a basic data table that shows the counts of all the different
# species of ground beetle in different management areas and in different years.

# Order the species according to Bousquet 2012:
ordered_spp <- c("Notiophilus_aeneus", "Sphaeroderus_canadensis", "Sphaeroderus_stenostomus",
                 "Scaphinotus_viduus", "Scaphinotus_imperfectus", "Carabus_goryi",
                 "Lophoglossus_scrutator", "Pterostichus_mutus", "Pterostichus_corvinus",
                 "Pterostichus_sayanus", "Pterostichus_coracinus", "Pterostichus_melanarius",
                 "Pterostichus_lachrymosus", "Pterostichus_stygicus", "Pterostichus_hamiltoni",
                 "Pterostichus_moestus", "Pterostichus_diligendus", "Pterostichus_rostratus",
                 "Pterostichus_adoxus", "Pterostichus_tristis", "Cyclotrachelus_fucatus",
                 "Cyclotrachelus_convivus", "Cyclotrachelus_sigillatus", "Chlaenius_emarginatus",
                 "Chlaenius_laticollis", "Dicaelus_politus", "Dicaelus_teter",
                 "Notiobia_nitidipennis", "Anisodactylus_harrisii", "Anisodactylus_melanopus",
                 "Anisodactylus_nigerrimus", "Amphasia_interstitialis", "Agonoleptus_thoracicus",
                 "Harpalus_spadiceus", "Trichotichnus_autumnalis", "Pseudamara_arenaria",
                 "Olisthopus_parmatus", "Agonum_ferreum", "Agonum_fidele",
                 "Agonum_retractum", "Platynus_decentis", "Platynus_tenuicollis",
                 "Platynus_angustatus", "Cymindis_limbata", "Cymindis_platicollis",
                 "Apenes_lucidula", "Galerita_bicolor")

# Make a table of the counts for each species and each year:

count_table <- data.frame(species = ordered_spp,
                          counts_2015_w = unname(colSums(g %>% filter(Year==2015, Treatment=="Windthrow") %>% select(all_of(ordered_spp)))),
                          counts_2015_s = unname(colSums(salvaged_2015[,ordered_spp])),
                          counts_2015_f = unname(colSums(forest_2015[,ordered_spp])),
                          counts_2022_excluding_last_2_int_w = unname(colSums(windthrow_2022[,ordered_spp])),
                          counts_2022_excluding_last_2_int_s = unname(colSums(salvaged_2022[,ordered_spp])),
                          counts_2022_excluding_last_2_int_f = unname(colSums(forest_2022[,ordered_spp])))

#write.csv(count_table, "Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts_summary.csv")

# Export data tables ###########################################################
# I will run the statistical models in another R script, so I need to export 
# the plot-level data tables

#write.csv(g_stdz_with_CWMs, file="Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts_by_plot_standardized.csv", row.names = F)

# Observation: the capture rate of ground beetles in 2015 varied widely between
# plots. One plot caught a ton of Chlaenius emarginatus and Pterostichus
# moestus


