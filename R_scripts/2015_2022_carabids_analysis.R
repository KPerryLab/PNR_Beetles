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

# Data standardization #########################################################

# First, calculate totals across all sampling intervals, for each year:

carab_species <- colnames(carab_0)[7:63] # IMPORTANT: make sure the species columns
# start at column 7

# For 2015:
carab_by_plot_2015_0 <- carab_0 %>% filter(Year==2015) %>% group_by(Plot) %>%
  summarise(num_nonmissing_intervals = n() - sum(is.na(Agonoleptus_thoracicus)),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))
# The column num_nonmissing intervals is the number of intervals, for a given
# plot, which have carabid data (the trap contents were not lost)

# Now join the 2015 carabid count data to the plot information dataframe:
carab_by_plot_2015_1 <- full_join(plot_locations, carab_by_plot_2015_0, by="Plot")

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
carab_by_plot_2022 <- tibble(full_join(plot_locations, carab_by_plot_2022_1, 
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

carab_by_plot_2022_stdz$total_count_stdz <- rowSums(carab_by_plot_2022_stdz[,carab_species])
mean(carab_by_plot_2022_stdz$total_count_stdz) # Each trap, on average, caught
# about 0.33 ground beetles per day in 2022

# Investigate number of species found #########################################

# Make a list of species found in each year:
species_2015 <- names(which((colSums(carab_by_interval_2015[,carab_species]) > 0)))
species_2022 <- names(which((colSums(carab_by_interval_2022[,carab_species]) > 0)))

# How many species were found in both years?
base::intersect(species_2015, species_2022)

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

# Chao estimator ##############################################################

carab_by_treatment_2015 <- carab_by_plot_2015 %>% group_by(Treatment) %>%
  summarize(across(all_of(carab_species), ~sum(.)))

carab_by_treatment_2022 <- carab_by_plot_2022 %>% group_by(Treatment) %>%
  summarize(across(all_of(carab_species), ~sum(.)))

forest_counts_2015 <- t(carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Forest", 
                                     carab_species])
windthrow_counts_2015 <- carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Windthrow", 
                                        carab_species]
salvaged_counts_2015 <- carab_by_treatment_2015[carab_by_treatment_2015$Treatment=="Salvaged", 
                                       carab_species]

forest_counts_2022 <- carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Forest", 
                                              carab_species]
windthrow_counts_2022 <- carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Windthrow", 
                                                 carab_species]
salvaged_counts_2022 <- carab_by_treatment_2022[carab_by_treatment_2022$Treatment=="Salvaged", 
                                                carab_species]

# Undisturbed forest 2015:
ChaoSpecies(forest_counts_2015, datatype = "abundance", k = 10, conf=0.95)

# Windthrow 2015:
ChaoSpecies(windthrow_counts_2015, datatype = "abundance", k = 10, conf=0.95)

# Salvaged 2015:
ChaoSpecies(salvaged_counts_2015, datatype = "abundance", k = 10, conf=0.95)

# Undisturbed forest 2022:
ChaoSpecies(forest_counts_2022, datatype = "abundance", k = 10, conf=0.95)

# Windthrow 2022:
ChaoSpecies(windthrow_counts_2022, datatype = "abundance", k = 10, conf=0.95)

# Salvaged 2022:
ChaoSpecies(salvaged_counts_2022, datatype = "abundance", k = 10, conf=0.95)

# Investigate accumulation of ground beetle species over the season ############

carab_by_interval_2015 <- carab_0 %>% filter(Year==2015) %>% group_by(Interval) %>%
  summarise(num_nonmissing_plots = n() - sum(is.na(Agonoleptus_thoracicus)),
            Set_date = first(Set_date), Collection_date = first(Collection_date),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

carab_by_interval_2022 <- carab_0 %>% filter(Year==2022) %>% group_by(Interval) %>%
  summarise(num_nonmissing_plots = n() - sum(is.na(Agonoleptus_thoracicus)),
            Set_date = first(Set_date), Collection_date = first(Collection_date),
            across(all_of(carab_species), ~ sum(.x, na.rm=TRUE)))

plot(specaccum(comm=carab_by_interval_2015[,carab_species], method="collector"),
     xlab="Interval") + title("2015 accumulation of species over the season") 

plot(specaccum(comm=carab_by_interval_2022[,carab_species], method="collector"),
     xlab="Interval") + title("2022 accumulation of species over the season") 

# Additional species were caught in 2022 in September:

September_only_species <- names(which((colSums(carab_by_interval_2022[7:8,carab_species]) > 0) & 
  (colSums(carab_by_interval_2022[1:6,carab_species]) == 0)))

carab_by_interval_2015[, September_only_species]
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

# Taxonomic alpha-diversity ####################################################

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

# Plot Shannon diversity as a function of treatment:
ggplot(data=carab_by_plot_2015_stdz, aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Shannon diversity") + xlab("Forest disturbance") + 
  ggtitle("2015 Ground Beetle Shannon diversity")

ggplot(data=carab_by_plot_2022_stdz, aes(x=Treatment, y=shannon_diversity)) +
  geom_jitter(width=0.05, height=0, alpha=0.5) +  ylim(0,15)+
  ylab("Shannon diversity") + xlab("Forest disturbance") + 
  ggtitle("2022 Ground Beetle Shannon diversity")





