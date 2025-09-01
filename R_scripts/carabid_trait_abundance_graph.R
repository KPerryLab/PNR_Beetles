# Aaron Tayal
# July 17, 2025
# Trait abundance graph

# I want to display the species found at a particular forest management treatment,
# plotted onto the first two PC axes. The size of the dot will indicate the 
# activity-abundance of that species.

library(ggplot2)
library(dplyr)
library(ggforce) # for the geom_arc_bar() command
library(tidyr)

# Read in the community data:
c2015 <- read.csv("Aaron_PNR_formatted_data/PNR2015_carabid_counts_by_plot_standardized.csv")
c2022 <- read.csv("Aaron_PNR_formatted_data/PNR2022_carabid_counts_by_plot_standardized.csv")

# Read in the trait data:
traits <- read.csv("Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv")

spp <- traits$Species
traits$species_short <- sub(".*_", "", traits$Species)

# Compute the total activity-abundance of each species (for each year)
traits$abundance_2015 <- colSums(c2015[,spp])
traits$abundance_2022 <- colSums(c2022[,spp])

# Compute the activity-abundance for each species, in each treatment:
traits$abundance_forest_2015 <- colSums(c2015[c2015$Treatment=="Forest", spp])
traits$abundance_windthrow_2015 <- colSums(c2015[c2015$Treatment=="Windthrow", spp])
traits$abundance_salvaged_2015 <- colSums(c2015[c2015$Treatment=="Salvaged", spp])

# For each species, compute the proportion that were caught in each treatment:
traits$rel_ab_2015_Forest <- traits$abundance_forest_2015 / traits$abundance_2015
traits$rel_ab_2015_Windthrow <- traits$abundance_windthrow_2015 / traits$abundance_2015
traits$rel_ab_2015_Salvaged <- traits$abundance_salvaged_2015 / traits$abundance_2015

# Same for 2022:
traits$abundance_forest_2022 <- colSums(c2022[c2022$Treatment=="Forest", spp])
traits$abundance_windthrow_2022 <- colSums(c2022[c2022$Treatment=="Windthrow", spp])
traits$abundance_salvaged_2022 <- colSums(c2022[c2022$Treatment=="Salvaged", spp])

# For each species, compute the proportion that were caught in each treatment:
traits$rel_ab_2022_Forest <- traits$abundance_forest_2022 / traits$abundance_2022
traits$rel_ab_2022_Windthrow <- traits$abundance_windthrow_2022 / traits$abundance_2022
traits$rel_ab_2022_Salvaged <- traits$abundance_salvaged_2022 / traits$abundance_2022


# Pivot the data from a "wide" to a "long" format. I need a single column that 
# has a relative abundance (for example, the abundance of Carabus goryi in
# the windthrow divided by the total abundance of Carabus goryi). Thus for each 
# species there will be three rows:
traits_long_2015 <- pivot_longer(traits, cols = c("rel_ab_2015_Forest",
                                                  "rel_ab_2015_Windthrow",
                                                  "rel_ab_2015_Salvaged"),
                                 names_to = "Treatment",
                                 values_to = "relative_abundance_2015",
                                 names_prefix = "rel_ab_2015_")

traits_long_2022 <- pivot_longer(traits, cols = c("rel_ab_2022_Forest",
                                                  "rel_ab_2022_Windthrow",
                                                  "rel_ab_2022_Salvaged"),
                                 names_to = "Treatment",
                                 values_to = "relative_abundance_2022",
                                 names_prefix = "rel_ab_2022_")

# Plot a graphs:

ggplot(traits_long_2015 %>% filter(Species != "Notiophilus_aeneus")) + geom_arc_bar(aes(x0 = PC1, y0 = PC2, r0 = 0, 
                                            r = sqrt(abundance_2015)/3, 
                                            amount = relative_abundance_2015, 
                                            fill = Treatment), stat = "pie", alpha=0.4) +
  coord_fixed() + theme_classic() + xlab("PC1") + ylab("PC2") +
  scale_y_continuous(breaks=seq(-3, 3, 1)) +
  scale_fill_manual(values = c("Forest" = "palegreen4", "Salvaged" = "brown4", 
                               "Windthrow" = "goldenrod2"))+
  geom_text(data=traits %>% filter(Species != "Notiophilus_aeneus"), aes(x=PC1,y=PC2,label=species_short, size=abundance_2015))

ggplot(traits_long_2022 %>% filter(Species != "Notiophilus_aeneus")) + geom_arc_bar(aes(x0 = PC1, y0 = PC2, r0 = 0, 
                                                                                        r = sqrt(abundance_2022)/3, 
                                                                                        amount = relative_abundance_2022, 
                                                                                        fill = Treatment), stat = "pie", alpha=0.4) +
  coord_fixed() + theme_classic() + xlab("PC1") + ylab("PC2") +
  scale_y_continuous(breaks=seq(-3, 3, 1)) +
  scale_fill_manual(values = c("Forest" = "palegreen4", "Salvaged" = "brown4", 
                               "Windthrow" = "goldenrod2"))+
  geom_text(data=traits %>% filter(Species != "Notiophilus_aeneus"), aes(x=PC1,y=PC2,label=species_short, size=abundance_2022))


ggplot(data = traits_long_2015 %>% filter(Species != "Notiophilus_aeneus")) + 
  geom_arc_bar(aes(x0 = rear_trochanter_length_standard, 
                                            y0 = eye_length_standard, r0 = 0, 
                                            r = sqrt(abundance_2015)/300, 
                                            amount = relative_abundance_2015, 
                                            fill = Treatment), stat = "pie", alpha=0.4) +
  coord_fixed() + theme_classic() + xlab("Standardized rear trochanter length") + ylab("Standardized eye length") +
  scale_y_continuous(breaks=seq(-3, 3, 1)) +
  scale_fill_manual(values = c("Forest" = "palegreen4", "Salvaged" = "brown4", 
                               "Windthrow" = "goldenrod2"))+geom_text(data=traits%>% filter(Species != "Notiophilus_aeneus"), aes(x=rear_trochanter_length_standard,y=eye_length_standard,label=Species, alpha=abundance_2015/2.3264))

ggplot(data = traits_long_2022 %>% filter(Species != "Notiophilus_aeneus")) + 
  geom_arc_bar(aes(x0 = rear_trochanter_length_standard, 
                   y0 = eye_length_standard, r0 = 0, 
                   r = sqrt(abundance_2022)/300, 
                   amount = relative_abundance_2022, 
                   fill = Treatment), stat = "pie", alpha=0.4) +
  coord_fixed() + theme_classic() + xlab("Standardized rear trochanter length") + ylab("Standardized eye length") +
  scale_y_continuous(breaks=seq(-3, 3, 1)) +
  scale_fill_manual(values = c("Forest" = "palegreen4", "Salvaged" = "brown4", 
                               "Windthrow" = "goldenrod2"))+geom_text(data=traits%>% filter(Species != "Notiophilus_aeneus"), aes(x=rear_trochanter_length_standard,y=eye_length_standard,label=Species, alpha=abundance_2022/2.3264))


# Graphs of flight capable vs. flight incapable beetles ########################



ggplot(data=traits %>% filter(Species != "Notiophilus_aeneus"),
       aes(x=PC1, y=PC2, size=abundance_2015,
           color=Flight_capability)) + geom_point(alpha=0.5) + 
  geom_text(aes(label = species_short))

ggplot(data=traits %>% filter(Species != "Notiophilus_aeneus"),
       aes(x=PC1, y=PC2, size=abundance_2022,
           color=Flight_capability)) + geom_point(alpha=0.5) + 
  geom_text(aes(label = species_short))













