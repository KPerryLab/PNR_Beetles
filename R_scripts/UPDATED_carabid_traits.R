# Aaron Tayal
# UPDATED VERSION, Oct 3, 2025

# Purpose is to:
# 1. Investigate correlation between the trait values
# 2. Divide traits by body length to standardize them
# 3. Run a principal components analysis to extract PC axes that explain 95% 
# of the variation
# 4. Create a distance matrix between species in trait space

library(tidyverse)
library(readxl)
library(ggplot2)
theme_set(theme_classic())
library(corrplot) # creates correlation visualizations
library(factoextra) # plot PCA results
library(FD) # for gower distance
library(geometry) # for dot product

traits_0 <- read_excel("PNR_Raw_Data/PNR_SpeciesTraits_2022.xlsx", sheet=2, na='NA')

# Which species were found between June and August? ############################

# We are removing species found in September 2022

species_2015_2022_1st_6 <- c("Agonum_ferreum", "Agonum_fidele", "Agonum_retractum", "Amphasia_interstitialis", 
                             "Anisodactylus_harrisii", "Anisodactylus_melanopus", "Anisodactylus_nigerrimus", 
                             "Carabus_goryi", "Chlaenius_emarginatus", "Chlaenius_laticollis", 
                             "Cyclotrachelus_convivus", "Cyclotrachelus_fucatus", "Cyclotrachelus_sigillatus", 
                             "Dicaelus_politus", "Dicaelus_teter", "Harpalus_spadiceus", "Notiobia_nitidipennis", 
                             "Notiophilus_aeneus", "Olisthopus_parmatus", "Platynus_angustatus", 
                             "Platynus_tenuicollis", "Pseudamara_arenaria", "Pterostichus_adoxus", 
                             "Pterostichus_coracinus", "Pterostichus_corvinus", "Pterostichus_diligendus", 
                             "Pterostichus_lachrymosus", "Pterostichus_melanarius", "Pterostichus_moestus", 
                             "Pterostichus_mutus", "Pterostichus_rostratus", "Pterostichus_stygicus", 
                             "Pterostichus_tristis", "Scaphinotus_viduus", "Sphaeroderus_canadensis", 
                             "Sphaeroderus_stenostomus", "Trichotichnus_autumnalis", "Agonoleptus_thoracicus", 
                             "Apenes_lucidula", "Cymindis_limbata", "Cymindis_platicollis", 
                             "Galerita_bicolor", "Lophoglossus_scrutator", "Platynus_decentis", 
                             "Pterostichus_hamiltoni", "Pterostichus_sayanus", "Scaphinotus_imperfectus")

traits <- traits_0 %>% filter(Species %in% species_2015_2022_1st_6)

unique(traits_0$Species)
unique(traits$Species)
setdiff(unique(traits_0$Species), unique(traits$Species)) # 8 species will get removed

traits$Rear_tarsi_length <- as.numeric(traits$Rear_tarsi_length)
# one carabid has missing rear tarsi
traits$Forest_affinity <- as.factor(traits$Forest_affinity)
traits$Water_affinity <- as.numeric(traits$Water_affinity)
traits$Flight_capability <- as.numeric(traits$Flight_capability)
traits$Genus <- as.factor(traits$Genus)
traits$Location_of_collection <- as.factor(traits$Location_of_collection)
table(traits$Location_of_collection) # Most of the traits were measured from beetles
# collected at Powdermill, but a few were measured from beetles collected in 
# Erie and Cuyahoga counties, OH
traits$Species <- as.factor(traits$Species)


# Make an abbreviation for Genus and species: ##################################

# Function to abbreviate
abbrev_species <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  genus <- parts[1]
  species <- parts[2]
  paste0(substr(genus, 1, 2), ".", substr(species, 1, 2))
}

traits$Species_char <- as.character(traits$Species)
traits <- traits %>% mutate(spp_abbrev = sapply(Species_char, abbrev_species))

# Calculate some traits that are sums or differences of measurements ###########

# Body length:
traits$body_length <- traits$Elytra_length + traits$Pronotum_length + traits$Head_length

# Eye protrusion
traits$eye_protrusion <- (traits$Outer_eye_distance - traits$Inner_eye_distance)

# Rear leg length
traits$rear_leg_length <- 
  ( traits$Rear_femur_length + traits$Rear_tibia_length + 
      traits$Rear_tarsi_length )

# Body length standardization ##################################################
# Now I'd like to calculate body length adjusted values for each specific trait:

traits$antenna_length_standard <- traits$Antenna_length / traits$body_length

traits$eye_protrusion_standard <- 
  (traits$eye_protrusion) / traits$body_length

traits$eye_length_standard <- traits$Eye_length / traits$body_length

traits$pronotum_width_standard <- traits$Pronotum_width / traits$body_length

traits$abdomen_width_standard <- traits$Abdomen_width / traits$body_length

traits$rear_leg_length_standard <- 
  ( traits$rear_leg_length) / traits$body_length

traits$rear_trochanter_length_standard <- traits$Rear_trochanter_length / 
  traits$body_length

trait_list <- c("body_length", "Antenna_length", "eye_protrusion",
            "Eye_length", "Pronotum_width",
            "Abdomen_width", "rear_leg_length", 
            "Rear_trochanter_length", "Water_affinity", 
            "Flight_capability")
trait_list_standard <- c("body_length", "antenna_length_standard", "eye_protrusion_standard",
                "eye_length_standard", "pronotum_width_standard",
                "abdomen_width_standard", "rear_leg_length_standard", 
                "rear_trochanter_length_standard", "Water_affinity", 
                "Flight_capability")
numeric_traits_standard <- c("body_length", "antenna_length_standard", "eye_protrusion_standard",
                    "eye_length_standard", "pronotum_width_standard",
                    "abdomen_width_standard", "rear_leg_length_standard", 
                    "rear_trochanter_length_standard")


# Examine correlation matrices ###############################################

cor_matrix <- cor(na.omit(traits[, trait_list]), method="pearson")
corrplot::corrplot(cor_matrix, method="number")

# Look to see if the correlation matrix with standardized traits has a reduction 
# in correlation
cor_matrix_standard <- cor(na.omit(traits[, trait_list_standard]), method="pearson")
corrplot::corrplot(cor_matrix_standard, method="number")
# Looks like eye protrusion standard is correlated with eye length standard,
# and rear leg length standard is correlated with antenna length standard. 

# Average the values of the (up to) six individuals ###########

traits_by_spp <- traits %>% group_by(Species) %>%
  summarize(num_indiv = n(),
            Genus = first(Genus),
            spp_abbrev = first(spp_abbrev),
            across(all_of(numeric_traits_standard), ~ mean(., na.rm=T)),
            Forest_affinity = first(Forest_affinity),
            Water_affinity = first(Water_affinity),
            Flight_capability = first(Flight_capability))
# Note: the ~ indicates a "lambda function" which allows me to tell the mean
# function to remove NA values when calculating the mean

# The only NA's in the dataset are for the water affinity or flight capability
# for Amerizus_unknown and Agonoleptus_thoracicus

# Look at the correlation matrix of the traits_by_spp table:
cor_matrix_by_spp <- cor(traits_by_spp[, numeric_traits_standard], method="pearson")
corrplot::corrplot(cor_matrix_by_spp, method="number")

# Create a data table with Notiophilus ommitted (see below):
traits_by_spp_wo_N <- traits_by_spp %>% filter(Species != "Notiophilus_aeneus")

# Run a principal components analysis #####################

pc <- prcomp(traits_by_spp[, c(numeric_traits_standard)], center = T, scale. = T)
# prcomp does a singular value decomposition of the (centered and
# scaled) data matrix

# visualize eigenvalues via a scree plot:
factoextra::get_eig(pc)
factoextra::fviz_pca_biplot(pc, axes=c(1,2))
factoextra::fviz_pca_biplot(pc, axes=c(2,3))
factoextra::fviz_pca_biplot(pc, axes=c(3,4))
pc$rotation

# Append the PC axes onto the traits_by_spp table:
traits_by_spp <- bind_cols(traits_by_spp, data.frame(pc$x))

# What is the correlation between PC1, PC2, PC3, PC4, Water_affinity, 
# and Flight_capability?
variables_for_distance_matrix <- c("PC1", "PC2", "PC3", "PC4", "Water_affinity", 
                                   "Flight_capability")
cor_matrix_pc <- cor(na.omit(traits_by_spp[, variables_for_distance_matrix]))
corrplot::corrplot(cor_matrix_pc, method="number")

# What is the variance of each PC axis?
var(traits_by_spp$PC1) # equal to the eigenvalue

# Graph scatter plots of the main findings of the PCA:

# PC1:
ggplot(data=traits_by_spp) + geom_point(aes(x=PC1, y=PC2))

# I'm worried that Notiophilus might be driving many of the observed patterns

# Rerun the PCA, but without Notiophilus ######################################

# The species Notiophilus aenius had just a few captures in the entire study.
# It has extremely big eyes, which I think might be having undue influence
# on the PC axes. Thus, I'll remove Notiophilus and see if the PC axes change.

pc_without_Notiophilus <- 
  prcomp(traits_by_spp_wo_N[, numeric_traits_standard], center = T, scale. = T)

plot(pc_without_Notiophilus)

factoextra::get_eig(pc_without_Notiophilus)
pc_without_Notiophilus$rotation

# Append the PC axes onto the traits_per_spp_without_Notiophilus table:
traits_by_spp_wo_N_1 <- 
  bind_cols(traits_by_spp_wo_N, data.frame(pc_without_Notiophilus$x))

# What is the correlation between PC1, PC2, PC3, Water_affinity, and Flight_capability?
cor_matrix2 <- cor(na.omit(traits_by_spp_wo_N_1[, variables_for_distance_matrix]), method="pearson")
corrplot::corrplot(cor_matrix2, method="number")
# PC2 is somewhat correlated with Flight capability, because smaller species
# with proportionally bigger eyes are likely to be flight capable.

# Make a finalized graph of the first two PC axes:
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(1,2), repel=T, 
                            geom.ind = "none", 
                            geom.var = c("arrow", "text"),
                            title="") + 
  coord_fixed() + theme_classic() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab("PC1 (38.7%)") + ylab("PC2 (26.4%)") +
  geom_text(data=traits_by_spp_wo_N_1,
            aes(x=PC1, y=PC2, label=spp_abbrev),
            size=3) + scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4),
                                         limits = c(-4.7, 4.7)) +
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3),
                     limits=c(-3.5,3.5))

# Now for PC1 vs PC3:
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(1,3), repel=T, 
                            geom.ind = "none", 
                            geom.var = c("arrow", "text"),
                            title="") + 
  coord_fixed() + theme_classic() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab("PC1 (38.7%)") + ylab("PC3 (13.6%)") +
  geom_text(data=traits_by_spp_wo_N_1,
            aes(x=PC1, y=PC3, label=spp_abbrev),
            size=3) + scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4),
                                         limits = c(-4.7, 4.7)) +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1,2,3), limits = c(-2.5,2.5))

# Now for PC1 vs PC4:
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(1,4), repel=T, 
                            geom.ind = "none", 
                            geom.var = c("arrow", "text"),
                            title="") + 
  coord_fixed() + theme_classic() + 
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab("PC1 (38.7%)") + ylab("PC4 (10.8%)") +
  geom_text(data=traits_by_spp_wo_N_1,
            aes(x=PC1, y=PC4, label=spp_abbrev),
            size=3) + scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4)) +
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3))

# Now I'd like to project Notiophilus onto the PC axes I found, so I can still
# include the species:
Notiophilus_traits_0 <- 
  traits_by_spp[traits_by_spp$Species == "Notiophilus_aeneus", 
                         c("Species", "num_indiv", "Genus", "spp_abbrev", 
                           numeric_traits_standard, "Forest_affinity",
                           "Water_affinity", "Flight_capability")]


numeric_traits_standard == names(pc_without_Notiophilus$center) # The order of traits is the same

# center and scale the Notiophilus data:
Notiophilus_traits <- Notiophilus_traits_0
Notiophilus_traits[numeric_traits_standard] <- 
  (Notiophilus_traits_0[numeric_traits_standard] - pc_without_Notiophilus$center) /
  pc_without_Notiophilus$scale

Notiophilus_traits$PC1 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC1,
              c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC2 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC2,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC3 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC3,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC4 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC4,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC5 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC5,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC6 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC6,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC7 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC7,
                c(t(Notiophilus_traits[numeric_traits_standard])))
Notiophilus_traits$PC8 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC8,
                c(t(Notiophilus_traits[numeric_traits_standard])))

# I want to keep the uncentered, unscaled trait data for the final data table,
# but use the projected PC values:
Notiophilus_traits_1 <- Notiophilus_traits
Notiophilus_traits_1[,trait_list_standard] <- 
  Notiophilus_traits_0[,trait_list_standard] # Moving the uncentered unscaled data
# back into the dataframe

# Join the Notiophilus data back onto the trait data table:
colnames(Notiophilus_traits_1) == colnames(traits_by_spp_wo_N_1) # column names are the same
combined_0 <- rbind(Notiophilus_traits_1, traits_by_spp_wo_N_1)

# rearrange the rows:
combined <- combined_0 %>% arrange(Species)

# Graph it for a reality check:
ggplot(data=combined, aes(x=PC1, y=PC2, label=spp_abbrev)) + 
  geom_text(alpha=0.5) + coord_fixed()
ggplot(data=combined, aes(x=PC1, y=PC3, label=spp_abbrev)) + 
  geom_text(alpha=0.5)+ coord_fixed()

# Export a table of the PC rotation (loading values):

rotation_df <- round(data.frame(pc_without_Notiophilus$rotation), 2)
rotation_df$trait <- rownames(rotation_df)

write.csv(rotation_df, "Aaron_PNR_formatted_data/PCA_traits_loadings.csv", row.names = F)

# Export a data table of traits ###############################################

#write.csv(combined, "Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv", row.names = F)

# Distance matrix incorporating water affinity and flight capability: #########

# Use Gower distance (see Swenson), incorporating PC1 through PC4, and 
# water affinity and Flight capability
combined$Water_affinity <- factor(combined$Water_affinity, order = T, 
                                     levels = c(0,0.5,1)) # Make this an ordered factor
combined$Flight_capability <- factor(combined$Flight_capability, order = T,
                                      levels = c(0,0.5,1))
gower_dist <- FD::gowdis(combined[, c(variables_for_distance_matrix)], ord="metric")

# try to visualize the distance matrix using a heatmap:
heatmap(as.matrix(gower_dist))
# I'm confused whether this is a symmetric matrix. It appears to be.

# What is the mean distance?
mean(gower_dist) # 0.27

# the min?
min(gower_dist) # 0.02
# the max?
max(gower_dist) # 0.80

# Export the distance matrix: #################################################
# Label the column names for each species:

gower_dist_dataframe <- as.data.frame(as.matrix(gower_dist))

combined$Species

colnames(gower_dist_dataframe) <- combined$Species

#write.csv(gower_dist_dataframe, file="Aaron_PNR_formatted_data/carabid_dist_in_trait_space.csv", row.names = F)

# Investigate trochanters of various species: ##################################

trochanter_ranking <- 
  combined %>% select(Species, rear_trochanter_length_standard) %>%
  arrange(rear_trochanter_length_standard)

# Species like Carabus goryi, Dicaelus teter, Sphaeroderus stenostomus,
# Notiophilus aeneus, Platynus angustatus, have short trochanters. Species like 
# Cyclotrachelus, Chlaenius, Harpalus spadiceous, Myas, have longer rear 
# trochanters.

# Investigate antenna length various species: ######################

antenna_ranking <- 
  combined %>% select(Species, antenna_length_standard) %>%
  arrange(antenna_length_standard)

hist(antenna_ranking$antenna_length_standard)

# Further investigate trait relationships: #####################################

ggplot(data=combined, aes(x=antenna_length_standard, y=rear_trochanter_length_standard,
       label=Species)) + 
  geom_point(alpha=0.5) + geom_text()

ggplot(data=combined %>% filter(spp_abbrev != "No.ae"), aes(x=antenna_length_standard, 
                                                            y=eye_length_standard)) + 
  geom_point(alpha=0.5)

ggplot(data=combined %>% filter(spp_abbrev != "No.ae"), aes(x=antenna_length_standard, 
                                                            y=eye_protrusion_standard)) + 
  geom_point(alpha=0.5)

