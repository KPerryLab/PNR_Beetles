# Aaron Tayal
# UPDATED VERSION, May 5, 2025

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

traits <- read_excel("PNR_Raw_Data/PNR_SpeciesTraits_2022.xlsx", sheet=2, na='NA')

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

# Calculate some traits that are sums or differences of measurements ###########

# Body length:
traits$body_length <- traits$Elytra_length + traits$Pronotum_length + traits$Head_length

# Eye protrusion
traits$eye_protrusion <- (traits$Outer_eye_distance - traits$Inner_eye_distance)

# Rear leg length
traits$rear_leg_length <- 
  ( traits$Rear_femur_length + traits$Rear_tibia_length + 
      traits$Rear_tarsi_length )

# Calculate some traits that are ratios of two measurements ####################

traits$eye_protrusion_ratio <- traits$eye_protrusion / traits$Eye_length

traits$antenna_rear_leg_ratio <- traits$Antenna_length / traits$rear_leg_length

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

# Look at a set of traits that includes eye protrusion ratio (a ratio to eye length)
# as well as antenna:rear_leg ratio:
trait_list_modified <- c("body_length", "eye_length_standard",
                         "eye_protrusion_ratio", "pronotum_width_standard",
                         "abdomen_width_standard", "rear_leg_length_standard", 
                         "antenna_rear_leg_ratio",
                         "rear_trochanter_length_standard", "Forest_affinity",
                         "Water_affinity", 
                         "Flight_capability")

numeric_traits_modified <- c("body_length", "eye_length_standard",
                                   "eye_protrusion_ratio", "pronotum_width_standard",
                                   "abdomen_width_standard", "rear_leg_length_standard", 
                                   "antenna_rear_leg_ratio",
                                   "rear_trochanter_length_standard")

cor_matrix_modified <- cor(na.omit(traits[, numeric_traits_modified]), method="pearson")
corrplot::corrplot(cor_matrix_modified, method="number")
# This modified set of traits appears better because the traits have lower
# Pearson correlation coefficients than before

# Make some scatter plots of traits against each other: ######################
# Making scatter plots helps to determine if traits are proportional to one another
# and helps to detect outliers (possible data entry mistakes)

ggplot(data=traits) + geom_point(aes(x=body_length, y=Antenna_length), alpha=0.5) +
  xlim(0,35) + ylim(0,20)

ggplot(data=traits) + geom_point(aes(x=body_length, y=eye_protrusion), alpha=0.5) +
  xlim(0,30) + ylim(0,2)
ggplot(data=traits) + geom_point(aes(x=body_length, y=eye_protrusion), alpha=0.5) +
  xlim(0,8) + ylim(0,2)
traits[which(traits$body_length < 6),"Species"]
# It appears that Notiophilus has a slightly bigger eye protrusion than what would 
# be expected given its body size

ggplot(data=traits) + geom_point(aes(x=body_length, y=Eye_length), alpha=0.5) +
  xlim(0,35) + ylim(0,2) 

# Does eye protrusion ratio need to be standardized to body length?
ggplot(data=traits) + geom_point(aes(x=body_length, y=eye_protrusion_ratio), alpha=0.5)
# No, the eye protrusion ratio doesn't seem to be tightly correlated with body length

ggplot(data=traits) + geom_point(aes(x=body_length, y=Pronotum_width), alpha=0.5)
traits[which(traits$body_length > 20 & traits$Pronotum_width < 5),"Species"] # scaphinotus andrewsii
# has a very narrow pronotum
traits[which(traits$body_length > 17 & traits$Pronotum_width < 4),"Species"] # as well
# as Galerita bicolor

ggplot(data=traits) + geom_point(aes(x=body_length, y=Pronotum_length), alpha=0.5)

ggplot(data=traits) + geom_point(aes(x=body_length, y=Abdomen_width), alpha=0.5)
traits[which(traits$Abdomen_width > 10),"Species"] # Scaphinotus viduus

ggplot(data=traits) + geom_point(aes(x=body_length, y=rear_leg_length), alpha=0.5)
traits[which(traits$body_length < 15 & traits$rear_leg_length > 13),"Species"] # Platynus angustatus
# has very long legs

ggplot(data=traits) + geom_point(aes(x=body_length, y=Rear_trochanter_length), alpha=0.5)

# The following graph shows how eye protrusion is correlated with eye length:
ggplot(data=traits) + geom_point(aes(x=eye_length_standard, y=eye_protrusion_standard), alpha=0.5)

# How does eye protrusion ratio relate to eye length?
ggplot(data=traits) + geom_point(aes(x=Eye_length, y=eye_protrusion_ratio), alpha=0.5) +
  xlim(0,2) + ylim(0,2)

# How does antenna:rear leg ratio relate to rear leg?
ggplot(data=traits) + geom_point(aes(x=rear_leg_length, y=antenna_rear_leg_ratio), alpha=0.5)
# Observation: the antennae measurements are likely more variable due to the 
# difficulties in measuring antenna length

# MODIFIED VERSION Average the values of the (up to) six individuals ###########

traits$Species <- as.factor(traits$Species)
traits_by_spp_modified <- traits %>% group_by(Species) %>%
  summarize(Genus = first(Genus),
            across(all_of(numeric_traits_modified), ~ mean(., na.rm=T)),
            Forest_affinity = first(Forest_affinity),
            Water_affinity = first(Water_affinity),
            Flight_capability = first(Flight_capability))
# Note: the ~ indicates a "lambda function" which allows me to tell the mean
# function to remove NA values when calculating the mean

# The only NA's in the dataset are for the water affinity or flight capability
# for Amerizus_unknown and Agonoleptus_thoracicus

# Look at the correlation matrix of the traits_by_spp table:
cor_matrix_by_spp_modified <- cor(na.omit(traits_by_spp_modified[, numeric_traits_modified]), method="pearson")
corrplot::corrplot(cor_matrix_by_spp_modified, method="number")
# All correlation coefficients are below 0.6, except for rear leg length against 
# pronotum width


# Create a data table with Notiophilus ommitted (see below):
traits_by_spp_modified_without_N <- traits_by_spp_modified %>%
  filter(Species != "Notiophilus_aeneus")

# MODIFIED TRAIT LIST Run a principal components analysis #####################

pc_modified <- prcomp(traits_by_spp_modified[, c(numeric_traits_modified)], center = T, scale. = T)
# prcomp does a singular value decomposition of the (centered and
# scaled) data matrix

# visualize eigenvalues via a scree plot:
factoextra::get_eig(pc_modified)
factoextra::fviz_pca_biplot(pc_modified, axes=c(1,2))
factoextra::fviz_pca_biplot(pc_modified, axes=c(2,3))
factoextra::fviz_pca_biplot(pc_modified, axes=c(3,4))
pc_modified$rotation

# Append the PC axes onto the traits_per_indiv_modified table
traits_by_spp_modified <- 
  bind_cols(traits_by_spp_modified, data.frame(pc_modified$x))

# What is the correlation between PC1, PC2, PC3, PC4, Water_affinity, 
# and Flight_capability?
variables_for_distance_matrix <- c("PC1", "PC2", 
                                   "PC3", "PC4", "Water_affinity", 
                                   "Flight_capability")
cor_matrix_pc_modified <- cor(na.omit(traits_by_spp_modified[, variables_for_distance_matrix]))
corrplot::corrplot(cor_matrix_pc_modified, method="number")
# PC2 is somewhat correlated with Flight capability, because smaller species
# with proportionally bigger eyes are likely to be flight capable.

# What is the variance of each PC axis?
var(traits_by_spp_modified$PC1) # equal to the eigenvalue

# Do species within a Genus tend to group together in the PCA?
factoextra::fviz_pca_biplot(pc_modified, axes=c(1,2), 
                            habillage = traits_by_spp_modified$Genus)
factoextra::fviz_pca_biplot(pc_modified, axes=c(2,3), 
                            habillage = traits_by_spp_modified$Genus)

# What is the variance of water affinity?
var(traits_by_spp_modified$Water_affinity, na.rm = T)
table(traits_by_spp_modified$Water_affinity)

# What is the variance of Flight capability?
var(traits_by_spp_modified$Flight_capability, na.rm = T)
table(traits_by_spp_modified$Flight_capability)

# Graph scatter plots of the main findings of the PCA:

# PC1:
ggplot(data=traits_by_spp_modified) + geom_point(aes(x=rear_leg_length_standard,
                                                     y=pronotum_width_standard))

# PC2:
ggplot(data=traits_by_spp_modified) + geom_point(aes(x=body_length,
                                                     y=eye_length_standard))

# PC3:
ggplot(data=traits_by_spp_modified) + geom_point(aes(x=antenna_rear_leg_ratio,
                                                     y=eye_protrusion_ratio))

# PC4
ggplot(data=traits_by_spp_modified) + geom_point(aes(x=abdomen_width_standard,
                                                     y=eye_length_standard))

# I'm worried that Notiophilus might be driving many of the observed patterns

# Rerun the modified PCA, but without Notiophilus #############################

# Remove Notiophilus and run another principal component analysis: #############

# The species Notiophilus aenius had just 1 capture in the entire dataset.
# It has extremely big eyes, which I think might be having undue influence
# on the PC axes. Thus, I'll remove Notiophilus and see if the PC axes change.

pc_without_Notiophilus <- 
  prcomp(traits_by_spp_modified_without_N[, numeric_traits_modified], 
         center = T, scale. = T)

plot(pc_without_Notiophilus)

factoextra::get_eig(pc_without_Notiophilus)
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(1,2), repel=T, 
                            geom.ind = c("point", "text"), 
                            geom.var = c("arrow", "text"),
                            title="PCA of ground beetle traits") + 
  coord_fixed() + theme_classic() + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(2,3), repel=T, 
                            geom.ind = c("point", "text"), 
                            geom.var = c("arrow", "text"),
                            title="PCA of ground beetle traits") + 
  coord_fixed() + theme_classic() + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3)) + 
  scale_y_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(3,4))
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(4,5))
factoextra::fviz_pca_biplot(pc_without_Notiophilus, axes=c(5,6))
pc_without_Notiophilus$rotation

# Append the PC axes onto the traits_per_spp_without_Notiophilus table:
traits_by_spp_modified_without_N_1 <- 
  bind_cols(traits_by_spp_modified_without_N, data.frame(pc_without_Notiophilus$x))

# What is the correlation between PC1, PC2, PC3, Water_affinity, and Flight_capability?
cor_matrix2 <- cor(na.omit(traits_by_spp_modified_without_N_1[, c(12:16)]), method="pearson")
corrplot::corrplot(cor_matrix2, method="number")
# PC2 is somewhat correlated with Flight capability, because smaller species
# with proportionally bigger eyes are likely to be flight capable.

# Now I'd like to project Notiophilus onto the PC axes I found, so I can still
# include the species:
Notiophilus_traits_0 <- 
  traits_by_spp_modified[traits_by_spp_modified$Species == "Notiophilus_aeneus", 
                         c("Species", "Genus", trait_list_modified)]

# center and scale:
Notiophilus_traits <- Notiophilus_traits_0
Notiophilus_traits[numeric_traits_modified] <- 
  (Notiophilus_traits_0[numeric_traits_modified] - pc_without_Notiophilus$center) /
  pc_without_Notiophilus$scale

Notiophilus_traits$PC1 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC1,
              c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC2 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC2,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC3 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC3,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC4 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC4,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC5 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC5,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC6 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC6,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC7 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC7,
                c(t(Notiophilus_traits[numeric_traits_modified])))
Notiophilus_traits$PC8 <- 
  geometry::dot(as.data.frame(pc_without_Notiophilus$rotation)$PC8,
                c(t(Notiophilus_traits[numeric_traits_modified])))

# I want to keep the uncentered, unscaled trait data for the final data table,
# but use the projected PC values:
Notiophilus_traits_1 <- Notiophilus_traits
Notiophilus_traits_1[,trait_list_modified] <- 
  Notiophilus_traits_0[,trait_list_modified] # Moving the uncentered unscaled data
# back into the dataframe

# Join the Notiophilus data back onto the trait data table:
colnames(Notiophilus_traits_1) == colnames(traits_by_spp_modified_without_N_1) # column names are the same
combined_0 <- rbind(Notiophilus_traits_1, traits_by_spp_modified_without_N_1)

# rearrange the rows:
combined <- combined_0 %>% arrange(Species)

# Graph it for a reality check:
ggplot(data=combined, aes(x=PC1, y=PC2, label=Species)) + geom_point() + 
  geom_text(alpha=0.5) + coord_fixed()
ggplot(data=combined, aes(x=PC2, y=PC3, label=Species)) + geom_point() + 
  geom_text(alpha=0.5)+ coord_fixed()
ggplot(data=combined, aes(x=PC3, y=PC4, label=Species)) + geom_point() + 
  geom_text(alpha=0.5)+ coord_fixed() # Looks good

# Add back in the standardized antenna length (because I want to look at the 
# CWM antenna length):
traits_by_spp_antenna <- traits %>% group_by(Species) %>% 
  summarise(Species=first(Species),
            antenna_length_standard = mean(antenna_length_standard))
traits_by_spp_antenna$Species == combined$Species

combined$antenna_length_standard <- traits_by_spp_antenna$antenna_length_standard

# Add back in the standardized eye protrusion (because I want to look at the 
# CWM eye protrusion):
traits_by_spp_eye_protrusion <- traits %>% group_by(Species) %>% 
  summarise(Species=first(Species),
            eye_protrusion_standard = mean(eye_protrusion_standard))
traits_by_spp_eye_protrusion$Species == combined$Species

combined$eye_protrusion_standard <- traits_by_spp_eye_protrusion$eye_protrusion_standard

# Export a table of the PC rotation (loading values):

#write.csv(pc_without_Notiophilus$rotation, "Aaron_PNR_formatted_data/PCA_traits_loadings.csv")

# Export a data table of traits ###############################################

#write.csv(combined, "Aaron_PNR_formatted_data/PNR_carabid_traits_by_spp.csv", row.names = F)

# Distance matrix incorporating water affinity and flight capability: #########

# Use Gower distance (see Swenson), incorporating PC1 through PC4, and 
# water affinity and Flight capability
combined$Water_affinity <- factor(combined$Water_affinity, order = T, 
                                     levels = c(0,0.5,1))
combined$Flight_capability <- factor(combined$Flight_capability, order = T,
                                      levels = c(0,0.5,1))
gower_dist <- FD::gowdis(combined[, c(variables_for_distance_matrix)], ord="metric")

# try to visualize the distance matrix using a heatmap:
heatmap(as.matrix(gower_dist))
# I'm confused whether this is a symmetric matrix. It appears to be.

# What is the mean distance?
mean(gower_dist) # 0.29

# the max?
max(gower_dist) # 0.69

# Export the distance matrix: #################################################
# Label the column names for each species:

gower_dist_dataframe <- as.data.frame(as.matrix(gower_dist))

colnames(gower_dist_dataframe) <- traits_by_spp_modified$Species

#write.csv(gower_dist_dataframe, file="Aaron_PNR_formatted_data/carabid_dist_in_trait_space.csv", row.names = F)

# Investigate trochanters of various species: ##################################

trochanter_ranking <- 
  traits_by_spp_modified %>% select(Species, rear_trochanter_length_standard) %>%
  arrange(rear_trochanter_length_standard)




