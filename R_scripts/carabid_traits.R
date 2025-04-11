# Aaron Tayal
# 4/10/2025
# Purpose is to:
# 1. Find any correlation between the trait values
# 2. Divide traits by body length to standardize them
# 3. Run a principal components analysis to extract PC axes that explain 95% 
# of the variation
# 4. Create a distance matrix between species in trait space

library(tidyverse)
library(readxl)

traits <- read_excel("PNR_Raw_Data/PNR_SpeciesTraits_2022.xlsx", sheet=2)

traits$Rear_tarsi_length <- as.numeric(traits$Rear_tarsi_length) 
# one carabid has missing rear tarsi
traits$Forest_affinity <- as.factor(traits$Forest_affinity)
traits$Water_affinity <- as.numeric(traits$Water_affinity)
traits$Flight_capability <- as.numeric(traits$Flight_capability)

library(corrplot) # creates correlation visualizations

cor_matrix <- cor(na.omit(traits[, c(7:19, 22:23)]))
corrplot::corrplot(cor_matrix, method="ellipse")

# Now I'd like to calculate adjusted values for each specific trait:

traits$body_length <- traits$Elytra_length + traits$Pronotum_length + traits$Head_length

traits$antenna_length_standard <- traits$Antenna_length / traits$body_length

traits$eye_protrusion_standard <- 
  (traits$Outer_eye_distance - traits$Inner_eye_distance) / traits$body_length

traits$eye_length_standard <- traits$Eye_length / traits$body_length

traits$pronotum_width_standard <- traits$Pronotum_width / traits$body_length

traits$abdomen_width_standard <- traits$Abdomen_width / traits$body_length

traits$rear_leg_length_standard <- 
  ( traits$Rear_femur_length + traits$Rear_tibia_length + 
      traits$Rear_tarsi_length ) / traits$body_length

traits$rear_trochanter_length_standard <- traits$Rear_trochanter_length / 
  traits$body_length

trait_list <- c("body_length", "antenna_length_standard", "eye_protrusion_standard",
                "eye_length_standard", "pronotum_width_standard",
                "abdomen_width_standard", "rear_leg_length_standard", 
                "rear_trochanter_length_standard", "Water_affinity", 
                "Flight_capability")

# Look to see if the correlation matrix has a reduction in correlation
cor_matrix1 <- cor(na.omit(traits[, trait_list]))
corrplot::corrplot(cor_matrix1, method="ellipse")
# Looks like eye protrusion standard is correlated with eye length standard,
# and rear leg length standard is correlated with antenna length standard. Not
# sure if I need to remove those traits.

# Average the values of the six individuals ##################################

traits$Species <- as.factor(traits$Species)
traits_per_indiv <- traits %>% group_by(Species) %>%
  summarize(across(all_of(trait_list), mean))

# Run a principal component analysis: ##########################################

#rownames(traits_per_indiv) <- traits_per_indiv$Species # change row names to the 
# species names

pc <- prcomp(traits_per_indiv[, c(trait_list)], center=T, scale. = T)
# prcomp does a singular value decomposition of the (centered and possibly
# scaled) data matrix
attributes(pc)

library(factoextra)
library(ggplot2)

# visualize eigenvalues via a scree plot:
factoextra::get_eig(pc)
factoextra::fviz_pca_biplot(pc, axes=c(1,2))
factoextra::fviz_pca_biplot(pc, axes=c(2,3))









