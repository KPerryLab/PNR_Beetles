# Aaron Tayal
# 3/21/2025
# Purpose is to create a preliminary NMDS ordination of the carabid trait values
# so I can understand whether there might be a distinction between different
# locomotion methods

library(tidyverse)

traits <- read.csv("Aaron_PNR_formatted_data/PRELIM_carabid_traits.csv")

traits$Rear_tarsi_length <- as.numeric(traits$Rear_tarsi_length) 
# one carabid has missing rear tarsi

library(corrplot) # creates correlation visualizations

cor_matrix <- cor(traits[, c(5:16)]) # exclude rear tarsi length
corrplot::corrplot(cor_matrix, method="ellipse")

# Now I'd like to calculate adjusted values for each specific trait:

traits$body_length <- traits$Elytra_length + traits$Pronotum_length + traits$Head_length

traits$antenna_length_standard <- traits$Antenna_length / traits$body_length

traits$eye_protrusion_standard <- 
  (traits$Outer_eye_distance - traits$Inner_eye_distance) / traits$body_length

traits$eye_length_standard <- traits$Eye_length / traits$body_length

traits$pronotum_width_standard <- traits$Pronotum_width / traits$body_length

traits$abdomen_width_standard <- traits$Abdomen_width / traits$body_length

traits$rear_trochanter_length_standard <- traits$Rear_trochanter_length / 
  traits$body_length

trait_list <- c("body_length", "antenna_length_standard", "eye_protrusion_standard",
                "eye_length_standard", "pronotum_width_standard",
                "abdomen_width_standard", "rear_trochanter_length_standard")
# rear leg length standardized will have to be added to this eventually

traits[,trait_list]

# Look to see if the correlation matrix has a reduction in correlation
cor_matrix1 <- cor(traits[, trait_list]) # exclude rear tarsi length
corrplot::corrplot(cor_matrix1, method="ellipse")

# Average the values of the six individuals ##################################

traits$Species <- as.factor(traits$Species)
traits_per_indiv <- traits %>% group_by(Species) %>%
  summarize(across(all_of(trait_list), mean))

# Now run an NMDS to see if the beetles cluster into groups: ##################
library(vegan)

dist_matrix <- dist(traits_per_indiv[, trait_list], method="euclidean")
nmds_traits <- vegan::metaMDS(dist_matrix, trymax = 500, k=2)
nmds_traits
stressplot(nmds_traits)
plot(nmds_traits)
text(nmds_traits, labels = traits_per_indiv$Species)

# Run a principal component analysis: ##########################################

#rownames(traits_per_indiv) <- traits_per_indiv$Species # change row names to the 
# species names

pc <- prcomp(traits_per_indiv[, c(trait_list)], center=T, scale. = T)
attributes(pc)

library(factoextra)
library(ggplot2)

# visualize eigenvalues via a scree plot:
factoextra::get_eig(pc)
factoextra::fviz_pca_biplot(pc, axes=c(1,2))
factoextra::fviz_pca_biplot(pc, axes=c(2,3))









