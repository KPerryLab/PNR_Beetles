# 10/1/2024
# Aaron Tayal
# I would like to import the carabid count data from Powdermill 2022 
# properly

library(ggplot2)
library(dplyr)

# Import the carabid counts into the dataframe "carab". I would like to keep
# all column classes as character, because in the csv, blank cells in a species 
# column indicate true zeros, whereas "NA" indicates the trap contents were
# lost or destroyed by raccoons.
carab <- read.csv("Aaron_PNR_formatted_csvs/PNR2022_carabid_counts.csv",
                  na.strings = "NA", colClasses = "character")

sum(is.na(carab)) # Looks like we have 336 NA values. These are due to 
# 7 rows * 48 columns where the NA values occur.
carab[is.na(carab$Agonum_fidele), 1:3] # This lists the samples that are missing
# for various reasons (raccoons, mislabelling, misplacing vial)

carab_species <- colnames(carab)[which(colnames(carab)=="Agonum_fidele"):dim(carab)[2]] 
# This creates a vector containing the names of the carabid species.

# Now I want to replace blank cells with zeros for the species count columns.
# To do this, iterate over species and over rows, and only write a zero in 
# the column if there is currently a blank value "" and not if there is NA
for (species in carab_species) {
  for (row in 1:nrow(carab)){
    if (carab[row, species] == "" && !is.na(carab[row, species])){
      carab[row, species] <- "0"
    }
  }
}

# Now change the column classes to integers:
for (species in carab_species) {
  carab[,species] <- as.integer(carab[,species])
}
carab$Total_Carabidae_during_pinning <- as.integer(carab$Total_Carabidae_during_pinning)
carab$Total_Carabidae_from_sum_of_species_counts <- as.integer(carab$Total_Carabidae_from_sum_of_species_counts)

# Change the columns that indicate collection interval and trap location into 
# factors
carab$Collection_interval <- as.factor(carab$Collection_interval)
carab$PNR_Code <- as.factor(carab$PNR_Code)

# Graph the number of carabids captured, with trap location in the x-axis:
ggplot(data=carab, mapping=aes(x=PNR_Code, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()
# Now with collection interval on the x-axis:
ggplot(data=carab, mapping=aes(x=Collection_interval, y=Total_Carabidae_from_sum_of_species_counts)) +
  geom_violin(alpha=0.4) +
  theme_classic()

# Import information about the trap locations, such as transect and area:
trap_locations <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("integer", "integer", "factor", "factor", "factor", "numeric", "numeric"))


