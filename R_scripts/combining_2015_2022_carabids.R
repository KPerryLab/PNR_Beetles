# Aaron Tayal
# 3/19/2024
# Powdermill windthrow, salvaging experiment - carabid beetles study
# Combining the 2015 and 2022 datasets for carabid beetles

# ******************THIS SCRIPT NO LONGER WORKS BECAUSE THE SPECIES COLUMNS ARE 
# NOT MATCHING UP. *********************

# Outline:
# 1. I have a list of species which were either found in 2015 or 2022 in the 
# plots 41-65

# 2. I will import the 2015 and 2022 datasets

# 3. I will rename plot 65 in 2015 to plot 63 (which is located nearby - 
# see the Google map)

# 4. CONSIDERING taking away the final two collection intervals for 2022, because
# collections ended in August for the 2015 data

# 6. I will combine the two datasets based on matching column names to the 
# overall species list

library(dplyr)
library(ggplot2)

dat2015 <- read.csv("Aaron_PNR_formatted_data/PNR2015_subset_carabid_counts.csv",
                    na.strings = "NA", colClasses = "character")
# The 2015 data was subset to only include species columns where at least one
# individual was caught in plots 41-65.
# See the script: "which_carabids_were_caught_in_2015_41_65.R"

dat2022 <- read.csv("Aaron_PNR_formatted_data/PNR2022_carabid_counts.csv",
                    na.strings = "NA", colClasses = "character")

all.equal(colnames(dat2022), colnames(dat2015)) # UPDATE 8/5/2025: for some reason
# the 2015 dataset only has the columns for species found in 2015. So this 
# script is not working now. It's OK because I already have the combined
# dataset as a csv file. I think I must have manually added the species 
# columns to make them match up.

carab_spp <- colnames(dat2022)[7:63]

# Bind data tables together ##################################################

carab <- bind_rows(dat2015, dat2022)

# Add true zeros  #########################################

# Now I want to replace blank cells with zeros for the species count columns.
# To do this, iterate over species and over rows, and only write a zero in 
# the column if there is currently a blank value "" and not if there is NA
for (species in carab_spp) {
  for (row in 1:nrow(carab)){
    if (carab[row, species] == "" && !is.na(carab[row, species])){
      carab[row, species] <- "0"
    }
  }
}

# Change the column classes to the correct classes: ##########################

carab$Year <- as.factor(carab$Year) # year is a factor in the model

library(lubridate)
carab$Set_date <- mdy(carab$Set_date)
carab$Collection_date <- mdy(carab$Collection_date)
carab$Plot <- as.integer(carab$Plot)
carab$Treatment <- as.factor(carab$Treatment)
carab$Interval <- as.integer(carab$Interval)

for (species in carab_spp) {
  carab[,species] <- as.integer(carab[,species])
}

# Rename plot 65 to 63 #####################################################
carab$Plot == 65
carab$Plot[carab$Plot == 65] <- 63
carab$Plot

# Write a data table: #######################################################

#write.csv(carab, "PNR2015_2022_carabid_counts.csv", row.names = FALSE)







