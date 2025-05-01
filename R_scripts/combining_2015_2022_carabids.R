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

all.equal(colnames(dat2022), colnames(dat2015))# good, the columns match

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

# Investigate final two collection intervals of 2022 ########################

# What carabid species were found in the final two collection intervals
# of 2022?
dat2022_final_2_int <- carab %>% filter(Year==2022 & Interval > 6)
colSums(dat2022_final_2_int[,carab_spp], na.rm=T)
# This makes me curious to investigate the seasonal trends of carabids more.

# What carabid species were found in the first six collection intervals
# of 2022?
dat2022_first_6_int <- carab %>% filter(Year==2022 & Interval < 7)
colSums(dat2022_first_6_int[,carab_spp], na.rm=T)

# It seems like different kinds of carabids are being caught in September,
# compared to the previous months. For instance: here is Platynus angustatus:

carab %>% filter(Year==2022) %>% ggplot(aes(x=Set_date, y=Platynus_angustatus)) +
  geom_point(alpha=0.2)
carab %>% filter(Year==2015) %>% ggplot(aes(x=Set_date, y=Platynus_angustatus)) +
  geom_point(alpha=0.2) # Because Platynus angustatus is mostly caught in 
# August and September, the fact that the 2022 sampling season extends into
# September means that 2015 and 2022 are not comparable.

# Write a data table: #######################################################

write.csv(carab, "PNR2015_2022_carabid_counts.csv", row.names = FALSE)







