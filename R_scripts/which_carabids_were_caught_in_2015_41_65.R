# Aaron Tayal
# 3/18/2025
# Which carabid species were caught in 2015, within plots 41-65 at Powdermill?
# This script also creates a subsetted data table for 2015 with only species columns
# with non-zero total count

library(readxl)
library(dplyr)

dat_2015 <- read_excel("PNR_Raw_Data/PNR2015_InvertebrateCommunity.xlsx",
                       sheet = 5, na = ".")

dat_2015_subset <- dat_2015 %>% filter(Quadrat > 40)
dat_2015_subset$Treatment <- as.factor(dat_2015_subset$Treatment)
table(dat_2015_subset$Treatment)
sum(is.na(dat_2015_subset$`Agonum ferreum`)) # looks like 8 rows with missing data
sum(!is.na(dat_2015_subset$`Agonum ferreum`))

species_init <- colnames(dat_2015_subset[,8:64])

# Find out which species have zero abundance:
colSums(dat_2015_subset[,species_init], na.rm = T)
nonzero_species <- 
  species_init[colSums(dat_2015_subset[,species_init], na.rm = T) > 0]
colSums(dat_2015_subset[,nonzero_species], na.rm = T)

front_columns <- colnames(dat_2015_subset)[1:7]

# Remove any species columns if their total count is zero:
columns_to_keep <- c(front_columns, nonzero_species)
dat_2015_subset_drop_columns <- dat_2015_subset %>% 
  select(all_of(columns_to_keep))

# write a csv file with the subsetted data table:
write.csv(dat_2015_subset_drop_columns, 
          file="PNR2015_subset_carabid_counts.csv", row.names=FALSE)




