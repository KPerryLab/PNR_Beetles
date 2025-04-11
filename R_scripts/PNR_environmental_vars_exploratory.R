# Aaron Tayal
# Sept 24, 2024
# Exploratory data visualization of environmental variables recorded in 2022 at 
# Powdermill Nature Reserve in the tornado/salvage areas

library(ggplot2)
library(dplyr)
library(lubridate)

# Soil moisture 2022 ################################################################

soil_moisture_0 <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_Soil_moisture_PNR_ENV_2022.csv")
soil_moisture_0$SoilMoisture_day <- as.factor(ymd(soil_moisture_0$SoilMoisture_Date))

# Information on the transect groups and areas of each trap
trap_locations <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("integer", "integer", "factor", "factor", "factor", "numeric", "numeric", "numeric", "numeric"))
soil_moisture <- right_join(trap_locations, soil_moisture_0, by="Plot")

summary(soil_moisture)
plot(soil_moisture$percent_Moisture_Avg)
hist(soil_moisture$percent_Moisture_Avg)
hist(soil_moisture$percent_MOISTURE1)
hist(soil_moisture$percent_MOISTURE2)
hist(soil_moisture$percent_MOISTURE3)

ggplot(data=soil_moisture, aes(x=SoilMoisture_day, 
                               y=percent_Moisture_Avg)) +
  geom_point(alpha=0.5) + theme_classic()

# Idea: make a Parallel coordinates plot

# How much variation in soil moisture is related to forest disturbance
# type (undisturbed, wind-disturbed, or salvage-logged)? How much variation
# is related to transect? How much variation is related to date of measurement?

# Calculate the variance of percent soil moisture (averaging the three readings):
overall_var <- var(soil_moisture$percent_Moisture_Avg)
overall_mean <- mean(soil_moisture$percent_Moisture_Avg)
x_values <- -25:125
# Prob. density function of normal distribution:
y_values <- dnorm(x_values, overall_mean, sqrt(overall_var))
histogram <- hist(soil_moisture$percent_Moisture_Avg, plot=FALSE, breaks=20)
plot(histogram)
lines(x=x_values, y=1000*y_values)

# I think a normal distribution fits the data fairly well (subjectively).
ggplot(data=soil_moisture, aes(x=Treatment, y=percent_Moisture_Avg)) +
         geom_jitter(height=0, width=0.1, alpha=0.5)
# However, the above graph suggests that the within-group variance of 
# wind-disturbed plots (with multiple sample dates as separate replicates)
# is greater than that of salvage-logged plots. 

# Verify my observation by looking at the variances:
soil_moisture_by_treatment <- soil_moisture %>% group_by(Treatment) %>% 
  summarise(mean_moisture = mean(percent_Moisture_Avg),
            var_moisture = var(percent_Moisture_Avg),
            n_moisture = n())
soil_moisture_by_treatment

# Does this pattern still hold after taking the mean soil moisture for each plot 
# over the whole season?
soil_moisture_by_plot <- soil_moisture %>% group_by(Plot) %>% 
  summarise(Treatment = first(Treatment),
            mean_moisture = mean(percent_Moisture_Avg),
            var_moisture = var(percent_Moisture_Avg),
            n_moisture = n())
ggplot(data=soil_moisture_by_plot, aes(x=Treatment, y=mean_moisture)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)
# Yes, the pattern still holds that wind-disturbed sites seem more heterogeneous
# in their soil moisture compared to salvaged sites.

# Ground cover 2022 #################################################################

cover <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_percent_cover_PNR_ENV_2022.csv")
# It seems that the ground cover was evaluated four times in 2022: in June, 
# July, August, and September

cover$Treatment <- as.factor(cover$Treatment)
cover$Plot <- as.factor(cover$Plot)
cover$Date1 <- (ymd(cover$Date_yyyy.mm.dd))
cover$day_of_year <- yday(cover$Date_yyyy.mm.dd)

# What happens if you add up all the cover values?
cover$total_cover <- cover$VegAvg + cover$LitterAvg + cover$BareGAvg +
  cover$FWDAvg + cover$CWDAvg + cover$RockAvg

cover$total_cover1 <- cover$Veg1 + cover$Litter1 + cover$BareG1 +
  cover$FWD1 + cover$CWD1 + cover$Rock1

cover$total_cover2 <- cover$Veg2 + cover$Litter2 + cover$BareG2 +
  cover$FWD2 + cover$CWD2 + cover$Rock2

# I found a typo in the data: 
cover[cover$Date == "Aug 11 2022" & cover$Plot == 52, "Veg2"] # says 660
# It seems like this should have been 60%, not 660%
typo_row <- which(cover$Date == "Aug 11 2022" & cover$Plot == 52)
cover[typo_row, "total_cover2"] # says 700
cover[typo_row, "Veg2"] <- 60
cover[typo_row, "Veg1"]
cover[typo_row, "VegAvg"] <- 50

cover$total_cover <- cover$VegAvg + cover$LitterAvg + cover$BareGAvg +
  cover$FWDAvg + cover$CWDAvg + cover$RockAvg

cover$total_cover1 <- cover$Veg1 + cover$Litter1 + cover$BareG1 +
  cover$FWD1 + cover$CWD1 + cover$Rock1

cover$total_cover2 <- cover$Veg2 + cover$Litter2 + cover$BareG2 +
  cover$FWD2 + cover$CWD2 + cover$Rock2

# Plot some histograms:
hist(cover$VegAvg, breaks=10) # vegetation cover ranges from 0 to 100%, with
# most measurements around 10-60%
hist(cover$LitterAvg, breaks=10) # Litter ranges from 0 to 90%, with measurements
# fairly evenly distributed in that range
hist(cover$BareGAvg, breaks=10) # Bare ground seems relatively rare overall
hist(cover$FWDAvg, breaks=10) # Fine woody debris ranges from 0 to 25%, with most
# measurements being about 5 to 15% cover
hist(cover$CWDAvg, breaks=10) # CWD seems rare overall
hist(cover$RockAvg, breaks=10) # Rock cover ranges from 0 to ~30%, with the 
# majority of measurements in the 0-10% range

# Average across sample date: 
cover_by_plot <- cover %>% group_by(Plot) %>%
  summarize(Treatment = first(Treatment),
            VegAvg = mean(VegAvg),
            LitterAvg = mean(LitterAvg),
            BareGAvg = mean(BareGAvg),
            FWDAvg = mean(FWDAvg),
            CWDAvg = mean(CWDAvg),
            RockAvg = mean(RockAvg),
            VegHtAvg = mean(VegHtAvg))

# How does vegetation cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=VegAvg)) +
  geom_point(alpha=0.5) 
# I notice that Plot 49 (Salvaged) has a very high vegetation cover
# It is on the very southwest side of the whole area (see map)

# How does litter cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=LitterAvg)) +
  geom_point(alpha=0.5) 

# How does bare ground cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=BareGAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# A few of the Salvaged and Windthrow plots appear to have slightly
# higher bare ground %s than all of the Forest plots

# How does fine woody debris cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=FWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# How does coarse woody debris cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=CWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Although the percent cover of CWD is low in general, 
# it is suprising to me that a few of the forest plots have
# higher CWD %s than the windthrow or salvaged plots

# How does rock cover vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=RockAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# Vegetation height 2022 ##############################################

# Did vegetation height increase over the growing season?
ggplot(data=cover, aes(x=day_of_year, y=VegHtAvg)) +
  geom_point(alpha=0.5)
veg_height_change <- lm(VegHtAvg ~ Date_yyyy.mm.dd, data=cover)
summary(veg_height_change)
# Vegetation height appears to slightly decrease over time, with
# lowest height in September

# How does vegetation height vary by treatment?
ggplot(data=cover_by_plot, aes(x=Treatment, y=VegHtAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Again, plot 49 has abnormally high vegetation height
# Some of the windthrow plots seem to have taller vegetation than the 
# salvaged or forest plots

# Canopy openness 2022 #############################################################

densi <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_canopy_openness_PNR_ENV_2022.csv")
densi$Treatment <- as.factor(densi$Treatment)

# How does canopy openness vary by treatment?
ggplot(data=densi, aes(x=Treatment, y=Densi.Total)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Plot 49 has a very open canopy

ggplot(data=densi, aes(x=Treatment, y=Densi.Total)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) +
  ylim(0,20)
# Also, some of the other salvaged and windthrow plots have slightly
# higher canopy openness than all the other forest plots. But it seems like 
# differences in canopy openness have largely disappeared.

# Soil moisture 2015 ###########################################################
library(readxl)
soil_moisture_temp_2015_0 <- read_excel("PNR_Raw_Data/PNR_EnvironmentalData.xlsx", 
                                      sheet = 9)
# IMPORTANT NOTE: It appears that soil moisture was measured at plot 65 instead
# of plot 63. Looking at the map shows 65 is a short distance away from 63 and 
# also in the Windthrow treatment

soil_moisture_temp_2015 <- soil_moisture_temp_2015_0 %>% 
  filter(Quadrat >= 41 & Quadrat <= 65) %>% select(-"Block ID")
# Looks like soil moisture was taken at 6 time points

soil_moisture_temp_2015$SoilMoist1 <- as.numeric(soil_moisture_temp_2015$SoilMoist1)
# SoilMoist1 has some amount of missing data esp. towards end of summer

soil_moisture_temp_2015 %>% mutate(SoilMoistAvrg = rowMeans(SoilMoist1, SoilMoist2, 
                                                        SoilMoist3, na.rm=T))
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^ having trouble

                                              

# Ground cover 2015 ############################################################

# Vegetation height 2015 #######################################################

# Canopy openness 2015 #########################################################






