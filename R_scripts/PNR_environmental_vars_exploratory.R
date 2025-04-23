# Aaron Tayal
# Sept 24, 2024 updated Apr 22, 2025
# Exploratory data visualization of environmental variables recorded at 
# Powdermill Nature Reserve in the tornado/salvage areas

library(ggplot2)
library(dplyr)
library(lubridate)

# Soil moisture 2022 ################################################################

soil_moisture_2022_0 <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_Soil_moisture_PNR_ENV_2022.csv")
soil_moisture_2022_0$SoilMoisture_day <- as.factor(ymd(soil_moisture_2022_0$SoilMoisture_Date))

# Information on the transect groups and areas of each trap
trap_locations <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("integer", "integer", "factor", "factor", "factor", "numeric", "numeric", "numeric", "numeric"))
soil_moisture_2022 <- right_join(trap_locations, soil_moisture_2022_0, by="Plot")

summary(soil_moisture_2022)
plot(soil_moisture_2022$percent_Moisture_Avg)
hist(soil_moisture_2022$percent_Moisture_Avg)
hist(soil_moisture_2022$percent_MOISTURE1)
hist(soil_moisture_2022$percent_MOISTURE2)
hist(soil_moisture_2022$percent_MOISTURE3)

ggplot(data=soil_moisture_2022, aes(x=SoilMoisture_day, 
                               y=percent_Moisture_Avg, group=Plot)) +
  geom_point(alpha=0.5) + geom_line(alpha=0.5) + theme_classic()

# How much variation in soil moisture is related to forest disturbance
# type (undisturbed, wind-disturbed, or salvage-logged)? How much variation
# is related to transect? How much variation is related to date of measurement?

# Calculate the variance of percent soil moisture (averaging the three readings):
overall_var <- var(soil_moisture_2022$percent_Moisture_Avg)
overall_mean <- mean(soil_moisture_2022$percent_Moisture_Avg)
x_values <- -25:125
# Prob. density function of normal distribution:
y_values <- dnorm(x_values, overall_mean, sqrt(overall_var))
histogram <- hist(soil_moisture_2022$percent_Moisture_Avg, plot=FALSE, breaks=20)
plot(histogram)
lines(x=x_values, y=1000*y_values)

# I think a normal distribution fits the data fairly well (subjectively).
ggplot(data=soil_moisture_2022, aes(x=Treatment, y=percent_Moisture_Avg)) +
         geom_jitter(height=0, width=0.1, alpha=0.5)
# However, the above graph suggests that the within-group variance of 
# wind-disturbed plots (with multiple sample dates as separate replicates)
# is greater than that of salvage-logged plots. 

# Verify my observation by looking at the variances:
soil_moisture_2022_by_treatment <- soil_moisture_2022 %>% group_by(Treatment) %>% 
  summarise(mean_moisture = mean(percent_Moisture_Avg),
            var_moisture = var(percent_Moisture_Avg),
            n_moisture = n())
soil_moisture_2022_by_treatment

# Does this pattern still hold after taking the mean soil moisture for each plot 
# over the whole season?
soil_moisture_2022_by_plot <- soil_moisture_2022 %>% group_by(Plot) %>% 
  summarise(Treatment = first(Treatment),
            mean_moisture = mean(percent_Moisture_Avg),
            var_moisture = var(percent_Moisture_Avg),
            n_moisture = n())
ggplot(data=soil_moisture_2022_by_plot, aes(x=Treatment, y=mean_moisture)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)
# Yes, the pattern still holds that wind-disturbed sites seem more heterogeneous
# in their soil moisture compared to salvaged sites.

# Ground cover 2022 #################################################################

cover_2022 <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_percent_cover_PNR_ENV_2022.csv")
# It seems that the ground cover was evaluated four times in 2022: in June, 
# July, August, and September

cover_2022$Treatment <- as.factor(cover_2022$Treatment)
cover_2022$Plot <- as.factor(cover_2022$Plot)
cover_2022$Date1 <- (ymd(cover_2022$Date_yyyy.mm.dd))
cover_2022$day_of_year <- yday(cover_2022$Date_yyyy.mm.dd)

# What happens if you add up all the cover values?
cover_2022$total_cover <- cover_2022$VegAvg + cover_2022$LitterAvg + 
  cover_2022$BareGAvg + cover_2022$FWDAvg + cover_2022$CWDAvg + cover_2022$RockAvg

cover_2022$total_cover1 <- cover_2022$Veg1 + cover_2022$Litter1 + cover_2022$BareG1 +
  cover_2022$FWD1 + cover_2022$CWD1 + cover_2022$Rock1

cover_2022$total_cover2 <- cover_2022$Veg2 + cover_2022$Litter2 + cover_2022$BareG2 +
  cover_2022$FWD2 + cover_2022$CWD2 + cover_2022$Rock2

# I found a typo in the data: 
cover_2022[cover_2022$Date == "Aug 11 2022" & cover_2022$Plot == 52, "Veg2"] # says 660
# It seems like this should have been 60%, not 660%
typo_row <- which(cover_2022$Date == "Aug 11 2022" & cover_2022$Plot == 52)
cover_2022[typo_row, "total_cover2"] # says 700
cover_2022[typo_row, "Veg2"] <- 60
cover_2022[typo_row, "Veg1"]
cover_2022[typo_row, "VegAvg"] <- 50

cover_2022$total_cover <- cover_2022$VegAvg + cover_2022$LitterAvg + cover_2022$BareGAvg +
  cover_2022$FWDAvg + cover_2022$CWDAvg + cover_2022$RockAvg

cover_2022$total_cover1 <- cover_2022$Veg1 + cover_2022$Litter1 + cover_2022$BareG1 +
  cover_2022$FWD1 + cover_2022$CWD1 + cover_2022$Rock1

cover_2022$total_cover2 <- cover_2022$Veg2 + cover_2022$Litter2 + cover_2022$BareG2 +
  cover_2022$FWD2 + cover_2022$CWD2 + cover_2022$Rock2

# Plot some histograms:
hist(cover_2022$VegAvg, breaks=10) # vegetation cover ranges from 0 to 100%, with
# most measurements around 10-60%
hist(cover_2022$LitterAvg, breaks=10) # Litter ranges from 0 to 90%, with measurements
# fairly evenly distributed in that range
hist(cover_2022$BareGAvg, breaks=10) # Bare ground seems relatively rare overall
hist(cover_2022$FWDAvg, breaks=10) # Fine woody debris ranges from 0 to 25%, with most
# measurements being about 5 to 15% cover
hist(cover_2022$CWDAvg, breaks=10) # CWD seems rare overall
hist(cover_2022$RockAvg, breaks=10) # Rock cover ranges from 0 to ~30%, with the 
# majority of measurements in the 0-10% range

# Average across sample date: 
cover_2022_by_plot <- cover_2022 %>% group_by(Plot) %>%
  summarize(Treatment = first(Treatment),
            VegAvg = mean(VegAvg),
            LitterAvg = mean(LitterAvg),
            BareGAvg = mean(BareGAvg),
            FWDAvg = mean(FWDAvg),
            CWDAvg = mean(CWDAvg),
            RockAvg = mean(RockAvg),
            VegHtAvg = mean(VegHtAvg))

# How does vegetation cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=VegAvg)) +
  geom_point(alpha=0.5) 
# I notice that Plot 49 (Salvaged) has a very high vegetation cover
# It is on the very southwest side of the whole area (see map)

# How does litter cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=LitterAvg)) +
  geom_point(alpha=0.5) 

# How does bare ground cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=BareGAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# A few of the Salvaged and Windthrow plots appear to have slightly
# higher bare ground %s than all of the Forest plots

# How does fine woody debris cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=FWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# How does coarse woody debris cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=CWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Although the percent cover of CWD is low in general, 
# it is suprising to me that a few of the forest plots have
# higher CWD %s than the windthrow or salvaged plots. Note that this is just
# the % cover, not the actual volume of CWD

# How does rock cover vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=RockAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# Vegetation height 2022 ##############################################

# Did vegetation height increase over the growing season?
ggplot(data=cover_2022, aes(x=day_of_year, y=VegHtAvg)) +
  geom_point(alpha=0.5)
veg_height_change_2022 <- lm(VegHtAvg ~ Date_yyyy.mm.dd, data=cover_2022)
summary(veg_height_change_2022)
# Vegetation height appears to slightly decrease over time, with
# lowest height in September

# How does vegetation height vary by treatment?
ggplot(data=cover_2022_by_plot, aes(x=Treatment, y=VegHtAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Again, plot 49 has abnormally high vegetation height
# Some of the windthrow plots seem to have taller vegetation than the 
# forest plots, but I haven't tested this with stats

# Canopy openness 2022 #############################################################

densi_2022 <- read.csv("Aaron_PNR_formatted_data/Aaron_formatted_canopy_openness_PNR_ENV_2022.csv")
densi_2022$Treatment <- as.factor(densi_2022$Treatment)

# How does canopy openness vary by treatment?
ggplot(data=densi_2022, aes(x=Treatment, y=Densi.Total)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Plot 49 has a very open canopy

ggplot(data=densi_2022, aes(x=Treatment, y=Densi.Total)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) +
  ylim(0,20)
# Also, some of the other salvaged and windthrow plots have slightly
# higher canopy openness than all the other forest plots. But it seems like any
# differences in canopy openness have mostly disappeared by 2022, except for 
# plot 49

# Soil moisture and temperature 2015 ##########################################

soil_moisture_temp_2015_0 <- 
  read.csv("Aaron_PNR_formatted_data/PNR_EnvironmentalData_SoilMoistureTemp_2015.csv")
# IMPORTANT NOTE: It appears that soil moisture was measured at plot 65 instead
# of plot 63. Looking at the map shows 65 is a short distance away from 63 and 
# also in the Windthrow treatment

soil_moisture_temp_2015_0$Day <- 
  as.factor(mdy(soil_moisture_temp_2015_0$Date)) 

soil_moisture_temp_2015 <- soil_moisture_temp_2015_0 %>% 
  filter(Plot >= 41 & Plot <= 65)
# Looks like soil moisture was taken at 6 time points

soil_moisture_temp_2015$SoilMoist1 <- as.numeric(soil_moisture_temp_2015$SoilMoist1)
# SoilMoist1 has some amount of missing data esp. towards end of summer

# I want to take the average of the three readings,
# but sometimes there are only 2 readings and the other one is listed as NA.
# But I still want to take the average of those two readings.
soil_moisture_temp_2015$SoilMoistAvrg <- rowMeans(data.frame(
  soil_moisture_temp_2015$SoilMoist1, soil_moisture_temp_2015$SoilMoist2,
  soil_moisture_temp_2015$SoilMoist3), na.rm=T)

# Same routine for soil temp:
soil_moisture_temp_2015$SoilTempAvrg <- rowMeans(data.frame(
  soil_moisture_temp_2015$SoilTemp1, soil_moisture_temp_2015$SoilTemp2,
  soil_moisture_temp_2015$SoilTemp3), na.rm=T)

# How did soil moisture change over the season in 2015?
ggplot(data=soil_moisture_temp_2015, aes(x=Day, 
                                    y=SoilMoistAvrg)) +
  geom_point(alpha=0.5) + theme_classic()
# Looks like the soil moisture strongly decreased over the season in 2015.
# Note: the spacing between measurement dates is not to scale

# How did soil temp change over the season in 2015?
ggplot(data=soil_moisture_temp_2015, aes(x=Day, 
                                         y=SoilTempAvrg)) +
  geom_point(alpha=0.5) + theme_classic()
# The soil temp was higher in August in 2015

# Summarize by plot:
soil_moisture_temp_2015_by_plot <- soil_moisture_temp_2015 %>% group_by(Plot) %>% 
  summarise(Treatment = first(Treatment),
            mean_moisture = mean(SoilMoistAvrg),
            var_moisture = var(SoilMoistAvrg),
            mean_temp = mean(SoilTempAvrg),
            var_temp = var(SoilTempAvrg),
            n_dates_measured = n())

# Graph soil moisture in 2015 by plot:
ggplot(data=soil_moisture_temp_2015_by_plot, aes(x=Treatment, y=mean_moisture)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)

# Graph soil temp in 2015 by plot:
ggplot(data=soil_moisture_temp_2015_by_plot, aes(x=Treatment, y=mean_temp)) +
  geom_jitter(height=0, width=0.1, alpha=0.5)
# Some of the windthrow and salvaged plots had higher mean temp for 2015.

# Ground cover 2015 ############################################################

cover_2015_0 <- read.csv("Aaron_PNR_formatted_data/PNR_EnvironmentalData_DensiVeg_2015.csv")
# It seems that the ground cover was evaluated twice in 2015: on 9 June, 
# and 7 July

# Note: the column "Densi_Date" probably refers to when the canopy openness 
# readings were taken, which were on 9-10 June and 5 August according to Kayla.
# Meanwhile the column "Date" refers to when the ground cover estimates were taken.


# Note: this data table also has canopy openness values (DENSI)

cover_2015 <- cover_2015_0 %>% filter(Plot >= 41 & Plot <= 65)

cover_2015$Treatment <- as.factor(cover_2015$Treatment)
cover_2015$Densi_Date <- as.factor(cover_2015$Densi_Date)
cover_2015$Date <- as.factor(cover_2015$Date)

# What happens if you add up all the cover values?

cover_2015$total_cover1 <- cover_2015$Veg1 + cover_2015$Litter1 + cover_2015$BareG1 +
  cover_2015$FWD1 + cover_2015$CWD1 + cover_2015$Rock1

cover_2015$total_cover2 <- cover_2015$Veg2 + cover_2015$Litter2 + cover_2015$BareG2 +
  cover_2015$FWD2 + cover_2015$CWD2 + cover_2015$Rock2
# Most add up to 100%, but one or two rows just add up to 40%, 95%, or 105%

# This time I need to take the average of the two readings, for each cover type:
cover_2015$VegAvg <- ( cover_2015$Veg1 + cover_2015$Veg2 ) / 2
cover_2015$LitterAvg <- ( cover_2015$Litter1 + cover_2015$Litter2 ) / 2
cover_2015$BareGAvg <- ( cover_2015$BareG1 + cover_2015$BareG2 ) / 2
cover_2015$FWDAvg <- ( cover_2015$FWD1 + cover_2015$FWD2 ) / 2
cover_2015$CWDAvg <- ( cover_2015$CWD1 + cover_2015$CWD2 ) / 2
cover_2015$RockAvg <- ( cover_2015$Rock1 + cover_2015$Rock2 ) / 2
# I'll do veg height while I'm at it:
cover_2015$VegHtAvg <- ( cover_2015$VegHt1 + cover_2015$VegHt2 ) / 2

# Plot some histograms:
hist(cover_2015$VegAvg, breaks=10) # vegetation cover ranges from 0 to 100%, with
# measurements fairly evenly disturbuted, but most plots below 50%
hist(cover_2015$LitterAvg, breaks=10) # Litter ranges from 0 to 90%, with measurements
# fairly evenly distributed in that range
hist(cover_2015$BareGAvg, breaks=10) # Bare ground ranges from 0 to 25%, with
# most in the 0-5% range but a few higher.
hist(cover_2015$FWDAvg, breaks=10) # Fine woody debris ranges from 0 to 30%, with most
# measurements being about 0-10% cover
hist(cover_2015$CWDAvg, breaks=10) # CWD seems rare overall, but does range up to
# 35%
hist(cover_2022$RockAvg, breaks=10) # Rock cover ranges from 0 to ~30%, with the 
# majority of measurements in the 0-10% range

# Average across sample date: 
cover_2015_by_plot <- cover_2015 %>% group_by(Plot) %>%
  summarize(Treatment = first(Treatment),
            VegAvg = mean(VegAvg),
            LitterAvg = mean(LitterAvg),
            BareGAvg = mean(BareGAvg),
            FWDAvg = mean(FWDAvg),
            CWDAvg = mean(CWDAvg),
            RockAvg = mean(RockAvg),
            VegHtAvg = mean(VegHtAvg),
            canopy_openness = mean(Densi.Total))

# How does vegetation cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=VegAvg)) +
  geom_point(alpha=0.5) 
# Salvaged plots have clearly higher veg cover than forest. Windthrow plots
# vary but some also have higher veg cover

# How does litter cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=LitterAvg)) +
  geom_point(alpha=0.5) 
# Forest clearly with higher (leaf?) litter values

# How does bare ground cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=BareGAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# How does fine woody debris cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=FWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# One forest plot has much more FWD than all other plots

# How does coarse woody debris cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=CWDAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 

# How does rock cover vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=RockAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# One forest plot has much higher rock cover than the rest of the plots

# Vegetation height 2015 #######################################################

# How did vegetation height change between June and July in 2015?
ggplot(data=cover_2015, aes(x=Date, y=VegHtAvg, group=Plot)) +
  geom_point(alpha=0.5) + geom_line(alpha=0.5) 

# How did vegetion height vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=VegHtAvg)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Unclear, but some of the windthrow and salvaged plots appear to have slightly
# higher vegetation height than the forest plots.

# Canopy openness 2015 #########################################################

# How did canopy openness vary between 9-10 June and 5 August, in 2015?
ggplot(data=cover_2015, aes(x=Densi_Date, y=Densi.Total, group=Plot)) +
  geom_point(alpha=0.5) + geom_line(alpha=0.5) # The values seem relatively consistent
# between months

# How did canopy openness vary by treatment in 2015?
ggplot(data=cover_2015_by_plot, aes(x=Treatment, y=canopy_openness)) +
  geom_jitter(alpha=0.5, width=0.05, height=0) 
# Five of six salvaged plots have extremely high canopy openness. One of six
# windthrow plots also has extremely high canopy openness. The rest of the 
# windthrow plots seem to have canopy openness values slightly higher than 
# those of the undisturbed forest.

# Investigating soil moisture, temperature, vegetation cover, veg height,
# and canopy openness in 2015:

# Kayla's pHD found that salvaged plots in 2015 either had high soil moisture
# and temp levels, or lower soil moisture and dense ground-level vegetation

# Join the moisture and temp data to the other data

env_2015_by_plot <- full_join(soil_moisture_temp_2015_by_plot, cover_2015_by_plot, 
                              by="Plot")
env_2015_by_plot$Treatment.x == env_2015_by_plot$Treatment.y

vars_list <- c("mean_moisture", "mean_temp", "canopy_openness", "VegAvg", "VegHtAvg")

# Run a PCA:
pc <- prcomp(env_2015_by_plot[, c(vars_list)], center=T, scale. = T)
attributes(pc)

library(factoextra)
factoextra::get_eig(pc)
factoextra::fviz_pca_biplot(pc, axes=c(1,2))
factoextra::fviz_pca_biplot(pc, axes=c(2,3))

View(pc$x)
pc$rotation

# This needs further study.


