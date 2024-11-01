# Aaron Tayal
# Sept 24, 2024
# Exploratory data visualization of environmental variables recorded at 
# Powdermill Nature Reserve in the tornado/salvage areas

library(ggplot2)
library(dplyr)
library(lubridate)

soil_moisture_0 <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_Soil_moisture_PNR_ENV_2022.csv")
soil_moisture_0$SoilMoisture_day <- as.factor(yday(soil_moisture_0$SoilMoisture_Date))

# Information on the transect groups and areas of each trap
trap_locations <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_PNR_PitfallTrapLocations_2015.csv",
                           colClasses = c("integer", "integer", "factor", "factor", "factor", "numeric", "numeric"))
soil_moisture <- right_join(trap_locations, soil_moisture_0, by="Plot")

summary(soil_moisture)
plot(soil_moisture$percent_Moisture_Avg)
hist(soil_moisture$percent_Moisture_Avg)
hist(soil_moisture$percent_MOISTURE1)
hist(soil_moisture$percent_MOISTURE2)
hist(soil_moisture$percent_MOISTURE3)

ggplot(data=soil_moisture, aes(x=SoilMoisture_day, 
                               y=percent_Moisture_Avg)) +
  geom_point()

# Idea: make a Parallel coordinates plot

# How much variation in soil moisture is related to forest disturbance
# type (undisturbed, wind-disturbed, or salvage-logged)? How much variation
# is related to transect? How much variation is related to date of measurement?

# Make sure I can calculate variance:
example <- c(1,4,9,7,6)
var(example)
mean(example)
squared_diffs <- vapply(example, FUN = function(x) (x - mean(example))^2,
                        FUN.VALUE = 1)
sum(squared_diffs) / (length(example) - 1)

# Now calculate the variance of percent soil moisture (averaging the three readings):
overall_var <- var(soil_moisture$percent_Moisture_Avg)
overall_mean <- mean(soil_moisture$percent_Moisture_Avg)
x_values <- -25:125
# Prob. density function of normal distribution:
y_values <- dnorm(x_values, overall_mean, sqrt(overall_var))
histogram <- hist(soil_moisture$percent_Moisture_Avg, plot=FALSE, breaks=20)
plot(histogram)
lines(x=x_values, y=1000*y_values)
qqline(soil_moisture$percent_Moisture_Avg)

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






