# Aaron Tayal
# Sept 24, 2024
# Exploratory data visualization of environmental variables recorded at 
# Powdermill Nature Reserve in the tornado/salvage areas

library(ggplot2)
library(dplyr)
library(lubridate)

soil_moisture <- read.csv("Aaron_PNR_formatted_csvs/Aaron_formatted_Soil_moisture_PNR_ENV_2022.csv")
soil_moisture$SoilMoisture_day <- yday(soil_moisture$SoilMoisture_Date)


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
