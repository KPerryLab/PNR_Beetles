# Import weather data from near Powdermill in order to understand precipitation
# and temperature patterns

# Aaron Tayal
# March 11, 2025

library(ggplot2)
library(lubridate)
library(dplyr)

# Import data
weather <- read.csv("PNR_Raw_Data/DONEGAL 2 NW, PA US (USC00362183).csv")
head(weather)

which(weather$Date == "2015-01-01") # Jan 1st, 2015
which(weather$Date == "2015-05-27") # Traps installed May 27-28, 2015
which(weather$Date == "2015-08-17") # Traps taken down Aug 17, 2015

# Convert character strings to dates:
weather$Date1 <- ymd(weather$Date)
weather$Year <- year(weather$Date1)

# Generate a column for cumulative precipitation starting from January 1st
weather$PRCP_cumulative_inches <- rep(0, length(weather[,"Date"]))
inches <- 0
for (day in 1:length(weather[,"Date"])){
  todays_date <- weather[day, "Date1"]
  if (month(todays_date) == 1 & day(todays_date) == 1) {
    inches <- 0
  }
  if (!is.na(weather[day, "PRCP..Inches."])) {
    inches <- inches + weather[day, "PRCP..Inches."]
  }
  weather$PRCP_cumulative_inches[day] <- inches
}

# Create a subset for 2015 weather
weather_2015 <- weather %>% filter(Year == 2015)

# Create a subset for 2015 weather (just from May 27 to Aug 17)
weather_2015_cutoff <- weather %>% filter(Date1 >= ymd("2015-05-27") & 
                                     Date1 < ymd("2015-08-17"))

# Plot temperature for 2015 (just from May 27 to Aug 17):
ggplot(data=weather_2015_cutoff, aes(x=Date1)) +
  geom_line(aes(y=TMAX..Degrees.Fahrenheit.), color="red") +
  geom_line(aes(y=TMIN..Degrees.Fahrenheit.), color="green")

# Plot temperature for 2015:
ggplot(data=weather_2015, aes(x=Date1)) +
  geom_line(aes(y=TMAX..Degrees.Fahrenheit.), color="red") +
  geom_line(aes(y=TMIN..Degrees.Fahrenheit.), color="green")

# Plot precipitation for 2015 (clipped at Aug 17):
ggplot(data=weather_2015_cutoff, aes(x=Date1)) +
  geom_line(aes(y=PRCP..Inches.), color="blue") +
  geom_line(aes(y=PRCP_cumulative_inches))

# Plot precipitation for 2015:
ggplot(data=weather_2015, aes(x=Date1)) +
  geom_line(aes(y=PRCP..Inches.), color="blue") +
  geom_line(aes(y=PRCP_cumulative_inches))

# Plot snowfall in 2015:
ggplot(data=weather_2015, aes(x=Date1)) +
  geom_line(aes(y=SNOW..Inches.)) +
  geom_line(aes(y=SNWD..Inches.), color="purple") # Snow depth

# Create a subset for 2022 weather
weather_2022 <- weather %>% filter(Year == 2022)

# Create a subset for 2022 weather (just from June 1 to Sept 6)
weather_2022_cutoff <- weather %>% filter(Date1 >= ymd("2022-06-01") & 
                                            Date1 < ymd("2022-09-6"))

# Plot temperature for 2022:
ggplot(data=weather_2022, aes(x=Date1)) +
  geom_line(aes(y=TMAX..Degrees.Fahrenheit.), color="red") +
  geom_line(aes(y=TMIN..Degrees.Fahrenheit.), color="green")

# Plot temperature for 2022 (clipped at Sept 6)
ggplot(data=weather_2022_cutoff, aes(x=Date1)) +
  geom_line(aes(y=TMAX..Degrees.Fahrenheit.), color="red") +
  geom_line(aes(y=TMIN..Degrees.Fahrenheit.), color="green")

# Plot precipitation for 2022:
ggplot(data=weather_2022, aes(x=Date1)) +
  geom_line(aes(y=PRCP..Inches.), color="blue") +
  geom_line(aes(y=PRCP_cumulative_inches))

# Plot precipitation for 2022 (clipped at Sept 6):
ggplot(data=weather_2022_cutoff, aes(x=Date1)) +
  geom_line(aes(y=PRCP..Inches.), color="blue") +
  geom_line(aes(y=PRCP_cumulative_inches))

# Plot snowfall in 2022:
ggplot(data=weather_2022, aes(x=Date1)) +
  geom_line(aes(y=SNOW..Inches.)) +
  geom_line(aes(y=SNWD..Inches.), color="purple") # Snow depth

# Create a subset of weather from 2012 through 2022:
weather_2012_2022 <- weather %>% filter(Date1 >= ymd("2012-01-01") & 
                     Date1 < ymd("2023-01-01"))

# plot precipitation from 2012 through 2022
ggplot(data=weather_2012_2022, aes(x=Date1)) +
  geom_line(aes(y=PRCP..Inches.), color="blue") +
  geom_line(aes(y=PRCP_cumulative_inches))

# Create a subset of Dec 31st of each year from 2012 to 2023:
weather_dec_31st <- weather %>% filter(Year>2011 & Year<2024 & day(Date1)==31
                                       & month(Date1)==12)

# now summarize the annual precipitation:
summary(weather_dec_31st$PRCP_cumulative_inches)



