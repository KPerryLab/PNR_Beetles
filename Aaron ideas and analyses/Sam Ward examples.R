
# the data --------------------------------------------------
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
head(ToothGrowth) # data set that is pre-loaded with base R
tail(ToothGrowth)
ToothGrowth$fac_dose <- as.factor(ToothGrowth$dose) # treating dose as a factor, which may or may not be justifiable
summary(ToothGrowth)
table(ToothGrowth$supp,ToothGrowth$dose) # balanced design
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




# randomizing treatments (e.g., 2X3 factorial layout) ------
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Note - this is just an example of how to randomize treatments; we will not 
# be using the dataframe we create in this section

# set up treatments
dose <- c("0.5", "1", "2")
supp <- c("OJ", "VC")
pairwise_combos <- expand.grid(x=dose, y=supp)
treatment_combos <- paste(pairwise_combos$x, pairwise_combos$y, sep=" x ")

# sample size
number_of_observations <- 60 # or number of experimental units


# with six treatments and 60 experimental units, we can have 10 reps
observations <- rep(treatment_combos, 10)
table(observations)


# shuffle the order of observations
set.seed(123) # ensures randomization is the same each time (for reproducibility)
randomized_obs <- sample(x= observations, size=length(observations), replace=F)

# let's say your field (or lab or whatever) layout is 9 rows of 10
matrix(data= randomized_obs, nrow=3, ncol=20)

# let's say your field (or lab or whatever) layout is 6 rows of 15 
matrix(data= randomized_obs, nrow=6, ncol=10)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







# Make a graph and export it -------------------------------
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# get the size the journal wants
# https://academic.oup.com/ee/pages/Manuscript_Preparation#Figures
# max width for single column: 82 mm 3.22835 inches
# max width for double column: 171 mm 6.73228 inches

plot(1:10) # activate plotting window

resize.win <- function(width=6, height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=width, height=height)
}

# for Mac users...
# quartz(width=width, height=height)

resize.win(width=6.73228, height=4)

# source for some code: http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
library(tidyverse)
ggplot(ToothGrowth, aes(x=fac_dose, y=len, fill=supp)) +
  geom_boxplot()+theme_classic()+ 
  theme(legend.position = c(0.8, 0.2)) +
  xlab("Dose level")+
  ylab("Length (units)")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





# Export the graph using R code-------------------------------
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getwd() # figure out R's working directory (and thus where it will save things by default)

# Step 1: Call the pdf command to start the plot
pdf(file = "C:/Users/ward.1792/Desktop/My Plot.pdf",   # The directory you want to save the file in
    width = 6.73228, # The width of the plot in inches
    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code
ggplot(ToothGrowth, aes(x=fac_dose, y=len, fill=supp)) +
  geom_boxplot()+theme_classic()+ 
  theme(legend.position = c(0.8, 0.2)) +
  xlab("Dose level")+
  ylab("Length (units)")

# Step 3: Run dev.off() to create the file!
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







# get stats from R into MS Word -----------------------------
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get means and standard errors for writing results
library(plotrix)
ToothGrowth %>%
  group_by(supp, fac_dose) %>%
  summarise(
    means = round(mean(len),2),
    se = round(std.error(len),2))

ToothGrowth %>%
  group_by(supp) %>%
  summarise(
    means = round(mean(len),2),
    se = round(std.error(len),2))


# get stats for a table
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}


fit1 <- lm(len ~ fac_dose + supp, data=ToothGrowth)
summary(fit1)

round(summary(fit1)$coefficients, 2)
round(summary(fit1)$coefficients, 4)

write.excel(round(summary(fit1)$coefficients, 2))
write.excel(round(summary(fit1)$coefficients, 4))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




