# Aaron Tayal
# 4/30/2025
# Carabids stats

# The purpose is to create linear mixed-effects models to investigate
# metrics of ground beetle activity-abundance, species richness, Shannon
# diversity, functional alpha diversity, and community-weighted means, and 
# whether these variables differ between plots that were salvaged, 
# just windthow-affected, or forest control. And whether they differ between 
# 2015 and 2022.

# I'll also look at some environmental variables such as soil moisture,
# percent vegetation cover, and canopy openness

library(ggplot2)
theme_set(theme_classic())
library(lme4) # linear mixed-effects models
library(lmerTest) # uses Satterthwaites approx. for degrees of freedom
library(DHARMa) # Can help test for overdispersion in a Poisson GLMM
library(car) # Used to get an ANOVA table for the Poisson GLMM
library(emmeans) # Estimated Marginal Means, aka Least Squares Means (enables 
# pairwise comparisons to be made)
library(ggbeeswarm) # an alternative to geom_jitter is geom_quasirandom()
library(ggpubr) # used for making pretty graphs
library(sjPlot) # used to plot models
library(gridExtra) # used for arranging graphs properly
library(plotrix) # for standard error function
library(dplyr)

dat <- read.csv("Aaron_PNR_formatted_data/PNR2015_2022_carabid_counts_by_plot_standardized.csv")

dat$Treatment <- as.factor(dat$Treatment) # The forest management treatment

dat$Year <- as.factor(dat$Year) # The year of sampling (need to treat it as 
# a factor variable)

dat$Transect <- as.factor(dat$Transect) # The random effect

dat$Year_Treatment <- interaction(dat$Year, dat$Treatment)

# Total activity-abundance models ###############################################

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=total_count_stdz)) + geom_quasirandom(width=0.1, alpha=0.5)

# take the log transform:
dat$log_total_count_stdz <- log(dat$total_count_stdz)

ggplot(dat, aes(x=Year_Treatment, y=log_total_count_stdz)) + geom_quasirandom(width=0.1, alpha=0.5)

# run the model (on the log-transformed data):
model_total <- lmerTest::lmer(log_total_count_stdz ~ Treatment + Year + 
                                Treatment*Year + (1|Transect), data = dat)
summary(model_total)

# plot the model:
plot_model(model_total, type = "pred", terms = c("Treatment", "Year")) +
  theme_classic()

# test the assumptions:
plot(model_total)
qqnorm(residuals(model_total))
qqline(residuals(model_total))
ranef_brood <- ranef(model_total)$Transect
hist(ranef_brood$`(Intercept)`, breaks=10)

# run the ANOVA test (type 3):
anova(model_total, type=3)

# run pairwise tests for treatment (simple effects, because the interaction term was significant)
emmeans(model_total, pairwise ~ Treatment | Year) # compares levels of Treatment within each level of Year

# Activity-abundance of open-habitat and eurytopic species models ###############

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=open_habitat_spp_stdz + eurytopic_spp_stdz)) + 
  geom_quasirandom(width=0.1, alpha=0.5)

# take the log transform (need to add 1 because one data point is still 0)
dat$log_oe <- log(dat$open_habitat_spp_stdz + dat$eurytopic_spp_stdz + 1)

ggplot(dat, aes(x=Year_Treatment, y=log_oe)) + geom_quasirandom(width=0.1, alpha=0.5)

# run the model (on the log(x+1) transformed data):
model_oe <- lmerTest::lmer(log_oe ~ Treatment + Year + 
                             Treatment*Year + (1|Transect), data = dat)
summary(model_oe)

# plot the model:
plot_model(model_oe, type = "pred", terms = c("Treatment", "Year")) +
  theme_classic()

# test assumptions:
plot(model_oe)
qqnorm(residuals(model_oe))
qqline(residuals(model_oe))

# run the ANOVA test:
anova(model_oe, type = 3)

# post-hoc comparisons between treatment levels:
emmeans(model_oe, pairwise~Treatment)

# Activity-abundance of forest species models ##################################

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=forest_specialist_spp_stdz)) + 
  geom_quasirandom(width=0.1, alpha=0.5)

# take the log transformation:
dat$log_f <- log(dat$forest_specialist_spp_stdz)

ggplot(dat, aes(x=Year_Treatment, y=log_f)) + geom_quasirandom(width=0.1, alpha=0.5)

# run the model:
model_f <- lmerTest::lmer(log_f ~ Treatment + Year + 
                             Treatment*Year + (1|Transect), data = dat)
summary(model_f)

# plot the model:
plot_model(model_f, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_f)
qqnorm(residuals(model_f))
qqline(residuals(model_f))

# run the ANOVA test:
anova(model_f, type = 3)

# Species richness models ######################################################

# See what the averages are:
mean(dat$sp_rich)

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=sp_rich)) + geom_quasirandom(width=0.1, alpha=0.5)

# Because I have a count response variable (the number of species), I need to
# use a Poisson response variable:

# The model that includes Transect gets a singular fit, so I need to run a regular
# GLM:
model_rich <- glm(sp_rich ~ Treatment + Year + Treatment*Year, data=dat, family="poisson")
summary(model_rich)

# plot the model:
plot_model(model_rich, type = "pred", terms = c("Treatment", "Year"))

# check assumptions:
#plot(model_rich)

# Run the ANOVA test (need to use the car package because it is a GLM):
car::Anova(model_rich, type="III")

# Post-hoc comparisons:
emmeans(model_rich, pairwise~Treatment)

# Shannon diversity models #####################################################

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=shannon_diversity)) + geom_quasirandom(width=0.1, alpha=0.5)

# run the model:
model_shannon <- lmerTest::lmer(shannon_diversity ~ Treatment + Year + 
                                  Treatment*Year + (1|Transect), data = dat)
summary(model_shannon)

# test assumptions:
plot(model_shannon)
qqnorm(residuals(model_shannon))
qqline(residuals(model_shannon))

# run the ANOVA test:
anova(model_shannon, type=3)

# Functional alpha diversity models ###########################################

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=mean_pairwise_distance)) + geom_quasirandom(width=0.1, alpha=0.5)

# run the model (running without Transect because it was a singular fit):
model_mpd <- lm(mean_pairwise_distance ~ Treatment + Year + 
                              Treatment*Year, data = dat)
summary(model_mpd)

# plot the model:
plot_model(model_mpd, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
#plot(model_mpd)

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_mpd)

# Community-weighted mean trait models ########################################

# PC1 #########################################################################

# First I'll do PC1. This variable is a trait syndrome associated with
# longer limbs (longer legs and probably longer antennae), shorter trochanter
# which is associated with less powerful pushing from legs, and a narrower
# pronotum, possibly to reach into small crevices.

# graph the data:
ggplot(dat, aes(x=Year_Treatment, y=PC1)) + geom_quasirandom(alpha=0.5, width=0.1)

# run the model:
model_PC1 <- lmerTest::lmer(PC1 ~ Treatment + Year + Treatment*Year +
                              (1|Transect), data = dat)
summary(model_PC1)

# plot the model:
plot_model(model_PC1, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_PC1)
qqnorm(residuals(model_PC1))
qqline(residuals(model_PC1))

# run the ANOVA test:
anova(model_PC1, type=3)

# run the post-hoc test for pairwise comparisons:
emmeans(model_PC1, pairwise~Treatment)

# Higher PC1 is associated with proportionally narrower pronotum, proportionally 
# longer rear legs, and proportionally shorter rear trochanter

# PC2 #########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=PC2)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_PC2 <- lmerTest::lmer(PC2 ~ Treatment + Year + Treatment*Year +
                            (1|Transect), data = dat)
summary(model_PC2)

# plot the model:
plot_model(model_PC2, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_PC2)
qqnorm(residuals(model_PC2))
qqline(residuals(model_PC2))

# run the ANOVA test:
anova(model_PC2, type=3)

# Now I need to do pairwise comparisons:
emmeans(model_PC2, pairwise ~ Treatment)

# PC3 ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=PC3)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_PC3 <- lmerTest::lmer(PC3 ~ Treatment + Year + Treatment*Year +
                              (1|Transect), data = dat)
summary(model_PC3)

# plot the model:
plot_model(model_PC3, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_PC3)
qqnorm(residuals(model_PC3))
qqline(residuals(model_PC3))

# run the ANOVA test:
anova(model_PC3, type=3)

# Now I need to do pairwise comparisons:
emmeans(model_PC3, pairwise ~ Treatment)

# Individual functional traits #################################################

# The traits are: body length, stdzd antenna length, antenna:rear leg ratio, stdzd rear leg length,
# stdzd eye length, stdzd eye protrusion, eye protrusion:eye length ratio, 
# stdzd pronotum width, stdzd abdomen width,  stdzd rear trochanter length, 
# flight capability, and water affinity

# Body length ##################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=body_length)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create model using linear model (LMM was singular fit)
model_body_length <- lm(body_length ~ Treatment + Year + Treatment*Year,
                              data = dat)
summary(model_body_length)

# plot the model:
plot_model(model_body_length, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
#plot(model_body_length)
qqnorm(residuals(model_body_length))
qqline(residuals(model_body_length))

# run the ANOVA test:
anova(model_body_length)

# Now I need to do pairwise comparisons:
emmeans(model_body_length, pairwise ~ Treatment)

# antenna_length_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=antenna_length_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_antenna_length_standard <- lmerTest::lmer(antenna_length_standard ~ Treatment + Year + Treatment*Year +
                              (1|Transect), data = dat)
summary(model_antenna_length_standard)

# plot the model:
plot_model(model_antenna_length_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_antenna_length_standard)
qqnorm(residuals(model_antenna_length_standard))
qqline(residuals(model_antenna_length_standard))

# run the ANOVA test:
anova(model_antenna_length_standard, type=3)

# rear_leg_length_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=rear_leg_length_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_rear_leg_length_standard <- lmerTest::lmer(rear_leg_length_standard ~ Treatment + Year + Treatment*Year +
                                                  (1|Transect), data = dat)
summary(model_rear_leg_length_standard)

# plot the model:
plot_model(model_rear_leg_length_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_rear_leg_length_standard)
qqnorm(residuals(model_rear_leg_length_standard))
qqline(residuals(model_rear_leg_length_standard))

# run the ANOVA test:
anova(model_rear_leg_length_standard, type=3)

# antenna_rear_leg_ratio ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=antenna_rear_leg_ratio)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear model (the LMM had a singular fit):
model_antenna_rear_leg_ratio <- lm(antenna_rear_leg_ratio ~ Treatment + Year + Treatment*Year, 
                                   data = dat)
summary(model_antenna_rear_leg_ratio)

# plot the model:
plot_model(model_antenna_rear_leg_ratio, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
#plot(model_antenna_rear_leg_ratio)
qqnorm(residuals(model_antenna_rear_leg_ratio))
qqline(residuals(model_antenna_rear_leg_ratio))

# run the ANOVA test:
anova(model_antenna_rear_leg_ratio)

# eye_length_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=eye_length_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_eye_length_standard <- lmerTest::lmer(eye_length_standard ~ Treatment + Year + Treatment*Year +
                                                 (1|Transect), data = dat)
summary(model_eye_length_standard)

# plot the model:
plot_model(model_eye_length_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_eye_length_standard)
qqnorm(residuals(model_eye_length_standard))
qqline(residuals(model_eye_length_standard))

# run the ANOVA test:
anova(model_eye_length_standard, type=3)

# pairwise comparisons:
emmeans(model_eye_length_standard, pairwise~Treatment)

# eye_protrusion_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=eye_protrusion_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_eye_protrusion_standard <- lmerTest::lmer(eye_protrusion_standard ~ Treatment + Year + Treatment*Year +
                                               (1|Transect), data = dat)
summary(model_eye_protrusion_standard)

# plot the model:
plot_model(model_eye_protrusion_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_eye_protrusion_standard)
qqnorm(residuals(model_eye_protrusion_standard))
qqline(residuals(model_eye_protrusion_standard))

# run the ANOVA test:
anova(model_eye_protrusion_standard, type=3)

# eye_protrusion_ratio ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=eye_protrusion_ratio)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_eye_protrusion_ratio <- lmerTest::lmer(eye_protrusion_ratio ~ Treatment + Year + Treatment*Year +
                                                  (1|Transect), data = dat)
summary(model_eye_protrusion_ratio)

# plot the model:
plot_model(model_eye_protrusion_ratio, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_eye_protrusion_ratio)
qqnorm(residuals(model_eye_protrusion_ratio))
qqline(residuals(model_eye_protrusion_ratio))

# run the ANOVA test:
anova(model_eye_protrusion_ratio, type=3)

# Test pairwise comparisons:
emmeans(model_eye_protrusion_ratio, pairwise~Treatment)

# pronotum_width_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=pronotum_width_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_pronotum_width_standard <- lmerTest::lmer(pronotum_width_standard ~ Treatment + Year + Treatment*Year +
                                              (1|Transect), data = dat)
summary(model_pronotum_width_standard)

# plot the model:
plot_model(model_pronotum_width_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_pronotum_width_standard)
qqnorm(residuals(model_pronotum_width_standard))
qqline(residuals(model_pronotum_width_standard))

# run the ANOVA test:
anova(model_pronotum_width_standard, type=3)

# abdomen_width_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=abdomen_width_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_abdomen_width_standard <- lmerTest::lmer(abdomen_width_standard ~ Treatment + Year + Treatment*Year +
                                                  (1|Transect), data = dat)
summary(model_abdomen_width_standard)

# plot the model:
plot_model(model_abdomen_width_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_abdomen_width_standard)
qqnorm(residuals(model_abdomen_width_standard))
qqline(residuals(model_abdomen_width_standard))

# run the ANOVA test:
anova(model_abdomen_width_standard, type=3)

# rear_trochanter_length_standard ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=rear_trochanter_length_standard)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear mixed-effects models:
model_rear_trochanter_length_standard <- lmerTest::lmer(rear_trochanter_length_standard ~ Treatment + Year + Treatment*Year +
                                                 (1|Transect), data = dat)
summary(model_rear_trochanter_length_standard)

# plot the model:
plot_model(model_rear_trochanter_length_standard, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_rear_trochanter_length_standard)
qqnorm(residuals(model_rear_trochanter_length_standard))
qqline(residuals(model_rear_trochanter_length_standard))

# run the ANOVA test:
anova(model_rear_trochanter_length_standard, type=3)

# pairwise comparisons:
emmeans(model_rear_trochanter_length_standard, pairwise~Treatment)

# Water_affinity ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=Water_affinity)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now create models using linear model (LMM was singular fit):
model_Water_affinity <- lm(Water_affinity ~ Treatment + Year + Treatment*Year, data = dat)
summary(model_Water_affinity)

# plot the model:
plot_model(model_Water_affinity, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
#plot(model_Water_affinity)
qqnorm(residuals(model_Water_affinity))
qqline(residuals(model_Water_affinity))

# run the ANOVA test:
anova(model_Water_affinity)

# Flight_capability ##########################################################################

# Graph the data:
ggplot(dat, aes(x=Year_Treatment, y=Flight_capability)) + geom_quasirandom(alpha=0.5, width=0.1)

# take the log(x+0.1) transform:
dat$log_Flight_capability <- log(dat$Flight_capability + 0.1)

ggplot(dat, aes(x=Year_Treatment, y=log_Flight_capability)) + geom_quasirandom(alpha=0.5, width=0.1)

# I'll now run the LMM model:
model_Flight_capability <- lmerTest::lmer(log_Flight_capability ~ Treatment + Year + 
                                            Treatment*Year + (1|Transect), data = dat)
summary(model_Flight_capability)

# plot the model:
plot_model(model_Flight_capability, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_Flight_capability)
qqnorm(residuals(model_Flight_capability))
qqline(residuals(model_Flight_capability))

# run the ANOVA test:
anova(model_Flight_capability, type = 3)

# Now test for pairwise differences within each year:
emmeans(model_Flight_capability, pairwise ~ Treatment | Year)

# Environmental variables: soil moisture #######################################

# Soil moisture must be analyzed separately for each year because a different
# sensor was used

# Import the environmental data:
env <- read.csv("Aaron_PNR_formatted_data/PNR2015_2022_environment_by_plot.csv")

# Change the following columns to factors:
env$Transect <- as.factor(env$Transect)
env$Treatment <- as.factor(env$Treatment)
env$Year <- as.factor(env$Year)

# Create a column for the Year Treatment concatenated variable
env$Year_Treatment <- interaction(env$Year, env$Treatment)

# Graph the data (2015):
ggplot(env %>% filter(Year==2015), aes(x=Treatment, y=mean_moisture)) + 
  geom_quasirandom(alpha=0.5, width=0.1)

# Graph the data (2022):
ggplot(env %>% filter(Year==2022), aes(x=Treatment, y=mean_moisture)) + 
  geom_quasirandom(alpha=0.5, width=0.1)

# Model (2015):
model_2015_soil_moisture <- 
  lmerTest::lmer(mean_moisture ~ Treatment + (1|Transect), data = env %>% 
                   filter(Year==2015))
summary(model_2015_soil_moisture)

# Plot the model (2015):
plot_model(model_2015_soil_moisture, type="pred", terms=c("Treatment"))

# Anova test (2015)
anova(model_2015_soil_moisture)

# Model (2022):
model_2022_soil_moisture <- 
  lmerTest::lmer(mean_moisture ~ Treatment + (1|Transect), data = env %>% 
                   filter(Year==2022))
summary(model_2022_soil_moisture)

# Plot the model (2022):
plot_model(model_2022_soil_moisture, type="pred", terms=c("Treatment"))

# Anova test (2022)
anova(model_2022_soil_moisture)

# Vegetation percentage cover ##################################################

# Graph the data:
ggplot(data=env, aes(x=Year_Treatment, y=VegAvg)) + 
  geom_quasirandom(alpha=0.5, width=0.1) + ylab("Mean % cover of vegetation (%)")

# Run the linear mixed-effects model:
model_VegAvg <- lmerTest::lmer(VegAvg ~ Treatment + Year + Treatment*Year +
                                 (1|Transect), data=env)
summary(model_VegAvg)

# Plot the model:
plot_model(model_VegAvg, type="pred", terms=c("Treatment", "Year"))

# Check assumptions:
plot(model_VegAvg)
qqnorm(residuals(model_VegAvg))
qqline(residuals(model_VegAvg))
hist(residuals(model_VegAvg), breaks=15) # Residuals have a long
# tail on the right side, but a short tail on the left side.

# Run the ANOVA test (type 3):
anova(model_VegAvg, type=3)

# Do pairwise comparisons by looking for simple effects (because 
# the Treatment*Year interaction was significant)
emmeans(model_VegAvg, pairwise ~ Treatment | Year)
emmeans(model_VegAvg, pairwise ~ Year | Treatment)

# leaf litter percentage cover ##################################################

# Graph the data:
ggplot(data=env, aes(x=Year_Treatment, y=LitterAvg)) + 
  geom_quasirandom(alpha=0.5, width=0.1) + ylab("Mean % cover of leaf litter (%)")

# Run the linear mixed-effects model:
model_LitterAvg <- lmerTest::lmer(LitterAvg ~ Treatment + Year + Treatment*Year +
                                 (1|Transect), data=env)
summary(model_LitterAvg)

# Plot the model:
plot_model(model_LitterAvg, type="pred", terms=c("Treatment", "Year"))

# Check assumptions:
plot(model_LitterAvg)
qqnorm(residuals(model_LitterAvg))
qqline(residuals(model_LitterAvg))
hist(residuals(model_LitterAvg), breaks=15) # Residuals have a long
# tail on the right side, but a short tail on the left side.

# Run the ANOVA test (type 3):
anova(model_LitterAvg, type=3)

# Do pairwise comparisons by looking for simple effects (because 
# the Treatment*Year interaction was significant)
emmeans(model_LitterAvg, pairwise ~ Treatment | Year)
emmeans(model_LitterAvg, pairwise ~ Year | Treatment)

# Canopy openness #############################################################

# Graph the data:
ggplot(data=env, aes(x=Year_Treatment, y=Densi.Total)) + 
  geom_quasirandom(alpha=0.5, width=0.1) + ylab("Canopy openness (%)")# + geom_text(aes(label=Plot))

# Take the logarithm to improve the homogeniety of variances assumption:
env$log_Densi.Total = log(env$Densi.Total)

# Graph again:
ggplot(data=env, aes(x=Year_Treatment, y=log_Densi.Total)) + 
  geom_quasirandom(alpha=0.5, width=0.1) + ylab("Canopy openness (%)")

# Run the model (linear mixed effects):
model_canopy_openness <- lmer(log_Densi.Total ~ Treatment + Year + 
                                          Treatment*Year + (1|Transect), data=env)
summary(model_canopy_openness)

# Plot the model:
plot_model(model_canopy_openness, type="pred", terms=c("Treatment", "Year"))

# Test assumptions:
plot(model_canopy_openness)
qqnorm(residuals(model_canopy_openness))
qqline(residuals(model_canopy_openness))

# ANOVA test:
anova(model_canopy_openness, type=3)

# pairwise comparisons (simple effects because the interaction term was significant)
emmeans(model_canopy_openness, pairwise ~ Treatment | Year)
emmeans(model_canopy_openness, pairwise ~ Year | Treatment)

# Pretty graphs ###############################################################

dat$Treatment <- factor(dat$Treatment, levels = c("Windthrow", "Salvaged", "Forest"))

abundance_graph <- ggplot(dat, aes(x = Treatment, y = total_count_stdz, fill = Year, group = Year)) + 
  stat_summary(fun = mean, geom = "bar", color="black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9)) +
  ylab("Number of individuals / plot") + theme(plot.title = element_text(size=18),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_text(size = 16, 
                                                                            margin = margin(r=20)),
                                                axis.text.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 14),
                                               legend.title = element_text(size = 14),
                                               legend.text = element_text(size = 14))+
  scale_fill_grey()
abundance_graph

richness_graph <- ggplot(dat, aes(x=Treatment, y=sp_rich, fill=Year, group=Year)) + 
  stat_summary(fun = mean, geom = "bar", color="black", position=position_dodge(width=0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position=position_dodge(width=0.9)) +
  ylab("Number of species / plot") + theme(plot.title = element_text(size=18),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(size = 16, 
                                                                margin = margin(r=15)),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14),
                                    legend.title = element_text(size = 14),
                                    legend.text = element_text(size = 14))+
  scale_fill_grey()
richness_graph

empty_graph <- ggplot() + theme_void()

ggarrange(abundance_graph, empty_graph, richness_graph,
          labels = c("A", "", "B"), ncol=3, nrow=1, widths = c(1, 0.1, 1), legend = "right")

# legend.position = "none"

# Make summary data tables for treatment means: ################################

response_vars <- c("total_count_stdz", "open_habitat_spp_stdz", 
                   "eurytopic_spp_stdz", "forest_specialist_spp_stdz",
                   "sp_rich", "shannon_diversity", "mean_pairwise_distance",
                   "PC1", "PC2", "PC3", "body_length","antenna_length_standard",
                   "rear_leg_length_standard", "antenna_rear_leg_ratio",
                   "eye_length_standard", "eye_protrusion_standard",
                   "eye_protrusion_ratio", "pronotum_width_standard", 
                   "abdomen_width_standard", "rear_trochanter_length_standard", 
                   "Water_affinity", "Flight_capability")

mean_concat_std_error <- function(x) {
  return(paste(round(mean(x), 3), "+-", round(std.error(x), 3)))
}

dat_by_treatment <- dat %>% group_by(Year, Treatment) %>%
  summarize(across(all_of(response_vars), ~ mean_concat_std_error(.)))

#write.csv(dat_by_treatment, "Aaron_PNR_formatted_data/PNR2015_2022_response_vars_by_treatment.csv", row.names = F)

env_vars <- c("mean_moisture", "VegAvg", "LitterAvg", "Densi.Total")

env_by_treatment <- env %>% group_by(Year, Treatment) %>%
  summarize(across(all_of(env_vars), ~ mean_concat_std_error(.)))

#write.csv(env_by_treatment, "Aaron_PNR_formatted_data/PNR2015_2022_environmental_vars_by_treatment.csv", row.names = F)


