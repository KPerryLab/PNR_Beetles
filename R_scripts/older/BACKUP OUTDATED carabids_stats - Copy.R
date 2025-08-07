# Aaron Tayal
# 4/30/2025
# Carabids stats

# The purpose is to create linear mixed-effects models to investigate
# metrics of ground beetle activity-abundance, species richness, Shannon
# diversity, and functional alpha diversity, and whether these variables
# differ between plots that were salvaged, windthow-affected, or forest control.

# I'm going to create separate models for 2015 and 2022.

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

dat_2015 <- read.csv("Aaron_PNR_formatted_data/PNR2015_carabid_counts_by_plot_standardized.csv")

dat_2022 <- read.csv("Aaron_PNR_formatted_data/PNR2022_carabid_counts_by_plot_standardized.csv")

dat_2015$Treatment <- as.factor(dat_2015$Treatment) # The predictor variable
dat_2022$Treatment <- as.factor(dat_2022$Treatment) 

dat_2015$Transect <- as.factor(dat_2015$Transect) # The response variable
dat_2022$Transect <- as.factor(dat_2022$Transect)

# Total activity-abundance models ###############################################

# graph the data:
ggplot(dat_2015, aes(x=Treatment, y=total_count_stdz, color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=total_count_stdz, color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_total <- lmerTest::lmer(total_count_stdz ~ Treatment + (1|Transect), 
                               data = dat_2015)
summary(model_2015_total)

model_2022_total <- lmerTest::lmer(total_count_stdz ~ Treatment + (1|Transect), 
                               data = dat_2022)
summary(model_2022_total)

anova(model_2015_total, type = 3)
anova(model_2022_total, type = 3)

plot(model_2015_total) # I'm seeing some heteroscedasticity
plot(model_2022_total) # The 2022 model seems OK in terms of variance of residuals

qqnorm(residuals(model_2015_total))
qqline(residuals(model_2015_total)) # Appears to be some over-dispersion -
# non-normality of residuals

qqnorm(residuals(model_2022_total))
qqline(residuals(model_2022_total))

ranef_brood_2015 <- ranef(model_2015_total)$Transect
hist(ranef_brood_2015$`(Intercept)`, breaks=10) # The estimated random intercepts
# do NOT seem to be normally distributed

ranef_brood_2022 <- ranef(model_2022_total)$Transect
hist(ranef_brood_2022$`(Intercept)`, breaks=10)

# A log transformation could help to address any heteroscedasticity.

dat_2015$log_total_count_stdz <- log(dat_2015$total_count_stdz)

# plot the log-transformed data:
ggplot(dat_2015, aes(x=Treatment, y=log_total_count_stdz, color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Now make a new model using the log-transformed response variable:
model_2015_total_2 <- lmerTest::lmer(log_total_count_stdz ~ Treatment + (1|Transect), 
                                   data = dat_2015)
summary(model_2015_total_2)
plot(model_2015_total_2) # There is no heterscedasticity, which is good
qqnorm(residuals(model_2015_total))
qqline(residuals(model_2015_total))
anova(model_2015_total_2, type=3)
car::Anova(model_2015_total_2, type = "III") # I'm getting a different p-value 
# from the Analysis of Deviance Table, versus the Analysis of Variance table

# I'm going to go with the value from the anova() command

# Now I need to run pairwise comparisons to see if there are differences
# between treatment groups in activity-abundance
emmeans(model_2015_total_2, pairwise~Treatment)

# Alternatively, I could try a Kruskal-Wallis test (although I can't include
# the random effect of transect in this one):
kruskal.test(total_count_stdz ~ Treatment, data=dat_2015)

# Activity-abundance of open-habitat and eurytopic species models ###############

# graph the data:
ggplot(dat_2015, aes(x=Treatment, y=open_habitat_spp_stdz + eurytopic_spp_stdz, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=open_habitat_spp_stdz + eurytopic_spp_stdz, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Create the linear mixed-effects models:
model_2015_oe <- lmerTest::lmer(open_habitat_spp_stdz + eurytopic_spp_stdz ~ Treatment + (1|Transect), 
                                   data = dat_2015)
summary(model_2015_oe)

model_2022_oe <- lmerTest::lmer(open_habitat_spp_stdz + eurytopic_spp_stdz ~ Treatment + (1|Transect), 
                                   data = dat_2022)
summary(model_2022_oe)

anova(model_2015_oe, type = 3)
anova(model_2022_oe, type = 3)

plot(model_2015_oe) # I'm seeing some heteroscedasticity
plot(model_2022_oe) # Also seeing some heteroscedasticity

qqnorm(residuals(model_2015_oe))
qqline(residuals(model_2015_oe)) # Appears to be a small amount of overdispersion

qqnorm(residuals(model_2022_oe))
qqline(residuals(model_2022_oe))

ranef_brood_2015_oe <- ranef(model_2015_oe)$Transect
hist(ranef_brood_2015_oe$`(Intercept)`, breaks=10) # The estimated random intercepts
# seem like they could be normally distributed

ranef_brood_2022_oe <- ranef(model_2022_oe)$Transect
hist(ranef_brood_2022_oe$`(Intercept)`, breaks=10) # not normally distributed

# To remedy some of the heteroscedasticity, I'll take the log. Because there
# are a couple of zeros in the 2015 data, I'll need to add 0.1 to the 2015 data before 
# taking the log, in order to prevent taking log(0):

dat_2015$log_open_eurytopic <- log(dat_2015$open_habitat_spp_stdz + 
                                     dat_2015$eurytopic_spp_stdz + 0.1)

dat_2022$log_open_eurytopic <- log(dat_2022$open_habitat_spp_stdz + 
                                     dat_2022$eurytopic_spp_stdz)

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=log_open_eurytopic, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=log_open_eurytopic, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Create the linear mixed-effects models:
model_2015_oe_2 <- lmerTest::lmer(log_open_eurytopic ~ Treatment + (1|Transect), 
                                data = dat_2015)
summary(model_2015_oe_2)

model_2022_oe_2 <- lmerTest::lmer(log_open_eurytopic ~ Treatment + (1|Transect), 
                                data = dat_2022)
summary(model_2022_oe_2)

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_2015_oe_2, type = 3)
anova(model_2022_oe_2, type = 3)

# Try to test the assumptions:
# Homoscedasticity of variances across different levels of the fixed effect:
plot(model_2015_oe_2)
plot(model_2022_oe_2)

# Normality of the residuals:
qqnorm(residuals(model_2015_oe_2))
qqline(residuals(model_2015_oe_2))

qqnorm(residuals(model_2022_oe_2))
qqline(residuals(model_2022_oe_2))

# Normality of the estimated random effects (for each transect):
ranef_brood_2015_oe_2 <- ranef(model_2015_oe_2)$Transect
hist(ranef_brood_2015_oe_2$`(Intercept)`, breaks=10)

ranef_brood_2022_oe_2 <- ranef(model_2022_oe_2)$Transect
hist(ranef_brood_2022_oe_2$`(Intercept)`, breaks=10)

# Because the 2015 model showed a significant effect of treatment, I'll need to
# do pairwise comparisons between treatment groups:
emmeans(model_2015_oe_2, pairwise~Treatment)

# Activity-abundance of forest species models ##################################

# graph the data:
ggplot(dat_2015, aes(x=Treatment, y=forest_specialist_spp_stdz, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)
# At plot 57 (transect E, salvaged) in 2015, there was a high activity-abundance
# of both forest specialists (such as Pterostichus moestus) and eurytopic species 
# (such as Chlaenius emarginatus). I wonder if Pterostichus moestus was also
# abundant in the salvaged plots in 2022

ggplot(dat_2022, aes(x=Treatment, y=forest_specialist_spp_stdz, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_f <- lmerTest::lmer(forest_specialist_spp_stdz ~ Treatment + (1|Transect), 
                                data = dat_2015) # Singular fit
summary(model_2015_f)
# Why might this be a singular fit (which I think means the random effects 
# could not be estimated)?

model_2022_f <- lmerTest::lmer(forest_specialist_spp_stdz ~ Treatment + (1|Transect), 
                                data = dat_2022)
summary(model_2022_f)

anova(model_2015_f, type = 3)
anova(model_2022_f, type = 3)

plot(model_2015_f) # I'm seeing some heteroscedasticity
plot(model_2022_f)

qqnorm(residuals(model_2015_f))
qqline(residuals(model_2015_f))

qqnorm(residuals(model_2022_f))
qqline(residuals(model_2022_f))

ranef_brood_2015_f <- ranef(model_2015_f)$Transect
hist(ranef_brood_2015_f$`(Intercept)`, breaks=10) 

ranef_brood_2022_f <- ranef(model_2022_f)$Transect
hist(ranef_brood_2022_f$`(Intercept)`, breaks=10)

# Take the log of the 2015 data to improve the heteroscedasticity:
dat_2015$log_forest <- log(dat_2015$forest_specialist_spp_stdz)

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=log_forest, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_f_2 <- lmerTest::lmer(log_forest ~ Treatment + (1|Transect),
                             data=dat_2015) # Still a singular fit
summary(model_2015_f_2)

anova(model_2015_f_2, type=3)

plot(model_2015_f_2) # the heteroscedasticity looks better

# Run a basic linear model for the 2015 data:
model_2015_f_3 <- lm(log_forest ~ Treatment, data=dat_2015)
summary(model_2015_f_3)
anova(model_2015_f_3)

#plot(model_2015_f_3)

# Species richness models ######################################################

# graph the data:
ggplot(dat_2015, aes(x=Treatment, y=species_richness, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=species_richness, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Because I have a count response variable (the number of species), I need to
# use a Poisson GLMM
model_2015_rich <- 
  lme4::glmer(species_richness ~ Treatment + (1|Transect),
              data=dat_2015, family="poisson")
summary(model_2015_rich)

model_2022_rich <- 
  lme4::glmer(species_richness ~ Treatment + (1|Transect),
              data=dat_2022, family="poisson") # Singular fit
summary(model_2022_rich)

plot(model_2015_rich)
plot(model_2022_rich)

# Now check for overdispersion in the Poisson GLMM:
# First run a simulation that simulates new values of the random effect intercepts:
model_2015_rich_sim_resid1 <- simulateResiduals(
  fittedModel = model_2015_rich, re.form = NA, plot = T)
# Second, run a simulation that conditions over the previously fitted values
# of the random effect intercepts:
model_2015_rich_sim_resid2 <- simulateResiduals(
  fittedModel = model_2015_rich, re.form = NULL, plot = T)

model_2022_rich_sim_resid1 <- simulateResiduals(
  fittedModel = model_2022_rich, re.form = NA, plot = T)
model_2022_rich_sim_resid1 <- simulateResiduals(
  fittedModel = model_2022_rich, re.form = NA, plot = T)

qqnorm(residuals(model_2015_rich))
qqline(residuals(model_2015_rich))

qqnorm(residuals(model_2022_rich))
qqline(residuals(model_2022_rich))

# Test the null hypothesis of no differences by treatment:
Anova(model_2015_rich, type="III")
Anova(model_2022_rich, type="III")

# Because the 2015 model was significant, I'll run pairwise comparisons:
emmeans(model_2015_rich, pairwise~Treatment)

# Try to run a simplified 2022 GLM without the random effect:
model_2022_rich_1 <- glm(species_richness ~ Treatment,
                         data=dat_2022, family="poisson")
summary(model_2022_rich_1)

# Shannon diversity models #####################################################

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=shannon_diversity, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)
  
ggplot(dat_2022, aes(x=Treatment, y=shannon_diversity, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# I'll now create models using linear mixed-effects models:
model_2015_shannon <- lmerTest::lmer(shannon_diversity ~ Treatment + (1|Transect), 
                                     data = dat_2015)
summary(model_2015_shannon)

model_2022_shannon <- lmerTest::lmer(shannon_diversity ~ Treatment + (1|Transect), 
                                     data = dat_2022)
summary(model_2022_shannon)

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_2015_shannon, type = 3)
anova(model_2022_shannon, type = 3)

# Try to test the assumptions:
# Homoscedasticity of variances across different levels of the fixed effect:
plot(model_2015_shannon)
plot(model_2022_shannon)

# Normality of the residuals:
qqnorm(residuals(model_2015_shannon))
qqline(residuals(model_2015_shannon))

qqnorm(residuals(model_2022_shannon))
qqline(residuals(model_2022_shannon))

# Normality of the estimated random effects (for each transect):
ranef_brood_2015_shannon <- ranef(model_2015_shannon)$Transect
hist(ranef_brood_2015_shannon$`(Intercept)`, breaks=10)

ranef_brood_2022_shannon <- ranef(model_2022_shannon)$Transect
hist(ranef_brood_2022_shannon$`(Intercept)`, breaks=10)

# Functional alpha diversity models ###########################################

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=mean_pairwise_distance, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=mean_pairwise_distance, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# I'll now create models using linear mixed-effects models:
model_2015_mpd <- lmerTest::lmer(mean_pairwise_distance ~ Treatment + (1|Transect), 
                                     data = dat_2015)
summary(model_2015_mpd)

model_2022_mpd <- lmerTest::lmer(mean_pairwise_distance ~ Treatment + (1|Transect), 
                                     data = dat_2022)
summary(model_2022_mpd)

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_2015_mpd, type = 3)
anova(model_2022_mpd, type = 3)

# Try to test the assumptions:
# Homoscedasticity of variances across different levels of the fixed effect:
plot(model_2015_mpd)
plot(model_2022_mpd)

# Normality of the residuals:
qqnorm(residuals(model_2015_mpd))
qqline(residuals(model_2015_mpd))

qqnorm(residuals(model_2022_mpd))
qqline(residuals(model_2022_mpd))

# Normality of the estimated random effects (for each transect):
ranef_brood_2015_mpd <- ranef(model_2015_mpd)$Transect
hist(ranef_brood_2015_mpd$`(Intercept)`, breaks=10)

ranef_brood_2022_mpd <- ranef(model_2022_mpd)$Transect
hist(ranef_brood_2022_mpd$`(Intercept)`, breaks=10)

# Community-weighted mean trait models ########################################

# PC1 #########################################################################

# First I'll do PC1. This variable is a trait syndrome associated with
# longer limbs (longer legs and probably longer antennae), shorter trochanter
# which is associated with less powerful pushing from legs, and a narrower
# pronotum, possibly to reach into small crevices.

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=PC1, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=PC1, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# I'll now create models using linear mixed-effects models:
model_2015_PC1 <- lmerTest::lmer(PC1 ~ Treatment + (1|Transect), 
                                 data = dat_2015)
summary(model_2015_PC1) # It is a singular fit

model_2022_PC1 <- lmerTest::lmer(PC1 ~ Treatment + (1|Transect), 
                                 data = dat_2022)
summary(model_2022_PC1) # It is a singular fit

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_2015_PC1, type = 3)
anova(model_2022_PC1, type = 3)

# Test the assumption of homogeneity of residuals:
plot(model_2015_PC1)
plot(model_2022_PC1)

qqnorm(residuals(model_2015_PC1))
qqline(residuals(model_2015_PC1))

qqnorm(residuals(model_2022_PC1))
qqline(residuals(model_2022_PC1))

ranef_brood_2015_PC1 <- ranef(model_2015_PC1)$Transect
hist(ranef_brood_2015_PC1$`(Intercept)`, breaks=10)

ranef_brood_2022_PC1 <- ranef(model_2022_PC1)$Transect
hist(ranef_brood_2022_PC1$`(Intercept)`, breaks=10)

# Try a basic linear model without the random effects just to make sure that 
# it doesn't impact the results:
model_2015_PC1_1 <- lm(PC1 ~ Treatment, data = dat_2015)
summary(model_2015_PC1_1)
stats::anova(model_2015_PC1_1)

model_2022_PC1_1 <- lm(PC1 ~ Treatment, data = dat_2022)
summary(model_2022_PC1_1)
stats::anova(model_2022_PC1_1)

# Now I need to do pairwise comparisons (use the models that exclude random effect)
emmeans(model_2015_PC1_1, pairwise ~ Treatment)
emmeans(model_2022_PC1_1, pairwise ~ Treatment)

# PC1 followup: pronotum width, rear leg length, rear trochanter length #######

# Because salvaged plot beetles tended to have lower PC1, I'd like to investigate
# specific traits that might contribute to this. Higher PC1 is associated with
# proportionally narrower pronotum, proportionally longer rear legs, and 
# proportionally shorter rear trochanter

# Graph the data (pronotum width standardized to body length):
ggplot(dat_2015, aes(x=Treatment, y=pronotum_width_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2015, aes(x=Treatment, y=log(pronotum_width_standard), 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5) # The log transform doesn't change
# the fact that there is a low "outlier"

ggplot(dat_2022, aes(x=Treatment, y=pronotum_width_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Graph the data (rear leg length standardized to body length):
ggplot(dat_2015, aes(x=Treatment, y=rear_leg_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=rear_leg_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Graph the data (rear trochanter length standardized to body length):
ggplot(dat_2015, aes(x=Treatment, y=rear_trochanter_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=rear_trochanter_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# Run models for standardized pronotum width: ##################################
model_2015_pronotum_width <- lmerTest::lmer(pronotum_width_standard ~ Treatment + 
                                              (1|Transect), data = dat_2015)
summary(model_2015_pronotum_width)

model_2022_pronotum_width <- lmerTest::lmer(pronotum_width_standard ~ Treatment + 
                                              (1|Transect), data = dat_2022)
summary(model_2022_pronotum_width) # singular fit

plot(model_2015_pronotum_width)
qqnorm(residuals(model_2015_pronotum_width))
qqline(residuals(model_2015_pronotum_width)) # One outlier in the pronotum 
# width seems to be breaking the assumption of normally distributed residuals

# Try the 2022 pronotum width model without the random effect:
model_2022_pronotum_width_1 <- lm(pronotum_width_standard ~ Treatment, data = dat_2022)
summary(model_2022_pronotum_width_1)
plot(model_2022_pronotum_width)
qqnorm(residuals(model_2022_pronotum_width))
qqline(residuals(model_2022_pronotum_width)) # the assumptions seem to be met
anova(model_2022_pronotum_width_1)
emmeans(model_2022_pronotum_width_1, pairwise ~ Treatment)

# Run models for standardized rear leg length ##################################

model_2015_leg <- lmerTest::lmer(rear_leg_length_standard ~ Treatment + (1|Transect),
                            data=dat_2015)
summary(model_2015_leg)

model_2022_leg <- lmerTest::lmer(rear_leg_length_standard ~ Treatment + (1|Transect),
                                 data=dat_2022)
summary(model_2022_leg)

# test assumptions: 
plot(model_2015_leg) # heteroscedasticity doesn't look too bad
plot(model_2022_leg) # looks homoscedastic

qqnorm(residuals(model_2015_leg))
qqline(residuals(model_2015_leg)) # The residuals do NOT appear to be super 
# normally distributed

qqnorm(residuals(model_2022_leg))
qqline(residuals(model_2022_leg)) # residuals look normally distributed

# Run the anova tests:
anova(model_2015_leg, type = 3) # I think this is running the ANOVA using 
# Satterthwaite's method of approximating degrees of freedom
anova(model_2022_leg, type = 3)

# Run models for standardized rear trochanter length ###########################

model_2015_trochanter <- lmerTest::lmer(rear_trochanter_length_standard ~ Treatment +
                                          (1|Transect), data=dat_2015)
summary(model_2015_trochanter)

model_2022_trochanter <- lmerTest::lmer(rear_trochanter_length_standard ~ Treatment +
                                          (1|Transect), data=dat_2022)
summary(model_2022_trochanter)

# Test assumptions:
plot(model_2015_trochanter)
plot(model_2022_trochanter)

qqnorm(residuals(model_2015_trochanter))
qqline(residuals(model_2015_trochanter)) 

qqnorm(residuals(model_2022_trochanter))
qqline(residuals(model_2022_trochanter)) 

# Run the anova tests:
anova(model_2015_trochanter, type = 3)
anova(model_2022_trochanter, type = 3)

# Pairwise test:
emmeans(model_2015_trochanter, pairwise ~ Treatment)

# PC2 #########################################################################

# Graph the data:
ggplot(dat_2015, aes(x=Treatment, y=PC2, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=PC2, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

# I'll now create models using linear mixed-effects models:
model_2015_PC2 <- lmerTest::lmer(PC2 ~ Treatment + (1|Transect), 
                                 data = dat_2015)
summary(model_2015_PC2)

model_2022_PC2 <- lmerTest::lmer(PC2 ~ Treatment + (1|Transect), 
                                 data = dat_2022)
summary(model_2022_PC2)

# Run the anova tests to test the null hypothesis of no treatment differences:
anova(model_2015_PC2, type = 3)
anova(model_2022_PC2, type = 3)

# Test the assumption of homogeneity of residuals:
plot(model_2015_PC2)
plot(model_2022_PC2)

qqnorm(residuals(model_2015_PC2))
qqline(residuals(model_2015_PC2))

qqnorm(residuals(model_2022_PC2))
qqline(residuals(model_2022_PC2))

ranef_brood_2015_PC2 <- ranef(model_2015_PC2)$Transect
hist(ranef_brood_2015_PC2$`(Intercept)`, breaks=10)

ranef_brood_2022_PC2 <- ranef(model_2022_PC2)$Transect
hist(ranef_brood_2022_PC2$`(Intercept)`, breaks=10)

# Now I need to do pairwise comparisons (just for 2015 because it was significant)
emmeans(model_2015_PC2, pairwise ~ Treatment)

# Because only the 2015 model was significant for PC2, I'll only need to 
# investigate the individual traits for 2015.

# PC2 followup: body length ####################################################

ggplot(dat_2015, aes(x=Treatment, y=body_length, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_body <- lmerTest::lmer(body_length ~ Treatment + (1|Transect), 
                                  data = dat_2015) # singular fit, need to
# try a regular linear model:
model_2015_body_1 <- lm(body_length ~ Treatment, data = dat_2015)
summary(model_2015_body_1)

# Test the assumptions of the regular linear model:
#plot(model_2015_body_1)

stats::anova(model_2015_body_1)
car::Anova(model_2015_body_1)

# PC2 followup: standardized eye length ########################################
ggplot(dat_2015, aes(x=Treatment, y=eye_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_eye_length <- lmerTest::lmer(eye_length_standard ~ Treatment + (1|Transect), 
                                        data = dat_2015)
summary(model_2015_eye_length)

# test assumptions:
plot(model_2015_eye_length)

qqnorm(residuals(model_2015_eye_length))
qqline(residuals(model_2015_eye_length))

# Run the anova test:
anova(model_2015_eye_length, type=3)

# Pairwise comparisons:
emmeans(model_2015_eye_length, pairwise ~ Treatment)

# PC2 followup: antenna to rear leg ratio ######################################

ggplot(dat_2015, aes(x=Treatment, y=antenna_rear_leg_ratio,
                     color=Transect, group=Transect)) +
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_antenna_rear_leg <- lmerTest::lmer(antenna_rear_leg_ratio ~ Treatment + (1|Transect), 
                                              data = dat_2015)
anova(model_2015_antenna_rear_leg, type=3)

# PC3 ##########################################################################

ggplot(dat_2015, aes(x=Treatment, y=PC3, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=PC3, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_PC3 <- lmerTest::lmer(PC3 ~ Treatment + (1|Transect), 
                                                   data = dat_2015)
summary(model_2015_PC3)

model_2022_PC3 <- lmerTest::lmer(PC3 ~ Treatment + (1|Transect), 
                                 data = dat_2022)
summary(model_2022_PC3)


anova(model_2015_PC3, type=3)
anova(model_2022_PC3, type=3)


# Water affinity ###############################################################

# Plot the data:
ggplot(dat_2015, aes(x=Treatment, y=Water_affinity, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=Water_affinity, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_water <- lmerTest::lmer(Water_affinity ~ Treatment + (1|Transect), 
                                                     data = dat_2015) # singular fit
summary(model_2015_water)
model_2015_water_1 <- lm(Water_affinity ~ Treatment, 
                           data = dat_2015)
summary(model_2015_water_1)
Anova(model_2015_water_1, type="III")
anova(model_2015_water_1)

model_2022_water <- lmerTest::lmer(Water_affinity ~ Treatment + (1|Transect), 
                                   data = dat_2022) # singular fit
model_2022_water_1 <- lm(Water_affinity ~ Treatment, 
                         data = dat_2022)
summary(model_2022_water_1)
anova(model_2022_water_1)


# Flight capability ############################################################

ggplot(dat_2015, aes(x=Treatment, y=Flight_capability, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=Flight_capability, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_flight <- lmerTest::lmer(Flight_capability ~ Treatment + (1|Transect),
                                    data=dat_2015)
summary(model_2015_flight)

model_2022_flight <- lmerTest::lmer(Flight_capability ~ Treatment + (1|Transect),
                                    data=dat_2022)
summary(model_2022_flight)

# Check assumptions:
plot(model_2015_flight)
plot(model_2022_flight)

qqnorm(residuals(model_2015_flight))
qqline(residuals(model_2015_flight))

qqnorm(residuals(model_2022_flight))
qqline(residuals(model_2022_flight)) # There is a long tail on the right
hist(residuals(model_2022_flight))

# Run the hypothesis tests:
anova(model_2015_flight, type=3)
anova(model_2022_flight, type=3)

# Pairwise tests:
emmeans(model_2015_flight, pairwise ~ Treatment)

# Standardized antenna length ##################################################

# graph:
ggplot(dat_2015, aes(x=Treatment, y=antenna_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

ggplot(dat_2022, aes(x=Treatment, y=antenna_length_standard, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_antenna <- lmerTest::lmer(antenna_length_standard ~ Treatment + (1|Transect),
                                     data=dat_2015)
summary(model_2015_antenna)
anova(model_2015_antenna, type=3)

model_2022_antenna <- lmerTest::lmer(antenna_length_standard ~ Treatment + (1|Transect),
                                     data=dat_2022) # singular fit
summary(model_2022_antenna)

# run model without the random effect:
model_2022_antenna_1 <- lm(antenna_length_standard ~ Treatment, data=dat_2022)
summary(model_2022_antenna_1)
anova(model_2022_antenna_1)


# Pretty graphs ###############################################################

abundance_graph_2015 <- ggplot(dat_2015, aes(x=Treatment, y=total_count_stdz)) + 
  stat_summary(fun = mean, geom = "bar", fill = "lightgrey", color="black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_quasirandom(width=0.05, alpha=0.5) +
  ylab("Ground beetles caught per day") + theme(plot.title = element_text(size=18),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_text(size = 16, 
                                                                            margin = margin(r=20)),
                                                axis.text.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 14))+
  ylim(c(0,2.5))+
  ggtitle("Activity-abundance in 2015")
abundance_graph_2015

abundance_graph_2022 <- ggplot(dat_2022, aes(x=Treatment, y=total_count_stdz)) + 
  stat_summary(fun = mean, geom = "bar", fill = "lightgrey", color="black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_quasirandom(width=0.05, alpha=0.5) +
  ylab("Ground beetles caught per day") + theme(plot.title = element_text(size=18),
                                                axis.title.x = element_blank(),
                                                axis.title.y = element_text(size = 16, 
                                                                            margin = margin(r=20)),
                                                axis.text.x = element_text(size = 14),
                                                axis.text.y = element_text(size = 14))+
  ylim(c(0,2.5))+
  ggtitle("Activity-abundance in 2022")
abundance_graph_2022

richness_graph_2015 <- ggplot(dat_2015, aes(x=Treatment, y=species_richness)) + 
  stat_summary(fun = mean, geom = "bar", fill = "lightgrey", color="black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_quasirandom(width=0.1, alpha=0.5) +
  ylab("Number of species") + theme(plot.title = element_text(size=18),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(size = 16, 
                                                                margin = margin(r=15)),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14))+
  ylim(0,20)+
  ggtitle("Species richness in 2015")
richness_graph_2015

richness_graph_2022 <- ggplot(dat_2022, aes(x=Treatment, y=species_richness)) + 
  stat_summary(fun = mean, geom = "bar", fill = "lightgrey", color="black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_quasirandom(width=0.1, alpha=0.5) +
  ylab("Number of species") + theme(plot.title = element_text(size=18),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(size = 16, 
                                                                margin = margin(r=15)),
                                    axis.text.x = element_text(size = 14),
                                    axis.text.y = element_text(size = 14))+
  ylim(0,20)+
  ggtitle("Species richness in 2022")
richness_graph_2022

ggarrange(abundance_graph_2015, abundance_graph_2022,
          richness_graph_2015, richness_graph_2022,
          labels = c("A","B","C","D"), ncol=2, nrow=2)


# Make summary data tables for treatment means:

response_vars <- c("total_count_stdz", "open_habitat_spp_stdz", 
                   "eurytopic_spp_stdz", "forest_specialist_spp_stdz",
                   "species_richness", "shannon_diversity", "mean_pairwise_distance",
                   "PC1", "PC2", "PC3", "pronotum_width_standard", "rear_leg_length_standard",
                   "rear_trochanter_length_standard", "eye_length_standard",
                   "body_length", "antenna_rear_leg_ratio", "Water_affinity",
                   "Flight_capability", "antenna_length_standard")

dat_2015_by_treatment <- dat_2015 %>% group_by(Treatment) %>%
  summarize(across(all_of(response_vars), mean))

#write.csv(dat_2015_by_treatment, "Aaron_PNR_formatted_data/PNR2015_response_vars_by_treatment.csv", row.names = F)

dat_2022_by_treatment <- dat_2022 %>% group_by(Treatment) %>%
  summarize(across(all_of(response_vars), mean))

#write.csv(dat_2022_by_treatment, "Aaron_PNR_formatted_data/PNR2022_response_vars_by_treatment.csv", row.names = F)

























