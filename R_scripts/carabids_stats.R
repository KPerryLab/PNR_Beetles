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

model_2022_f <- lmerTest::lmer(forest_specialist_spp_stdz ~ Treatment + (1|Transect), 
                                data = dat_2022)
summary(model_2022_f)

anova(model_2015_f, type = 3)
anova(model_2022_f, type = 3)

plot(model_2015_f) # I'm seeing some heteroscedasticity
plot(model_2022_f) # Also seeing some heteroscedasticity

qqnorm(residuals(model_2015_f))
qqline(residuals(model_2015_f))

qqnorm(residuals(model_2022_f))
qqline(residuals(model_2022_f))

ranef_brood_2015_f <- ranef(model_2015_f)$Transect
hist(ranef_brood_2015_f$`(Intercept)`, breaks=10) 

ranef_brood_2022_f <- ranef(model_2022_f)$Transect
hist(ranef_brood_2022_f$`(Intercept)`, breaks=10)

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
              data=dat_2022, family="poisson")
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

# Shannon diversity models

# Functional alpha diversity models

# Community-weighted mean trait models










