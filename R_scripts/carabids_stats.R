# Aaron Tayal
# 4/30/2025
# Carabids stats

# The purpose is to create linear mixed-effects models to investigate
# metrics of ground beetle activity-abundance, species richness, Shannon
# diversity, and functional alpha diversity

library(ggplot2)
theme_set(theme_classic())
library(lme4) # linear mixed-effects models
library(lmerTest) # uses Satterthwaites approx. for degrees of freedom

dat_2015 <- read.csv("Aaron_PNR_formatted_data/PNR2015_carabid_counts_by_plot_standardized.csv")

dat_2022 <- read.csv("Aaron_PNR_formatted_data/PNR2022_carabid_counts_by_plot_standardized.csv")

dat_2015$Treatment <- as.factor(dat_2015$Treatment)
dat_2022$Treatment <- as.factor(dat_2022$Treatment)
dat_2015$Transect <- as.factor(dat_2015$Transect)
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
plot(model_2022_total)

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

# A log transformation of 2015 data would likely not improve the model 
# meeting assumptions very much, because the range of the standardized total 
# count is between 0.04 and 2.2, and the log function is *fairly* linear
# over that domain

# Alternatively, I could try a Kruskal-Wallis test:
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
# (such as Chlaenius emarginatus)

ggplot(dat_2022, aes(x=Treatment, y=forest_specialist_spp_stdz, 
                     color=Transect, group=Transect)) + 
  geom_point(alpha=0.5) + geom_line(alpha=0.5)

model_2015_f <- lmerTest::lmer(forest_specialist_spp_stdz ~ Treatment + (1|Transect), 
                                data = dat_2015)
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












