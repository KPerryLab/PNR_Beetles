# Aaron Tayal
# 4/30/2025
# Carabids stats

# The purpose is to create linear mixed-effects models to investigate
# metrics of ground beetle activity-abundance, species richness, Shannon
# diversity, and functional alpha diversity, and whether these variables
# differ between plots that were salvaged, windthow-affected, or forest control.
# And whether they differ between 2015 and 2022.

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

# I'll now create models using linear mixed-effects models:
model_antenna_rear_leg_ratio <- lmerTest::lmer(antenna_rear_leg_ratio ~ Treatment + Year + Treatment*Year +
                                                   (1|Transect), data = dat)
summary(model_antenna_rear_leg_ratio)

# plot the model:
plot_model(model_antenna_rear_leg_ratio, type = "pred", terms = c("Treatment", "Year"))

# test assumptions:
plot(model_antenna_rear_leg_ratio)
qqnorm(residuals(model_antenna_rear_leg_ratio))
qqline(residuals(model_antenna_rear_leg_ratio))

# run the ANOVA test:
anova(model_antenna_rear_leg_ratio, type=3)

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

# pairwise comparisons:
emmeans(model_eye_protrusion_ratio, pairwise~Treatment)




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

























