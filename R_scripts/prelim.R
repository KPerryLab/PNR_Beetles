## Preliminary ground beetle analysis
## Powdermill Nature Reserve 2022
## 13 March 2024
# Written by Kayla Perry, edited or modified by Aaron Tayal

library(readxl)
library(hillR)
library(ggplot2)

# load the data
dat <- read_excel("PNR_Raw_Data/PNR2022_InvertebrateCommunity.xlsx", sheet = 2)

# let's focus on data from July
# pull out the intervals we want and create a new data frame

str(dat)
dat$Collection_interval <- as.factor(dat$Collection_interval)
levels(dat$Collection_interval)

july1 <- dat[which(dat$Collection_interval == "2"),]
july2 <- dat[which(dat$Collection_interval == "3"),]

dat2 <- rbind(july1, july2)

# now we need to replace the NAs with zeros (for the abundances)
dat2[,13:49][is.na(dat2[,13:49])] <- 0

# looks like we may have some missing samples since they are not represented with a zero
# let's remove those
dat3 <- dat2[!is.na(dat2$Total_Carabidae),]

# now let's calculate richness
library(hillR)
dat3$Treatment <- as.factor(dat3$Treatment)
levels(dat3$Treatment)

dat3$rich <- hill_taxa(dat3[,13:49], q = 0, MARGIN = 1)

dotchart(dat3$rich, group = dat3$Treatment, pch = 19)
hist(dat3$rich)
boxplot(dat3$rich ~ Treatment, data = dat3)
stripchart(dat3$rich ~ Treatment, data = dat3, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

#Create a species accumulation curve:
forest_july <- dat3[dat3$Treatment=="F",]
wind_july <- dat3[dat3$Treatment=="W",]
salvage_july <- dat3[dat3$Treatment=="S",]
sp.forest_july <- specaccum(forest_july[13:50], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.wind_july <- specaccum(wind_july[13:50], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.salvage_july <- specaccum(salvage_july[13:50], method = "rarefaction", permutations = 100, gamma = "jack2")

plot(sp.forest_july, col = "#481567FF", xvar = c("individuals"), lty = 4, lwd = 4,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 150), ylim = c(0, 25))
plot(sp.wind_july, add = TRUE, xvar = c("individuals"), lty = 1, lwd = 4, col = "#FDE725FF")
plot(sp.salvage_july, add = TRUE, xvar = c("individuals"), lty = 2, lwd = 4, col = "#73D055FF")

legend("bottomright", legend = c("Forest", "Tornado", "Salvaged"),
       pch = c(16, 17, 15, 18), lty = c(1,2,3,4), cex = 1, bty = "n", lwd = 4,
       col = c("#481567FF", "#FDE725FF", "#73D055FF"))


# Shannon diversity
dat3$div <- hill_taxa(dat3[,13:49], q = 1, MARGIN = 1)

dotchart(dat3$div, group = dat3$Treatment, pch = 19)
hist(dat3$div)
boxplot(dat3$div ~ Treatment, data = dat3)
stripchart(dat3$div ~ Treatment, data = dat3, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)


# now let's assess species composition
# let's try presence/absence
library(vegan)

colSums(dat3[,13:49])

dat4 <- dat3[,13:49]
dat4[dat4 > 0] <- 1
dat4$Treatment <- dat3$Treatment
str(dat4)
summary(dat4)

#remove columns (aka species) which weren't found at all in july
dat5 <- dat4[, colSums(dat4 !=0) > 0]
colSums(dat5[,1:23])

dis.matrix <- vegdist(dat5[,1:23], method = "jaccard")

nmds <- metaMDS(dis.matrix, trymax = 500, autotransform = TRUE, k = 2)
nmds
stressplot(nmds)
plot(nmds)

ordiplot(nmds, disp = "sites", type = "n", xlim = c(-5, 2), ylim = c(-1.5, 1.5))
points(nmds, dis = "sites", select = which(dat5$Treatment=="F"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds, dis = "sites", select = which(dat5$Treatment=="W"), pch = 18, cex = 2, col = "#481567FF")
points(nmds, dis = "sites", select = which(dat5$Treatment=="S"), pch = 15, cex = 2, col = "#2D708EFF")
orditorp(nmds, "sites")


