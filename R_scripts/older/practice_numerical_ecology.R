# Numerical Ecology in R, by Borcard, Gillet, and 
# Legendre

# Aaron Tayal Feb 5 2025

library(ade4)

data(doubs, package="ade4")

spe <- doubs$fish
View(spe)

hist(unlist(spe))

spa <- doubs$xy

plot(spa)
lines(spa)
text(spa, row.names(spa), cex=0.8, col="red")

# Compute classical diversity indices (page 17):
library(vegan)
?diversity
N0 <- rowSums(spe > 0) # species richness
H <- diversity(spe, index="shannon") # shannon entropy
N1 <- exp(H) # shannon diversity number

# Plot a bubble diagram of shannon diversity:
plot(spa, main="Shannon diversity number exp(H)", pch=21, cex=2*N1/max(N1)) + 
  lines(spa)

J <- H/log(N0)
plot(spa[2:30,], main="Pielou evenness", pch=21, cex=J[2:30]/max(J[2:30])) + 
  lines(spa) # The Pielou evenness seems not to vary that much

E1 <- N1/N0 # Shannon evenness
plot(spa[-8], main="Shannon evenness (Hill ratio)", pch=21, cex=E1[-8]/max(E1[-8])) + 
  lines(spa)
# Because the abundances are on a 1-5 scale of commonness, rather
# than being actual numerical abundance, then it is not 
# easy to interpret values such as evenness

# Question: what is the minimum possible value that Shannon diversity exp(H)
# can reach?
test <- c(0.9999, 0.0001)
sum(test)
diversity(test)





