
library(ggplot2)

setwd("C:/Users/Aaron/OneDrive/Documents/Entomology MS/Thesis research/PNR_Beetles/Aaron ideas and analyses")
dat <- read.csv("dummy_for_graph.csv")
ggplot(data=dat, aes(x=Time, y=dummy))+
  geom_point(color="white")+
  theme_classic()+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  scale_x_continuous(name="", breaks=c(2012.5, 2013.5, 2015.5, 2022.5), 
                     labels=c("June 2012", "June 2013", "June 2015", "June 2022"), 
                     limits=c(2011, 2024))
