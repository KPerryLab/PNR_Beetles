dat$Treatment <- factor(dat$Treatment, levels = c("Windthrow", "Salvaged", "Forest"))
treatment_colors = c("Forest" = "palegreen3", "Salvaged" = "goldenrod2", "Windthrow" = "brown4")

abundance_graph <- ggplot(dat, aes(x = Treatment, y = total_count_stdz, fill = Year, group = Year)) + 
  stat_summary(fun = mean, geom = "bar", color="black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9)) +
  ylab("Number of individuals / plot") + theme(plot.title = element_text(size=18),
                                               axis.title.x = element_blank(),
                                               axis.title.y = element_text(size = 16, 
                                                                           margin = margin(r=20)),
                                               axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                                               axis.text.y = element_text(size = 14),
                                               legend.position = "none")+
  scale_fill_grey() + coord_cartesian(ylim = c(0,110)) # Using coord_cartesian 
# instead of scale_y_continuous allows the graph to be "zoomed in" without omitting
# any data points in the summary statistics
abundance_graph

abundance_oe_graph <- ggplot(dat, aes(x = Treatment, 
                                      y = eurytopic_spp_stdz + open_habitat_spp_stdz, fill = Year, group = Year)) + 
  stat_summary(fun = mean, geom = "bar", color="black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9)) +
  ylab("Open-habitat or habitat-generalist \nindividuals / plot") + theme(plot.title = element_text(size=18),
                                                                          axis.title.x = element_blank(),
                                                                          axis.title.y = element_text(size = 16, 
                                                                                                      margin = margin(r=20)),
                                                                          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                                                                          axis.text.y = element_text(size = 14),
                                                                          legend.position = "none")+
  scale_fill_grey() + coord_cartesian(ylim = c(0,110))
abundance_oe_graph

abundance_forest_graph <- ggplot(dat, aes(x = Treatment, 
                                          y = forest_specialist_spp_stdz, fill = Year, group = Year)) + 
  stat_summary(fun = mean, geom = "bar", color="black", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.9)) +
  ylab("Forest-specialist \nindividuals / plot") + theme(plot.title = element_text(size=18),
                                                         axis.title.x = element_blank(),
                                                         axis.title.y = element_text(size = 16, 
                                                                                     margin = margin(r=20)),
                                                         axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
                                                         axis.text.y = element_text(size = 14),
                                                         legend.title = element_text(size = 14),
                                                         legend.text = element_text(size = 14))+
  scale_fill_grey() + coord_cartesian(ylim = c(0,110))
abundance_forest_graph
#angle = 45, hjust = 1 
empty_graph <- ggplot() + theme_void()

ggarrange(abundance_graph, empty_graph, abundance_oe_graph, empty_graph, abundance_forest_graph,
          labels = c("A", "", "B", "", "C"), ncol=5, nrow=1, widths = c(0.65, 0.1, 0.70, 0.1, 1))

# legend.position = "none"

# Richness graph ###############################################################

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