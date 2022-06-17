### diversity by time plot ###

 rich_time <- 
  ggplot(TB_time_df, aes(x = Date, y = richness, group = data_type)) +
  geom_point(aes(colour = data_type),size = 1.5) +
  geom_line(aes(colour = data_type),size = 1.2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  scale_y_continuous(name = "Species Richness")+
  scale_colour_manual(name = "", values = viridis(4)[c(2,3)]) +
  theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(fill = "transparent", colour = NA),
        # axis.title.x = element_blank(), axis.text.x = element_blank(),
        legend.spacing.y = unit(-0.1, 'cm'));rich_time

 jaccard_time <-
  ggplot(beta_df, aes(x = Date, y = value, group = partition)) + 
   geom_point(aes(colour = partition), size = 1.5) +
  geom_line(aes(colour = partition), size = 1.2) +
  scale_y_continuous(name = "Jaccard dissimilarity", limits = c(0,1.1), expand = c(0.01,0.01))+
  scale_colour_manual(name = "", values = viridis(7)[c(2,4,7)], labels = c("jaccard", "nestedness","turnover")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  theme(legend.background = element_rect(colour = NA, fill = 'transparent'),
        legend.position = c(0,1), legend.justification = c(0,1), axis.title.x = element_blank(),
        legend.spacing.y = unit(-0.2, 'cm'))

 tiff("./figures/TB_temporal_diversity.tiff", res = 400, width = 5, height = 7.5, units = "in", compression = "lzw")
 grid.arrange(rich_time, jaccard_time, ncol = 1)
 dev.off()
 