### diversity by time plot ###

rich_time <- ggplot(TB_time_df, aes(x = Date, y = richness, group = data_type)) +
  geom_point(aes(colour = data_type),size = 3) + 
  geom_line(aes(colour = data_type),size = 2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  scale_y_continuous(name = "Species Richness")+
  scale_colour_manual(name = "", values = viridis(4)[c(2,3)]) +
  theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.background = element_rect(fill = "transparent", colour = NA),
        axis.title.x = element_blank())


