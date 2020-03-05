# working figures #
TB_core_day = TB_core_summ %>%
  mutate(yday = yday(Date)) %>%
  group_by(yday) %>%
  mutate(day = as.Date(yday-1,  origin = "2020-01-01"))
hydro_plot = ggplot(hydro_summ, aes(x = DATE)) + 
  geom_line(aes(y = temp_C_na.rm_mean), colour = 'blue', size = 1.5)+
  geom_line(aes(y = sal_psu_na.rm_mean), colour = 'lightblue', size = 1.5)+
  geom_ribbon(aes(ymin = temp_C_min, ymax = temp_C_max), fill = 'blue', alpha = 0.5)+
  geom_ribbon(aes(ymin= sal_psu_min, ymax = sal_psu_max), fill = 'lightblue', alpha = 0.5) +
  scale_y_continuous(name = expr("Temperature ("*degree*"C) | Salinity (PSU)"))+
  scale_x_date(breaks = "1 months", date_labels = "%b", limits = c(as.Date("2020-01-01"),as.Date("2020-12-31")), expand = c(0,0.5)) +
  geom_vline(data = TB_core_day, aes(xintercept = day), size = 2) +
  theme(panel.grid.minor.x = element_line(colour = 'lightgrey'),
        panel.grid.major.x = element_line(colour = 'lightgrey'));hydro_plot

do_plot = ggplot(do_summ, aes(x = DATE)) +
  geom_line(aes(y = na.rm_mean), colour = 'red', size = 1.5) +
  geom_ribbon(aes(ymin = min, ymax = max), fill = 'red', alpha = 0.5)+
  scale_y_continuous(name = 'Dissolved Oxygen (mg/L)', limits = c(0,NA))+
  scale_x_date(breaks = "1 months", date_labels = "%b", limits = c(as.Date("2020-01-01"),as.Date("2020-12-31")), expand = c(0,0.5)) +
  geom_vline(data = TB_core_day, aes(xintercept = day), size = 2)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid.minor.x = element_line(colour = 'lightgrey'),
        panel.grid.major.x = element_line(colour = 'lightgrey'));do_plot

tiff(file = "./figures/env_plot.tiff", res = 400, width = 5, height = 4, units = "in", compression = "lzw")
grid.draw(gtable_rbind(ggplotGrob(do_plot), ggplotGrob(hydro_plot)))
dev.off()
## core counts by site 
core_site_plot = ggplot(TB_core_summ, aes(x = Station, y = Date)) + 
  geom_point(aes(fill = core_count,size = core_count), shape = 21) +
  # geom_point(aes(size = ))
  scale_size_continuous(name = "# cores") +
  scale_fill_viridis(name = "# cores")+
  coord_flip() + guides(size = FALSE)+
  scale_y_date(name = "Date",date_breaks = "3 months", 
               limits = c(as.Date(min(TB_core_summ$Date))-30,as.Date(max(TB_core_summ$Date))+30),
               date_labels = "%b-%y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

tiff(file = "./figures/cores_by_site.tiff", res = 400, width = 5, height = 4, units = "in", compression = "lzw")
core_site_plot
dev.off()


####            
