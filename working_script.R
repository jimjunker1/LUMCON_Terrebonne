# working script for updating metadata and sample counts for each sampling date.

TB_2018_core_meta = TB_2018_core_meta %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
TB_2019_core_meta = TB_2019_core_meta %>% rename(Date = "Date Taken") %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

TB_2018_core_summ = TB_2018_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())
TB_2019_core_summ = TB_2019_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())

TB_core_summ = bind_rows(TB_2018_core_summ, TB_2019_core_summ)

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
                       