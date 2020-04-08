## TB benthic core analysis

# working script for updating metadata and sample counts for each sampling date.

TB_2018_core_meta = TB_2018_core_meta %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
TB_2019_core_meta = TB_2019_core_meta %>% rename(Date = "Date Taken") %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

TB_2018_core_summ = TB_2018_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())
TB_2019_core_summ = TB_2019_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())

TB_core_summ = bind_rows(TB_2018_core_summ, TB_2019_core_summ)

#get number of cores processed
shallow_summ <- shallow %>% rownames_to_column("site") %>% 
  mutate(site = str_remove(site,"-.+")) %>%
  group_by(site) %>%
  summarise(N = n())

