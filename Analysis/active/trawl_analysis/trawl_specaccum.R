
# Full species accumulation curve (SAC)

##create taxa by site matrix
TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

sp1 <- vegan::specaccum(TB_trawl_commonsite %>% select(-year, -month, -date_id))
spp_mod.df <- data.frame(trawls = sp1[["sites"]], richness = sp1[["richness"]], richness_sd = sp1[["sd"]])
# plot(sp1, ci.type = "poly", col = "nlue", lwd = 2, ci.lty = 0)
# sp2 <- vegan::specaccum(TB_trawl_taxasite %>% select(-year, -month, -date_id), 'random')
# boxplot(sp2, col = "yellow", add = TRUE, pch = "+")

source("./figures/SAC_full_plot.R");SAC_full_plot()

## SAC for each year
TB_trawl_commondate <- TB_trawl_data %>%
  select(year, Date, Common_name) %>%
  group_by(Date, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()


TB_trawl_taxasite_yr <- TB_trawl_commondate %>%
  group_split(year) %>%
  map(., ~.x %>%  select(-year, -Date) %>% vegan::specaccum(.) %>%
       rlist::list.remove(names(.) %ni% c("sites","richness","sd")) %>%
        bind_cols) %>% 
  map2(., TB_trawl_commondate %>% group_split(year), ~..1 %>% 
         mutate(year = ..2 %>% select(year) %>% unlist %>% as.factor)) %>%
  bind_rows

source("./figures/SAC_yr_plot.R");SAC_yr_plot()

#create time-accumulation df from trawl data
TB_trawl_data %>%
  filter(!duplicated(Common_name)) %>%
  mutate(richness = 1:n(), 
         Date = as.Date(Date),
         richness = sort(richness))-> TB_time_accum

TB_annual_accum <- TB_trawl_data %>%
  mutate(pres = 1) %>% group_by(year) %>%
  filter(!duplicated(Common_name)) %>%
  summarise(richness = sum(pres)) %>%
  mutate(Date = as.Date(as.character(year), format = "%Y", origin = "2000-01-01"))

# partition beta diversity into nestedness and turnover #
# uses `betapart` package #
TB_time_df <- TB_time_accum %>%
  select(Date, year, richness) %>%
  mutate(data_type = 'continuous') %>%
  bind_rows(TB_annual_accum %>% mutate(data_type = "annual"))

source("./figures/diversity_time_plots.R")


#create moving jaccard calc by year-month

date_lists = TB_trawl_commonsite %>%
  group_split(date_id)

unique_comb


ggplot(TB_time_df, aes(x = Date, y = richness, group = data_type)) +
  geom_point(aes(colour = data_type),size = 3) + 
  geom_line(aes(colour = data_type),size = 2) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  scale_y_continuous(name = "Species Richness")+
  scale_colour_manual(name = "", values = viridis(4)[c(2,3)]) +
  theme(legend.position = c(0.1,0.8), axis.title.x = element_blank())


