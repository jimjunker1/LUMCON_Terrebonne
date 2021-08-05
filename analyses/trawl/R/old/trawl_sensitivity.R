source("datascript.R")

## trawl sampling sensitivity analsis ##
# Questions are:
# 1) how many trawls should be done over the year?
# 2) how should those trawls be distributed?
# 3) evenly? clumped? etc.

# import the data
TB_trawl_taxasite <- TB_trawl_data %>% 
dplyr::select(Date, Common_name, Abundance) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup()

TB_trawl_richness_yr <- TB_trawl_taxasite %>%
  mutate(year = year(Date),
         month = month(Date)) %>% 
  filter(year != 2020) %>%
  filter(month %ni% c(4,5)) %>%
  dplyr::select(Date, year, month, everything()) %>%
  junkR::named_group_split(year)

TB_trawl_richness <- TB_trawl_data %>% na.omit %>% mutate(year = year(Date)) %>%
  filter(year != 2020) %>% 
  dplyr::select(year, Common_name) %>% distinct %>%
  group_by(year) %>% count %>% bind_rows(data.frame(year= 2020, n = 75))

richness_2020 <- TB_trawl_taxasite %>% na.omit %>% mutate(year = year(Date)) %>%
  filter(year == 2020) %>% dplyr::select(-year) %>%
  pivot_longer(cols = -Date, names_to = "Common_name", values_to = "Abundance") %>%
  filter(Abundance > 0 & !is.na(Abundance)) %>%
  dplyr::select(Date, Common_name) %>% distinct %>%
  group_by(Date) %>% count %>% 
  rename(trawl_richness = n) %>% ungroup %>%
  mutate(trawls = 1:n()) %>% dplyr::select(-Date) %>%
  mutate(rel_richness = trawl_richness/75*100)


sample_n <- as.list(1:36)

set.seed(42)
trawl_sample_simulation <- lapply(sample_n, function(x) lapply(TB_trawl_richness_yr, dplyr::sample_n, size = x, replace = TRUE))

for(i in seq_along(trawl_sample_simulation)){
  trawl_sample_simulation[[i]] <- trawl_sample_simulation[[i]] %>%
    bind_rows %>% mutate(trawls = i) %>% dplyr::select(Date, year, month, trawls, everything())
} 

trawl_sample_simulation <- trawl_sample_simulation %>% bind_rows %>% 
  dplyr::select(Date, year, month, trawls, everything()) %>%
    pivot_longer(cols = -Date:-trawls, names_to = "Common_name", values_to = "Abundance") %>%
  filter(Abundance > 0 & !is.na(Abundance)) %>% dplyr::select(-Date, -month, -Abundance) %>%  group_by(year, trawls) %>%
  distinct %>% count %>% rename(trawl_richness = n) %>% left_join(TB_trawl_richness) %>%
  mutate(rel_richness = trawl_richness/n*100)

trawl_estimates <- trawl_sample_simulation %>%
  filter(trawls %in% c(9, 15)) %>% group_by(trawls) %>% 
  summarise(mean_rel = mean(rel_richness),
            rel_90 = quantile(rel_richness, 0.9),
            rel_10 = quantile(rel_richness, 0.1))

ggplot(trawl_sample_simulation, aes(x = trawls, y = rel_richness, group = factor(year))) + 
  geom_path(aes(color = factor(year)), size = 1.2) +
  scale_y_continuous(name = "Relative Spp. Richness (%)") +
  scale_x_continuous(name = "# of trawls")+
  scale_color_viridis(discrete = TRUE) +
  geom_vline(xintercept = c(3,9,15), linetype = c("dashed", rep("solid",2))) +
  annotate('text', x = 2.9, y = 0, label = "Current #", hjust = 1) +
  annotate("text", x = 8.8, y = 0, label = "1 trawl/month", hjust = 1)+
  annotate("text", x = 15.1, y =0, label = "2 trawls/month", hjust = 0) +
  annotate('text', x = 20, y = 35, label = "Total species 'recovery'", hjust = 0, fontface = "bold")+
  annotate('text', x = 20, y = 30, label = paste0("1 trawl/m = ",round(trawl_estimates[1,"mean_rel"],1)," [ ",
                                                   round(trawl_estimates[1,"rel_10"],1)," to ",
                                                   round(trawl_estimates[1,"rel_90"],1)," ]"), hjust = 0) +
  annotate('text', x = 20, y = 25, label = paste0("2 trawl/m = ",round(trawl_estimates[2,"mean_rel"],1)," [ ",
                                                    round(trawl_estimates[2,"rel_10"],1)," to ",
                                                    round(trawl_estimates[2,"rel_90"],1)," ]"), hjust = 0) +
  ggtitle("2020 outlook projections") +
  theme(legend.position = 'none')
