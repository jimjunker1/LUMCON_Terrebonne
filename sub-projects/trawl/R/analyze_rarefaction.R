source(here::here('datascript.R'))
here::i_am("sub-porjects/trawl/R/analyze_rarefaction.R")

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  dplyr::filter(year(Date) %ni% c("2020","2021")) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>% group_by(Date) %>%
  # remove trawls with less than 10 individuals 
  dplyr::filter(sum(Abundance) >=10) %>%
  group_by(Date, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::mutate(julian_date = julian(Date, origin = as.Date("2007-01-05")),
                biweekly = ceiling(julian_date/30),
                month = month(Date),
                year = year(Date)) %>%
  dplyr::select(julian_date, Date,year, month, biweekly, everything())


TB_trawl_sampling = TB_trawl_data %>%
  dplyr::filter(year(Date) %ni% c("2020","2021")) %>%
  dplyr::select(Date, year, Common_name, Abundance) %>%
  group_by(year) %>%
  dplyr::summarise(trawl_days = length(unique(Date)),
                   S = length(unique(Common_name)),
                   individuals = sum(Abundance, na.rm = TRUE))

min(TB_trawl_sampling$year, na.rm = TRUE)

TB_trawl_sampling %>% slice_min(trawl_days) %>% dplyr::select(year) %>% unlist
 ## Explore the (dis)similarities at different timescales



# trawl_rarefaction analysis
r = 2
trawl_sum = rowSums(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly))
n_ref <- round(min(r * trawl_sum[trawl_sum > 0]))


full_stats = mobr::calc_biodiv(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly), 
                  groups = TB_trawl_taxasite$year,
                  index = c("S_n","S_asymp","S_PIE","S","N"),
                  effort = n_ref,
                  extrapolate = T,
                  return_NA = F)

## Run the mobr to estimate w

## Estimated mobr for biweekly  for all communities
TB_trawl_biweekly = TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% c("2020","2021")) %>%
  dplyr::mutate(julian_date = julian(Date, origin = as.Date("2007-01-05")),
                biweekly = ceiling(julian_date/14)) %>%
  group_by(year,biweekly, Common_name) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  group_by(year,biweekly) %>%
  pivot_wider( names_from = 'Common_name', values_from = 'Abundance', values_fill = 0) %>%
  ungroup 


trawl_mob_in = mobr::make_mob_in(comm = TB_trawl_biweekly %>% dplyr::select(-biweekly),
                                 plot_attr = TB_trawl_biweekly %>% dplyr::select(year,biweekly),
                                 coord_names = "biweekly",
                                 latlong = FALSE)


x = mobr::plot_rarefaction(trawl_mob_in, group_var = 'year', method = 'sSBR')
  scale_color_manual(values = viridis(n = 12, option ='viridis'))

trawl_mob_stats = mobr::get_mob_stats(trawl_mob_in, group_var = 'year', index = c("N", "S", "S_asymp", "S_n", "S_PIE", "PIE"), effort_samples = n_year, extrapolate = TRUE, boot_groups = TRUE, n_perm = 199)
# trawl_mob_stats = mobr::get_mob_stats(trawl_mob_in, group_var = 'year', index = c("N", "S", "S_asymp", "S_n", "S_PIE", "PIE"), extrapolate = TRUE)
x = trawl_mob_stats$groups_stats

trawl_mob_stats$groups_stats %>% dplyr::filter(index == "S_PIE") %>% dplyr::select(group, index, median) %>%
  dplyr::mutate(Date = as.Date(paste(group,"01","01", sep = "-"), format = "%Y-%m-%d")) %>%
  ggplot() +
  geom_point(aes(x = Date, y = median))+
  geom_path(aes(x = Date, y = median))+
  scale_y_continuous(limits =c(0,10))

trawl_mob_stats$groups_stats %>% dplyr::filter(index == "S_asymp") %>% dplyr::select(group, index, median) %>%
  dplyr::mutate(Date = as.Date(paste(group,"01","01", sep = "-"), format = "%Y-%m-%d")) %>%
  ggplot() +
  geom_point(aes(x = Date, y = median))+
  geom_path(aes(x = Date, y = median))

trawl_mob_stats$groups_stats %>% dplyr::filter(index == "S_n") %>% dplyr::select(group, index, effort, median, lower, upper) %>%
  dplyr::mutate(Date = as.Date(paste(group,"01","01", sep = "-"), format = "%Y-%m-%d")) %>%
  ggplot() +
  geom_ribbon(aes(x = Date, ymin = lower, ymax = upper, group = effort))+
  geom_point(aes(x = Date, y = median, group = effort))+
  geom_path(aes(x = Date, y = median, group = effort))


mobr::plot_rarefaction(trawl_mob_in,'year', method = "IBR", lwd = 3, leg_loc = 'topright' )
mobr::plot_rarefaction(trawl_mob_in,'year', method = "SBR", lwd = 3, leg_loc = 'topright' )
mobr::plot_rarefaction(trawl_mob_in,'year', method = "nsSBR", lwd = 3, leg_loc = 'topright' )
mobr::plot_rarefaction(trawl_mob_in,'year', method = "sSBR", lwd = 3, leg_loc = 'topright' )



