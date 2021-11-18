# full analysis
source(here::here('datascript.R'))
source(here::here("sub-projects/trawl/R/Lorenz.R"))
source(here::here("sub-projects/trawl/R/Evenness.R"))

nperm = 99

year_removals = c("2007","2020","2021")

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::filter(year(Date) %ni% year_removals) %>%
  dplyr::mutate(Date = as.Date(paste0(as.character(year),"-",as.character(month),"-01"))) %>%
  dplyr::select(Date, date_id, Common_name, Abundance) %>%
  group_by(Date, date_id, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>% group_by(Date, date_id, Common_name) %>%
  # # remove trawls with less than 10 individuals 
  # dplyr::filter(sum(Abundance) >=10) %>%
  # group_by(Date, date_id, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::select(Date, date_id, everything())

TB_trawl_taxasite_long = TB_trawl_taxasite %>%
  pivot_longer(-Date:-date_id, names_to = "Common_name", values_to = "Abundance")


## Goals estimate the rarified richness through time
## Show changes in evenness and Hill #s
## Look at compositon 
## Show who is contributing
## Show patterns of most important contributors
# estimate the contribution of each species to diversity changes

#### ---- 

## Annual statistics for first paragraph outline
TB_trawl_sampling_yr = TB_trawl_data %>%
  dplyr::filter(year(Date) %ni% year_removals) %>%
  dplyr::select(Date, year, Common_name, Abundance) %>%
  group_by(year) %>%
  dplyr::summarise(trawl_days = length(unique(Date)),
                   S = length(unique(Common_name)),
                   individuals = sum(Abundance, na.rm = TRUE))

# when do we have the fewest samples
TB_trawl_sampling_yr %>% slice_min(trawl_days) %>% dplyr::select(year) %>% unlist
# 2011
# when do we have the most sampels
TB_trawl_sampling_yr %>% slice_max(trawl_days) %>% dplyr::select(year) %>% unlist

total_richness = length(unique(TB_trawl_taxasite_long$Common_name))

## Annual summaries of richness, etc.
TB_trawl_taxasite_yr <- TB_trawl_data %>%
  # exclude 2007, 2020 & 2021
  dplyr::filter(year %ni% year_removals) %>%
  dplyr::select(year, Common_name, Abundance) %>%
  group_by(year, Common_name) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  na.omit %>% #gropu_by(year) %>%
  # remove trawls with less than 20 individuals 
  # dplyr::filter(sum(Abundance) >=20) %>%
  group_by(year, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::select( year, everything())

## yearly stats
yr_chao1 = TB_trawl_taxasite_yr[,-1] %>%
  apply(., 1, mobr::calc_chao1) %>%
  bind_cols(TB_trawl_taxasite_yr$year,.) %>%
  setNames(., c("year","S_asymp")) %>%
  dplyr::mutate(Date = as.Date(paste0(year,"-01-01"), format = "%Y-%m-%d"))

S_n_yr = min(TB_trawl_sampling_yr$individuals, na.rm = TRUE)

yr_Sn = TB_trawl_taxasite_yr[,-1] %>%
  apply(., 1, function(x) mobr::rarefaction(x, effort = S_n_yr, method = "IBR", extrapolate = TRUE)) %>%
  bind_cols(TB_trawl_taxasite_yr$year,.) %>%
  setNames(., c("year","S_n")) %>%
  dplyr::mutate(Date = as.Date(paste0(year,"-01-01"), format = "%Y-%m-%d"))

yr_pct_rare =TB_trawl_taxasite_yr[,-1] %>%
  apply(., 1, function(x) pct_rare(x, 0.05)) %>%
  bind_cols(TB_trawl_taxasite_yr$year,.) %>%
  setNames(., c("year","pct_rare")) %>%
  dplyr::mutate(Date = as.Date(paste0(year,"-01-01"), format = "%Y-%m-%d"))

yr_group_stats = left_join(yr_chao1, yr_Sn) %>%
  left_join(yr_pct_rare) %>%
  dplyr::select(Date, year, everything()) %>%
  group_by(Date, year) %>%
  pivot_longer(c(S_n,S_asymp,pct_rare),names_to = "index", values_to = "median" ) %>%
  ungroup
#
TB_trawl_sampling_mth = TB_trawl_taxasite_long %>%
  dplyr::select(Date, date_id, Common_name, Abundance) %>%
  dplyr::filter(Abundance > 0) %>%
  group_by(Date, date_id) %>%
  dplyr::summarise(S = length(unique(Common_name)),
                   individuals = sum(Abundance, na.rm = TRUE))

S_n_mth = min(TB_trawl_sampling_mth$individuals, na.rm = TRUE)
## monthly stats
mth_chao1 = TB_trawl_taxasite[,-c(1,2)] %>%
  apply(., 1, mobr::calc_chao1) %>%
  bind_cols(TB_trawl_taxasite$date_id,.) %>%
  setNames(., c("date_id","S_asymp")) %>%
  left_join(TB_trawl_taxasite[,c(1,2)],.,by = 'date_id')

mth_Sn = TB_trawl_taxasite[,-c(1,2)] %>%
  apply(., 1, function(x) mobr::rarefaction(x, effort = S_n_mth, method = "IBR", extrapolate = TRUE)) %>%
  bind_cols(TB_trawl_taxasite$date_id,.) %>%
  setNames(., c("date_id","S_n")) %>%
  left_join(TB_trawl_taxasite[,c(1,2)],.,by = 'date_id')

mth_pct_rare =TB_trawl_taxasite[,-c(1,2)] %>%
  apply(., 1, function(x) pct_rare(x, 0.05)) %>%
  bind_cols(TB_trawl_taxasite$date_id,.) %>%
  setNames(., c("date_id","pct_rare")) %>%
  left_join(TB_trawl_taxasite[,c(1,2)],.,by = 'date_id')

mth_group_stats = left_join(mth_chao1, mth_Sn) %>%
  left_join(mth_pct_rare) %>%
  dplyr::select(Date, date_id, everything()) %>%
  group_by(Date, date_id) %>%
  pivot_longer(c(S_n,S_asymp, pct_rare),names_to = "index", values_to = "median" ) %>%
  ungroup

# trawl_rarefaction analysis

## Estimated mobr for monthly
# trawl_mob_in = mobr::make_mob_in(comm = TB_trawl_trawl %>% dplyr::select(-date_id, -Date),
#                                  plot_attr = TB_trawl_taxasite %>% dplyr::select(date_id, Date),
#                                  coord_names = "date_id",
#                                  latlong = FALSE)
# 
# trawl_mob_stats = mobr::get_mob_stats(trawl_mob_in, group_var = NULL, index = c("N", "S", "S_asymp", "S_n", "S_PIE", "PIE", "pct_rare"), effort_samples = S_n_mth, extrapolate = TRUE, boot_groups = FALSE, n_perm = nperm)

## Figure of the temporal trends of 

S_n_lab = paste0("S[",S_n,"]")

# rich_time <-
mth_group_stats %>%
  dplyr::filter(index == 'S_n') %>%
  ggplot(aes( x= Date, y = median))+
  scale_y_continuous(name = expression("Rarefied richness ("*italic(S[n])*")"), limits = c(0,12))+
  scale_x_date(date_breaks = "year", date_labels = "%Y")+
  geom_point()+
  geom_line()


## quantify loss and appearance of species in the timeseries
trawl_mth_turnover = purrr::map(c('appearance','disappearance','total'),
                                ~codyn::turnover(df = TB_trawl_taxasite_long %>% dplyr::select(-Date),
                                                 time.var = 'date_id',
                                                 species.var = 'Common_name',
                                                 abundance.var = 'Abundance',
                                                 metric = .x) %>%
                                  data.frame %>%
                                  dplyr::mutate(date_id = date_id,
                                                variable = .x) %>%
                                  dplyr::rename(value = !!names(.[1])) %>%
                                  bind_rows(data.frame(value = NA_real_, variable =.x))) %>%
  bind_rows %>%
  left_join(TB_trawl_taxasite[,c(1,2)],.,by = 'date_id')

gains = trawl_mth_turnover %>%
  dplyr::filter(variable == "appearance")

losses = trawl_mth_turnover %>%
  dplyr::filter(variable == "disappearance")

total = trawl_mth_turnover %>%
  dplyr::filter(variable == "total")
mean(gains$value)
mean(losses$value)
trawl_mth_turnover %>%
  na.omit %>%
  ggplot()+
  geom_line(aes(x = Date, y = value, group = variable, color = variable), size = 1.1)+
  scale_color_manual(values = viridis(n = 3, option = 'viridis'))+
  scale_y_continuous(name = 'Species turnover',limits = c(0,1), expand = c(0.01,0.01))+
  scale_x_date(date_breaks = 'year', limits = c(as.Date("2009-01-01"),as.Date("2019-12-31")),
               date_labels = "%Y")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'transparent', color = NA),
        legend.title = element_blank(),
        axis.title.x = element_blank())

# 
library(mgcv)
library(brms)
library(bayesplot)
library(bayestestR)
### monthly rarefied richness patterns
mth_Sn = mth_Sn %>%
  dplyr::mutate(month = month(Date),
                year = year(Date))

sn_bayes_full_mod = brm(bf(S_n ~ s(month, bs = 'cc', k = 12)+ s(date_id)),
                     data = mth_Sn, family = gaussian(),
                     cores = 4, seed = 42,
                     iter = 5000, warmup = 1000, thin = 10,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE),
                     file = here::here("sub-projects/trawl/derived-data/sn_full.rds"))

summary(bayes_full_mod)
plot(bayes_full_mod)
pp_check(bayes_full_mod)

sn_bayes_s_mod = brm(bf(S_n ~ s(month, bs = 'cc', k = 12)),
                     data = mth_Sn, family = gaussian(),
                     cores = 4, seed = 42,
                     iter = 5000, warmup = 1000, thin = 10,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE),
                     file = here::here("sub-projects/trawl/derived-data/sn_s.rds"))

## monthly turnover, gains, and losses
mth_turnover_wide = trawl_mth_turnover %>%
  na.omit %>%
  group_by(Date, date_id) %>%
  pivot_wider( names_from = variable, values_from = value, values_fill = NA) %>%
  dplyr::mutate(month = month(Date))

turn_bayes_full_mod = brm(bf(total ~ s(month, bs = 'cc', k = 12)+ s(date_id)),
                          data = mth_turnover_wide, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 5000, warmup = 1000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/turnover_full.rds"))
summary(turn_bayes_full_mod)

turn_bayes_s_mod = brm(bf(total ~ s(month, bs = 'cc', k = 12)),
                          data = mth_turnover_wide, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 5000, warmup = 1000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/turnover_s.rds"))

summary(turn_bayes_s_mod)

turn_bayes_null_mod = brm(bf(total ~ 1), 
                          data = mth_turnover_wide, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 5000, warmup = 1000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/turnover_null.rds"))
#### patterns of Hill diversity through time

set.seed(123)
hillpart0 = hill_taxa_parti(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0)
hillpart1 = hill_taxa_parti(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 1)
hillpart2 = hill_taxa_parti(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 2)

# pairwise hill diversity 
hill_pairs0 = hill_taxa_parti_pairwise(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0, pairs = 'unique')
hill_pairs1 = hill_taxa_parti_pairwise(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 1, pairs = 'unique')
hill_pairs2 = hill_taxa_parti_pairwise(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 2, pairs = 'unique')

# create df to test differences based on time
hill_timediff = vctrs::vec_c(hill_pairs0,hill_pairs1,hill_pairs2) %>%
  left_join(TB_trawl_taxasite %>% dplyr::select(time1 = "date_id") %>%
              dplyr::mutate(site1 = as.character(1:n())), by = 'site1') %>%
  left_join(TB_trawl_taxasite %>% dplyr::select(time2 = "date_id") %>%
              dplyr::mutate(site2 = as.character(1:n())), by = 'site2') %>%
  dplyr::mutate(diff = time2-time1)

hill_timediff %>% dplyr::filter(q == 2) %>%
  ggplot()+
  geom_point(aes(x = diff, y = TD_alpha))
  
## hill 

hill_time0 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0) %>%
  data.frame %>%  dplyr::mutate(q = "q0") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time1 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 1)%>%
  data.frame %>% dplyr::mutate(q = "q1") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time2 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 2)%>%
  data.frame %>% dplyr::mutate(q = "q2") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)

hill_time = vctrs::vec_c(hill_time0,hill_time1,hill_time2) %>%
  group_by(Date, date_id) %>%
  pivot_wider(names_from = 'q', values_from = 'value') %>%
  left_join(mth_group_stats %>% group_by(Date, date_id) %>%
              pivot_wider(names_from = "index", values_from = "median")) %>%
  dplyr::mutate(q0_eff = q0/q0,
                q1_eff = q1/q0,
                q2_eff = q2/q0)

hill_time %>%
  # dplyr::filter(q == "2") %>%
  ggplot(aes(x = Date, y = q2))+
  geom_point()+
  geom_line()+
  scale_y_continuous(name = expression("Order of diversity ("~italic(''^q*D)~")"), limits = c(0,NA))+
  scale_x_date(date_breaks = 'year', limits = c(as.Date("2009-01-01"),as.Date("2019-12-31")),
               date_labels = "%Y")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_rect(fill = 'transparent', color = NA),
        legend.title = element_blank(),
        axis.title.x = element_blank())  

mth_even = TB_trawl_taxasite %>%
  dplyr::select(-Date, -date_id) %>%
  apply(., 1, gini_even) %>%
  t() %>%
  data.frame() %>%
  bind_cols(TB_trawl_taxasite[,c('Date','date_id')], .) %>%
  setNames(., c("Date","date_id", "gini_raw", "gini_norm")) %>%
  dplyr::mutate(month = month(Date)) %>%
  # pivot_longer(-Date:-date_id, names_to = "measure", values_to = "median") %>%
  ungroup

mth_even %>%
  ggplot(aes(x = Date, y = gini_norm))+
  geom_point()+
  geom_line()+
  scale_y_continuous(name = 'Normalized Gini coefficient', limits = c(0,0.5))+
  scale_x_date(date_breaks = 'year', limits = c(as.Date("2009-01-01"),as.Date("2019-12-31")),
               date_labels = "%Y")



even_bayes_full_mod = brm(bf(gini_norm ~ s(month, bs = 'cc', k = 12)+ s(date_id)),
                          data = mth_even, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 5000, warmup = 1000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/even_full.rds"))

summary(even_bayes_full_mod)
plot(even_bayes_full_mod)


?
q2_boots = boot_sample_groups(TB_trawl_taxasite %>% dplyr::select(-Date), 
                              group_var = 'date_id',
                              qD = 2,
                              n_perm = nperm)
set.seed(123)
q1_boots = boot_sample_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly), 
                              group_var = 'year',
                              qD = 1,
                              n_perm = nperm)

qD_df = q2_boots %>% 
  # summarise q2 by year
  group_by(year) %>%
  dplyr::summarise(across(q2, list(mean = ~mean(.x, na.rm = TRUE),
                                   median = ~quantile(.x, 0.5, na.rm = TRUE),
                                   quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                   quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                   quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                   quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
  dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"),
                q_order = "q2") %>%
  bind_rows(
    q1_boots %>%
      # summarise q1 by year
      group_by(year) %>%
      dplyr::summarise(across(q1, list(mean = ~mean(.x, na.rm = TRUE),
                                       median = ~quantile(.x, 0.5, na.rm = TRUE),
                                       quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                       quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                       quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                       quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
      dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"),
                    q_order = "q1")
  ) 
# 
# annual_evenness_boots = boot_gini_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly), 
#                                 group_var = 'year',
#                                 n_perm = nperm)
# 
# annual_evenness_df = annual_evenness_boots %>%
#   dplyr::select(-boot_id) %>%
#   pivot_longer(-year, names_to = 'measure', values_to = 'value') %>%
#   group_by(year,measure) %>%
#   dplyr::summarise(across(value, list(mean = ~mean(.x, na.rm = TRUE),
#                                       median = ~quantile(.x, 0.5, na.rm = TRUE),
#                                       quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
#                                       quant25= ~quantile(.x, 0.25, na.rm = TRUE),
#                                       quant75= ~quantile(.x, 0.75, na.rm = TRUE),
#                                       quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
#   dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d")) %>%
#   ungroup

annual_evenness_df = TB_trawl_taxasite_yr %>%
  dplyr::select(-year) %>%
  apply(., 1, gini_even) %>%
  t() %>%
  data.frame() %>%
  bind_cols(TB_trawl_taxasite_yr$year, .) %>%
  setNames(., c("year", "gini_raw", "gini_norm")) %>%
  pivot_longer(-year, names_to = "measure", values_to = "median") %>%
  dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d")) %>%
  ungroup

group_stats %>%
  dplyr::filter(index == "pct_rare") -> rarity




###### TEMPORAL PATTERNS OF RICHNESS -----------
#create time-accumulation df from trawl data
total_richness = length(unique(TB_trawl_taxasite_long$Common_name));
#138L

# create time data frame to calculate accumulation of different species abundances over time
TB_trawl_data %>%
  # exclude 2007, 2020 & 2021
  dplyr::filter(year %ni% year_removals) %>%
  dplyr::filter(!duplicated(Common_name, nmax = total_richness )) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  arrange(Date) %>%
  dplyr::mutate(richness = 1:n())-> TB_time_accum

# 
TB_annual_accum <- TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% year_removals) %>%
  mutate(pres = 1) %>% group_by(year) %>%
  filter(!duplicated(Common_name)) %>%
  summarise(richness = sum(pres)) %>%
  mutate(Date = as.Date(paste(as.character(year),"01","01", sep = "-"), format = "%Y-%m-%d", origin = "2000-01-01"))


# partition beta diversity into nestedness and turnover #
# uses `betapart` package #
TB_time_df <- TB_time_accum %>%
  dplyr::select(Date, year, richness) %>%
  mutate(data_type = 'continuous') %>%
  bind_rows(TB_annual_accum %>% mutate(data_type = "annual"))

#create moving jaccard calc by year-month
# create a dataframe for time series of distance compared to first sampling 
unique_combinations <- expand.grid(date_A = unique(levels(as.factor(TB_trawl_taxasite$Date))), 
                                   Date = unique(levels(as.factor(TB_trawl_taxasite$Date)))) %>%
  filter(date_A == as.character(min(TB_trawl_taxasite$Date)))

TB_date_bray = vegdist((TB_trawl_taxasite %>% dplyr::select(-julian_date,-year:-biweekly) %>% 
                          column_to_rownames("Date")), method = "bray") %>%
  as.matrix %>% as.data.frame %>% rownames_to_column("date_A") %>% 
  pivot_longer(cols = 2:dim(.)[2], names_to = "Date", values_to = "bray") %>%
  inner_join(unique_combinations) %>% filter(date_A != Date) %>% 
  mutate(Date = gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', Date),
         Date = as.Date(paste0(Date,"-01"), format = "%Y-%m-%d"))

ggplot(TB_date_bray, aes(x = Date, y = bray)) + 
  geom_point() +
  geom_line() +
  scale_y_continuous()

# partitioning beta diversity metrics from trawl dataset
# create named list of each sampling date
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}

date_lists = TB_trawl_taxasite %>%
  named_group_split(Date) %>%
  map(., ~.x %>% select(-year, -month,-date_id, ))

beta_lists <- lapply(date_lists, function(x,r){ betapart::beta.temp(x = x, y = r, index.family = 'jaccard')}, x = date_lists[[1]]) 
beta_dates <- gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', names(beta_lists))
beta_df <- bind_rows(beta_lists) %>% mutate(Date = as.Date(paste0(beta_dates,"-01"), "%Y-%m-%d")) %>% #select(Date, everything()) %>%
  pivot_longer(-Date, names_to = "partition", values_to = "value")

#### TURNOVER ------------
## quantify loss and appearance of species in the timeseries
TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% year_removals) %>%
  group_by(year, Common_name) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  group_by(year) %>%
  pivot_wider( names_from = 'Common_name', values_from = 'Abundance', values_fill = 0) %>%
  ungroup %>%
  # dplyr::mutate(time.var = 1:n()) %>%
  pivot_longer(c(-year), names_to = 'Common_name', values_to = 'Abundance')-> TB_trawl_yr

TB_trawl_taxasite_yr = TB_trawl_yr %>%
  group_by(year) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance) %>%
  dplyr::select(-Minnow, -`Sygnathus Scovelli`)
  
trawl_yr_turnover_list = purrr::map(c('appearance','disappearance','total'),
                                     ~codyn::turnover(df = TB_trawl_yr,
                                                      time.var = 'year',
                                                      species.var = 'Common_name',
                                                      abundance.var = 'Abundance',
                                                      metric = .x) %>%
                                       data.frame %>%
                                       dplyr::mutate(Date = as.Date(paste0(unique(TB_trawl_yr$year)[-1],"-01-01"), format = "%Y-%m-%d"),
                                                     variable = .x) %>%
                                       dplyr::rename(value = !!names(.[1])) %>%
                                       bind_rows(data.frame(value = NA_real_, time.var = NA_real_, Date = as.Date("2008-01-01"), variable =.x))) %>%
                                     bind_rows

trawl_yr_turnover_list %>%
  ggplot()+
  geom_line(aes(x = Date, y = value, group = variable, color = variable), size = 1.1)+
  scale_color_manual(values = viridis(n = 3, option = 'viridis'))
  

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  dplyr::filter(year(Date) %ni% year_removals) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>% group_by(Date) %>%
  # remove trawls with less than 10 individuals 
  dplyr::filter(sum(Abundance) >=10) %>%
  group_by(Date, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::mutate(julian_date = julian(Date, origin = as.Date("2009-01-01")),
                biweekly = ceiling(julian_date/30),
                month = month(Date),
                year = year(Date)) %>%
  dplyr::select(julian_date, Date,year, month, biweekly, everything())

## Estimated qD = 2 for all communities
TB_trawl_biweekly = TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% year_removals) %>%
  dplyr::mutate(julian_date = julian(Date, origin = as.Date("2009-01-01")),
                biweekly = ceiling(julian_date/14)) %>%
  group_by(biweekly, Common_name) %>%
  dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  group_by(biweekly) %>%
  pivot_wider( names_from = 'Common_name', values_from = 'Abundance', values_fill = 0) %>%
  ungroup 

###### Biweekly 2D patterns
TB_trawl_biweekly_2D = TB_trawl_biweekly %>%
  boot_sample_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly),
                     group_var = 'biweekly',
                     qD = 2,
                     n_perm = nperm)

TB_trawl_biweekly_2D %>%
  group_by(biweekly) %>%
  dplyr::summarise(across(q2, list(mean = ~mean(.x, na.rm = TRUE),
                                   median = ~quantile(.x, 0.5, na.rm = TRUE),
                                   quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                   quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                   quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                   quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}"))-> TB_trawl_biweekly_2D_summ

  
TB_trawl_biweekly_1D = TB_trawl_biweekly %>%
  boot_sample_groups(. %>% dplyr::select(-biweekly),
                     group_var = 'biweekly',
                     qD = 1,
                     n_perm = nperm)

TB_trawl_biweekly_1D %>%
  group_by(biweekly) %>%
  dplyr::summarise(across(q1, list(mean = ~mean(.x, na.rm = TRUE),
                                   median = ~quantile(.x, 0.5, na.rm = TRUE),
                                   quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                   quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                   quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                   quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) -> TB_trawl_biweekly_1D_summ
  ggplot(TB_trawl_biweekly_1D_summ) +
  geom_line(aes(x = biweekly, y = median), color = 'blue') +
  geom_line(data = TB_trawl_biweekly_2D_summ, aes(x = biweekly, y = median), color = 'red')

w.q2 = WaveletComp::analyze.wavelet(TB_trawl_biweekly_1D_summ,
                                    my.series = 2)
c.q2 = WaveletComp::analyze.coherency(TB_trawl_biweekly_1D_summ)
w.q2.image = WaveletComp::wt.image(w.q2)
######

TB_trawl_biweekly_qD = TB_trawl_biweekly %>%
  split(., seq(1:nrow(.))) %>%
  purrr::map2(., c(0,1,2), ~hillR::hill_taxa(.x %>% dplyr::select(-year), q = .y)) %>% 
  setNames(., nm = c("q0","q1","q2")) %>%
  bind_cols() %>%
  dplyr::mutate(biweekly = TB_trawl_biweekly$biweekly)


TB_trawl_yr_qD = TB_trawl_taxasite_yr %>%
  split(nrow(.)) %>%
  purrr::map2(., c(0,1,2), ~hillR::hill_taxa(.x %>% dplyr::select(-year), q = .y)) %>% 
  setNames(., nm = c("q0","q1","q2")) %>%
  bind_cols() %>%
  dplyr::mutate(year = TB_trawl_taxasite_yr$year)

ggplot(trawl_qD, aes(x =Date, y =., group = qD, color = qD)) +
  geom_point(size = 2, alpha = 0.5)+
  geom_line(alpha = 0.5)+
  geom_smooth( se = FALSE, span = 0.1)+
  scale_y_continuous(limits = c(0,10))

q2_boots = boot_sample_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly), 
                   group_var = 'year',
                   qD = 2,
                   n_perm = nperm)

q1_boots = boot_sample_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly), 
                              group_var = 'year',
                              qD = 1,
                              n_perm = nperm)

qD_df = q2_boots %>% 
  # summarise q2 by year
  group_by(year) %>%
  dplyr::summarise(across(q2, list(mean = ~mean(.x, na.rm = TRUE),
                                   median = ~quantile(.x, 0.5, na.rm = TRUE),
                                   quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                   quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                   quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                   quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
  dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"),
                q_order = "q2") %>%
  bind_rows(
    q1_boots %>%
      # summarise q1 by year
      group_by(year) %>%
      dplyr::summarise(across(q1, list(mean = ~mean(.x, na.rm = TRUE),
                                       median = ~quantile(.x, 0.5, na.rm = TRUE),
                                       quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                       quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                       quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                       quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
      dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"),
                    q_order = "q1")
  )

qD_df %>% dplyr::filter(q_order == 'q2') %>%
  dplyr::summarise(across(median, ~mean(.x ,na.rm = TRUE))) %>% round(2)

qD_df %>%
  ggplot()+
  geom_ribbon(aes(x = Date, ymin = quant2.5, ymax = quant97.5, group = q_order, fill = q_order), alpha = 0.5)+
  geom_line(aes(x = Date, y = median, group = q_order, color = q_order), size = 1.2)+
  geom_point(aes(x = Date, y= mean, group = q_order), shape = 6)+
  scale_y_continuous(name = expression("Order of diversity ("~italic(''^q*D)~")"),limits = c(0,11))+
  scale_color_manual(labels = c("1","2"), values = viridis(n =2, option = 'viridis'))+
  scale_fill_manual(labels = c("1","2"),values = viridis(n = 2, alpha = 0.5, option = 'viridis'))+
  theme(legend.title = element_blank(),
        legend.position = c(0.01,0.998),
        legend.justification = c(0,1),
        )

qeven_df = list(q1_boots,q2_boots) %>%
  map(~.x %>% group_by(year) %>% slice_sample(n = 1e5, replace = TRUE) %>% dplyr::select(-boot_id)) %>% bind_cols %>%
  dplyr::select(year = 'year...1', q1,q2,-3) %>%
  dplyr::mutate(qeven = q1/q2) %>%
  group_by(year) %>%
  dplyr::summarise(across(qeven, list(mean = ~mean(.x, na.rm = TRUE),
                                      median = ~quantile(.x, 0.5, na.rm = TRUE),
                                      quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                      quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                      quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                      quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
    dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"))

qeven_df %>%
ggplot()+
  geom_ribbon(aes(x = Date, ymin = quant2.5, ymax = quant97.5), alpha = 0.5)+
  geom_line(aes(x = Date, y = median), size = 1.2)+
  scale_y_continuous(name = expression("Order of diversity ("~italic(''^q*D)~")"),limits = c(0,3))#+
  # scale_color_manual(labels = c("1","2"), values = viridis(n =2, option = 'viridis'))+
  # scale_fill_manual(labels = c("1","2"),values = viridis(n = 2, alpha = 0.5, option = 'viridis'))+
  # theme(legend.title = element_blank(),
  #       legend.position = c(0.01,0.998),
  #       legend.justification = c(0,1),
  # )

TB_trawl_biweekly_qD = TB_trawl_biweekly %>%
  split(., seq(1:nrow(.))) %>%
  purrr::map2(., c(0,1,2), ~hillR::hill_taxa(.x %>% dplyr::select(-year), q = .y)) %>% 
  setNames(., nm = c("q0","q1","q2")) %>%
  bind_cols() %>%
  dplyr::mutate(biweekly = TB_trawl_biweekly$biweekly)


### Calculate the diversity profile of the annual community -----------

TB_trawl_yr_evenness = TB_trawl_taxasite_yr %>%
  split(., seq(1:nrow(.))) %>%
  purrr::map(~new_fun(.x, seq(0,4,0.1), "E3")) %>%
  setNames(., nm = TB_trawl_taxasite_yr$year) %>%
  purrr::map(function(a) a %>% t %>% data.frame %>%
               setNames(., nm = seq(0,4,0.1))) %>%
          bind_rows(.id = 'year') %>%
  pivot_longer(-year, names_to = 'q', values_to = 'eff_spp') %>%
  dplyr::mutate(Date = as.Date(paste0(year,"-01-01"), format = "%Y-%m-%d")) %>%
  dplyr::mutate(across(q, as.numeric))

TB_trawl_q1_plot =
  TB_trawl_yr_evenness %>%
  dplyr::filter(q == 1) %>%
  ggplot()+
  geom_line(aes(x = Date, y = eff_spp))+
  scale_y_continuous(name = 'Eff. Species (%)', limits = c(0,0.25))+
  scale_x_date(breaks = as.Date(c("2009-01-01","2013-01-01","2019-01-01"),format = "%Y-%m-%d"),
               date_labels = "%Y")+
  theme(axis.title.x = element_blank())+
  annotate('text', x = as.Date(Inf), y = as.Date(Inf), label = expression(italic(""^1*D)), size =6,
           hjust = 1, vjust = 1)

TB_trawl_q2_plot =
  TB_trawl_yr_evenness %>%
  dplyr::filter(q == 2) %>%
  ggplot()+
  geom_line(aes(x = Date, y = eff_spp))+
  scale_y_continuous(name = 'Eff. Species (%)', limits = c(0,0.18))+
  scale_x_date(breaks = as.Date(c("2009-01-01","2013-01-01","2019-01-01"),format = "%Y-%m-%d"),
               date_labels = "%Y")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  annotate('text', x = as.Date(Inf), y = as.Date(Inf), label = expression(italic(""^2*D)), size =6,
           hjust = 1, vjust = 1)

TB_trawl_yr_evenness %>%
  ggplot() +
  geom_path(aes(x = q, y = eff_spp, group = year, color = year), size = 1.1)+
  scale_y_continuous(name = 'Effective species (%)', limits = c(0,1), expand = c(0.01,0.01))+
  scale_x_continuous(name = expression(italic(""^q*D)), limits = c(0,3)) +
  scale_color_manual(values = viridis(n = 12, option = 'viridis')) +
  annotation_custom(
    ggplotGrob(TB_trawl_q1_plot),
    xmin = 0.75, xmax = 1.9, ymin = 0.5, ymax = 1
  )+
  annotation_custom(
    ggplotGrob(TB_trawl_q2_plot),
    xmin = 2, xmax = 3, ymin = 0.5, ymax = 1
  )+
  theme(legend.title = element_blank(),
        legend.position = 'top')+
  guides(color=guide_legend(nrow=1,byrow=TRUE))
  



vp1 <- viewport(width = 0.4, height = 0.4, x = 0.33, y = 0.8)


vp2 <- viewport(width = 0.4, height = 0.4, x = 0.66, y = 0.8)

### calculate the diversity and similarity of orders of diversity ----------

# debugonce(hillR::hill_taxa_parti)
TB_trawl_taxasite_yr_q1 =  hillR::hill_taxa_parti_pairwise(comm = TB_trawl_taxasite_yr[,-1],
                         q = 1, base = exp(1), rel_then_pool = TRUE, pairs = 'unique') %>%
  dplyr::filter(site1 == 1) %>%
  bind_cols(Date = as.Date(paste0(TB_trawl_taxasite_yr$year[-1],"-01-01"), format = "%Y-%m-%d"),.)
  

TB_trawl_taxasite_yr_q2 =  hillR::hill_taxa_parti_pairwise(comm = TB_trawl_taxasite_yr[,-1],
                         q = 2, base = exp(1), rel_then_pool = TRUE, pairs = 'unique') %>%
  dplyr::filter(site1 == 1) %>%
  bind_cols(Date = as.Date(paste0(TB_trawl_taxasite_yr$year[-1],"-01-01"), format = "%Y-%m-%d"),.)
  
TB_trawl_taxasite_yr_q0 =  hillR::hill_taxa_parti_pairwise(comm = TB_trawl_taxasite_yr[,-1],
                                                             q = 0, base = exp(1), rel_then_pool = TRUE, pairs = 'unique') %>%
  dplyr::filter(site1 == 1) %>%
    bind_cols(Date = as.Date(paste0(TB_trawl_taxasite_yr$year[-1],"-01-01"), format = "%Y-%m-%d"),.)

TB_trawl_taxasite_yr_q0.5 =  hillR::hill_taxa_parti_pairwise(comm = TB_trawl_taxasite_yr[,-1],
                                                           q = 0.5, base = exp(1), rel_then_pool = TRUE, pairs = 'unique') %>%
  dplyr::filter(site1 == 1) %>%
  bind_cols(Date = as.Date(paste0(TB_trawl_taxasite_yr$year[-1],"-01-01"), format = "%Y-%m-%d"),.)

TB_trawl_taxasite_yr_qD = bind_rows(TB_trawl_taxasite_yr_q0, TB_trawl_taxasite_yr_q0.5) %>%
  bind_rows(TB_trawl_taxasite_yr_q1) %>%
  bind_rows(TB_trawl_taxasite_yr_q2) 

TB_trawl_taxasite_yr_qD %>% 
  dplyr::mutate(q = factor(q)) %>%
  ggplot() + 
  geom_point(aes(x = Date, y = local_similarity, group = q, color = q)) + 
  geom_smooth(aes(x = Date, y = local_similarity, group = q, color =q),span = 0.9, se = FALSE)


## tutorial
View(FD::dummy)


### create ordering based on relative abundances -------------
TB_trawl_order =  TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  dplyr::filter(year(Date) %ni% year_removals) %>%
  group_by(Date) %>%
  dplyr::mutate(rel_abun = Abundance/sum(Abundance, na.rm = TRUE)) %>%
  dplyr::select(-Abundance) %>%
  group_by(Common_name) %>%
  dplyr::summarise(rel_abun = mean(rel_abun, na.rm = TRUE)) %>%
  arrange(desc(rel_abun)) %>%
  na.omit %>%
  dplyr::select(Common_name) %>%
  unlist %>% unname

TB_trawl_order_names = fuzzySim::spCodes(TB_trawl_order, nchar.gen = 4, nchar.sp = 6)

TB_trawl_taxasite_yr_dis = TB_trawl_taxasite_yr %>%
  pivot_longer(-year,names_to = 'Common_name', values_to = 'Abundance') %>%
  pivot_wider(names_from = 'year', values_from = 'Abundance') %>%
  column_to_rownames('Common_name')

TB_trawl_q2_dis = dis1(TB_trawl_taxasite_yr_dis, q = 2, type = 'tax') %>%
  .[,TB_trawl_order] %>%
  data.frame %>%
  rownames_to_column('dis_type') %>%
  dplyr::filter(dis_type == "UqN") %>%
  pivot_longer(-dis_type, names_to = 'Common_name', values_to = 'value') %>%
  dplyr::mutate(Common_name = factor(Common_name, make.names(TB_trawl_order)))


TB_trawl_q1_dis = dis1(TB_trawl_taxasite_yr_dis, q = 1, type = 'tax') %>%
  .[,TB_trawl_order] %>%
  data.frame %>%
  rownames_to_column('dis_type') %>%
  dplyr::filter(dis_type == "UqN") %>%
  pivot_longer(-dis_type, names_to = 'Common_name', values_to = 'value') %>%
  dplyr::mutate(Common_name = factor(Common_name, make.names(TB_trawl_order)))

TB_trawl_q0_dis = dis1(TB_trawl_taxasite_yr_dis, q = 0, type = 'tax') %>%
  .[,TB_trawl_order] %>%
  data.frame %>%
  rownames_to_column('dis_type') %>%
  dplyr::filter(dis_type == "UqN") %>%
  pivot_longer(-dis_type, names_to = 'Common_name', values_to = 'value') %>%
  dplyr::mutate(Common_name = factor(Common_name, make.names(TB_trawl_order)))

q0_dis_plot = TB_trawl_q0_dis %>%
  ggplot()+
  geom_col(aes(x= Common_name, y = value));q0_dis_plot

q1_dis_plot = TB_trawl_q1_dis %>%
  ggplot()+
  geom_col(aes(x= Common_name, y = value));q1_dis_plot

q2_dis_plot = TB_trawl_q2_dis %>%
  ggplot()+
  geom_col(aes(x= Common_name, y = value));q2_dis_plot

### evenness -----------

TB_trawl_evenness = TB_trawl_taxasite %>%
  split(., seq(1:nrow(.))) %>%
  purrr::map(~.x[-c(1:5)] %>% 
                   unlist %>%
                   gini_even(.)) %>% 
  bind_rows() %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(julian_date:biweekly),.)

# gini_plot = 
  
  TB_trawl_evenness %>%
  ggplot()+
  geom_line(aes(x = Date, y = `Non-normalized Gini`))

# debugonce(boot_gini_groups)
annual_evenness_boots = boot_gini_groups(TB_trawl_taxasite %>% dplyr::select(-julian_date,-Date,-month,-biweekly), 
                                group_var = 'year',
                                n_perm = nperm)

annual_evenness_df = annual_evenness_boots %>%
  dplyr::select(-boot_id) %>%
  pivot_longer(-year, names_to = 'measure', values_to = 'value') %>%
  group_by(year,measure) %>%
  dplyr::summarise(across(value, list(mean = ~mean(.x, na.rm = TRUE),
                                                median = ~quantile(.x, 0.5, na.rm = TRUE),
                                                quant2.5= ~quantile(.x, 0.025, na.rm = TRUE),
                                                quant25= ~quantile(.x, 0.25, na.rm = TRUE),
                                                quant75= ~quantile(.x, 0.75, na.rm = TRUE),
                                                quant97.5= ~quantile(.x, 0.975, na.rm = TRUE)), .names = "{.fn}")) %>%
  dplyr::mutate(Date = as.Date(paste(year,"01","01", sep = "-"), format = "%Y-%m-%d"))

annual_evenness_df %>%
  dplyr::filter(measure == 'gini_raw') %>%
  ggplot()+
  geom_ribbon(aes(x = Date, ymin = quant2.5, ymax = quant97.5), alpha = 0.5) +
  geom_line(aes(x = Date, y = median), size = 1.1) +
  geom_point(aes(x = Date, y = mean), shape = 22, size = 3)

  # dplyr::mutate(gini_even = purrr::map(.,~.x %>% 
  #                                            dplyr::select(-julian_date:-biweekly) %>%
  #                                            unlist %>%
  #                                            gini_even(.)
  #                                            ))


library(bsts)

y = mth_Sn$S_n
ss <- AddLocalLevel(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 4)

sn_model = bsts(y, state.specification = ss, niter = 5000, seed = 42)
burn = SuggestBurn(0.2, sn_model)
summary(sn_model)

plot(sn_model)


components <- cbind.data.frame(
  colMeans(sn_model$state.contributions[-(1:burn),'trend',]),
  colMeans(sn_model$state.contributions[-(1:burn),'seasonal.4.1',]),
  Date = as.Date(mth_Sn$Date)) %>%
  setNames(.,c("Trend", "Seasonality", "Date")) %>%
  pivot_longer(-Date, names_to = 'variable', values_to = 'values')


ggplot(components, aes(x = Date, y = values)) +
  geom_line()+
  facet_wrap(~variable, scales = 'free_y', ncol = 1)

summary(sn_model, burn = burn)

p = predict.bsts(sn_model, horizon = 12, burn = burn, quantiles = c(0.025,0.975))

d2 = data.frame(
  c(as.numeric(-colMeans(sn_model$one.step.prediction.errors[-(1:burn),])+y),
    as.numeric(p$mean)),
  c(as.numeric(mth_Sn$S_n),rep(NA,12)),
  c(as.Date(mth_Sn$Date), rep(NA,12))
) %>%
  setNames(., c("Fitted", "Actual", "Date"))

ggplot(d2, aes( Fitted, Actual))+geom_point()

