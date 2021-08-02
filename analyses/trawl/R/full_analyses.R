# full analysis
source(here::here('datascript.R'))
library(hillR)


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

## check sample-level statistics
TB_trawl_sampling = TB_trawl_data %>%
  dplyr::filter(year(Date) %ni% c("2020","2021")) %>%
  dplyr::select(Date, year, Common_name, Abundance) %>%
  group_by(year) %>%
  dplyr::summarise(trawl_days = length(unique(Date)),
                   S = length(unique(Common_name)),
                   individuals = sum(Abundance, na.rm = TRUE))

min(TB_trawl_sampling$year, na.rm = TRUE)
# 2007

# when do we have the fewest samples
TB_trawl_sampling %>% slice_min(trawl_days) %>% dplyr::select(year) %>% unlist
# 2011

#create time-accumulation df from trawl data
total_richness = length(unique(TB_trawl_data$Common_name));
#138L

# create time data frame to calculate accumulation of different species abundances over time
TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% c("2020","2021")) %>%
  dplyr::filter(!duplicated(Common_name, nmax = total_richness )) %>%
  dplyr::mutate(Date = as.Date(Date)) %>%
  arrange(Date) %>%
  dplyr::mutate(richness = 1:n())-> TB_time_accum

# 
TB_annual_accum <- TB_trawl_data %>%
  # exclude 2020 & 2021
  dplyr::filter(year %ni% c("2020","2021")) %>%
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

date_lists = TB_trawl_commonsite %>%
  named_group_split(date_id) %>%
  map(., ~.x %>% select(-year, -month,-date_id))

beta_lists <- lapply(date_lists, function(x,r){ betapart::beta.temp(x = x, y = r, index.family = 'jaccard')}, x = date_lists[[1]]) 
beta_dates <- gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', names(beta_lists))
beta_df <- bind_rows(beta_lists) %>% mutate(Date = as.Date(paste0(beta_dates,"-01"), "%Y-%m-%d")) %>% #select(Date, everything()) %>%
  pivot_longer(-Date, names_to = "partition", values_to = "value")


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

## Estimated qD = 2 for all communities
biweekly_sampling = TB_trawl_taxasite %>%
  group_by(biweekly) %>%
  dplyr::summarise(n = n()) %>%
  right_join(TB_trawl_taxasite %>% dplyr::select(Date, ))

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
