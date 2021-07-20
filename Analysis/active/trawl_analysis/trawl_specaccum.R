
# Full species accumulation curve (SAC)

##create taxa by site matrix
TB_trawl_commonsite <- TB_trawl_data %>%
  dplyr::select(year, month, date_id, Common_name) %>%
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

#create moving jaccard calc by year-month
# create a dataframe for time series of distance compared to first sampling 
unique_combinations <- expand.grid(date_A = unique(levels(as.factor(TB_trawl_commonsite$date_id))), 
                                   Date = unique(levels(as.factor(TB_trawl_commonsite$date_id)))) %>%
  filter(date_A == as.character(min(TB_trawl_commonsite$date_id)))

TB_date_jaccard = vegdist((TB_trawl_commonsite %>% select(-year, -month) %>% 
                             column_to_rownames("date_id")), method = "jaccard") %>%
  as.matrix %>% as.data.frame %>% rownames_to_column("date_A") %>% 
  pivot_longer(cols = 2:dim(.)[2], names_to = "Date", values_to = "jaccard") %>%
  inner_join(unique_combinations) %>% filter(date_A != Date) %>% 
  mutate(Date = gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', Date),
         Date = as.Date(paste0(Date,"-01"), format = "%Y-%m-%d"))

ggplot(TB_date_jaccard, aes(x = Date, y = jaccard)) + 
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
                       
source("./figures/diversity_time_plots.R")



  