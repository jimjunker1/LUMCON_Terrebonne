source("datascript.R")

devtools::install_github("weecology/LDATS")
library(LDATS)
data(rodents)

r_LDATS <- LDA_TS(rodents, topics = 2:5, nseeds = 2, formulas = ~1, nchangepoints = 0:1, timename = "newmoon")

# give a go on trawl data set.
# need to keep catches in raw abundance. 
# the TS models use the total catch information 
# as weights for the quality of information

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() 

seasonal_time <- data.frame(day = 1:365, sin_year = sin(1:365), cos_year = cos(1:365))

# write.csv(TB_trawl_data %>% dplyr::select(Common_name, species_mod) %>% group_by(Common_name, species_mod) %>% unique, file = "./data/spp_abbrev.csv", quote = FALSE, row.names = FALSE)

#create covariate table
TB_trawl_covariates <- TB_trawl_taxasite %>% dplyr::select(Date) %>% mutate(Date = as.Date(Date ,origin = min(Date)),
                                                                     time = as.integer(Date, origin = min(Date)),
                                                                     day = lubridate::yday(Date)) %>%
  left_join(seasonal_time) %>% data.frame

TB_trawl_taxasite %>% dplyr::select(-Date) %>% data.frame -> TB_trawl_taxa

# combine taxa data and covariate data to list, name, and conform to LDATS inputs
TB_LDATS <- list(TB_trawl_taxa, TB_trawl_covariates) %>%
  setNames(.,nm = c("document_term_table","document_covariate_table")) %>%
  LDATS::conform_LDA_TS_data(quiet = FALSE)

# debugonce(check_LDA_TS_inputs)
check_LDA_TS_inputs(TB_LDATS, topics = 2:5, nseeds = 2, formulas = ~time, nchangepoints = 0:1, timename = "time", 
                    weights = TRUE, control = list(quiet = FALSE))
#should return NULL if everything is okay

lda_model_set <- LDA_set(document_term_table = TB_LDATS$document_term_table,
                         topics = c(2:7),
                         nseeds = 5,
                         control = list(quiet = TRUE))

# saveRDS(lda_model_set, file = "./data/lda_model_set.RDS")
lda_model_set <- readRDS(file = "./data/lda_model_set.rds")
selected_lda_set <- select_LDA(lda_model_set)
plot(selected_lda_set[[1]])

TB_trawl_ts <- TS_on_LDA(LDA_models = selected_lda_set, 
          document_covariate_table = TB_LDATS$document_covariate_table,
          formulas = ~sin_year + cos_year,
          nchangepoints = c(1:5),
          timename = "time",
          weights = document_weights(TB_LDATS$document_term_table),
          control = list(nit = 1000))

# saveRDS(TB_trawl_ts, file = "./data/TB_trawl_ts.rds")
TB_trawl_ts <- readRDS(file = "./data/TB_trawl_ts.rds")

selected_ts <- select_TS(TB_trawl_ts)
# debugonce(check_changepoints)
plot(selected_ts)

selected_ts$nchangepoints


trawl_LDATS <- LDA_TS(TB_LDATS,
                      topics = 4,
                      nseeds = 10,
                      formulas = ~sin_year + cos_year + time,
                      nchangepoints = 1:2,
                      timename = "time")
plot(trawl_LDATS)

saveRDS(trawl_LDATS, file = "./data/TB_trawl_LDATS.rds")
trawl_LDATS <- readRDS(file = "./data/TB_trawl_LDATS.rds")

plot(trawl_LDATS, option = "B")

trawl_LDATS_nochange <- LDA_TS(data = TB_LDATS,
                               topics = 5:7,
                               nseeds = 10,
                               formulas = ~sin_year + cos_year, 
                               nchangepoints = 0:4,
                               timename = "time",
                               control = list(nit = 1000))

saveRDS(trawl_LDATS_nochange, file = "./data/TB_trawl_LDATS_nochange.rds")
trawl_LDATS_nochange <- readRDS(file = "./data/TB_trawl_LDATS_nochange.rds")
plot(trawl_LDATS_nochange)

trawl_LDATS_noseason <- LDA_TS(data = TB_LDATS,
                               topics = 5:7,
                               nseeds = 10,
                               formulas = ~time,
                               nchangepoints = 1:4,
                               timename = 'time',
                               control = list(nit = 1000))

# saveRDS(trawl_LDATS_noseason, file = "./data/TB_trawl_LDATS_noseason.rds")
trawl_LDATS_noseason <- readRDS(file = "./data/TB_trawl_LDATS_noseason.rds")
plot(trawl_LDATs_noseason)

trawl_LDATS_noseason_nochange <- LDA_TS(data = TB_LDATS,
                               topics = 5:7,
                               nseeds = 10,
                               formulas = ~time,
                               nchangepoints = 0:4,
                               timename = 'time',
                               control = list(nit = 1000))

# saveRDS(trawl_LDATS_noseason_nochange, file = "./data/TB_trawl_LDATS_noseason_nochange.rds")
plot(trawl_LDATS_noseason_nochange)


##mapping across parameter space ##
params_list <-data.frame(expand.grid(topics = 2:8, nchangepoints = 0:2)) %>%
  split(., seq(nrow(.)))

# null
trawl_LDATS_null <- 
  lapply(params_list, function(x){ LDA_TS(data = TB_LDATS,
                                          topics = x$topics,
                                          nchangepoints = x$nchangepoints,
                                          nseeds = 10,
                                          formulas = ~1,
                                          timename = 'time',
                                          control = list(nit = 1000))})
# saveRDS(trawl_LDATS_null, "./data/trawl_LDATS_null.rds")
trawl_LDATS_null <- readRDS(file = "./data/trawl_LDATS_null.rds")

# time
trawl_LDATs_time <- lapply(params_list, function(x){ LDA_TS(data = TB_LDATS,
                                                            topics = x$topics,
                                                            nchangepoints = x$nchangepoints,
                                                            nseeds = 10,
                                                            formulas = ~time,
                                                            timename = 'time',
                                                            control = list(nit = 1000))})
# saveRDS(trawl_LDATs_time, file = "./data/trawl_LDATS_timelist.rds")
trawl_LDATS_time <- readRDS(file = "./data/trawl_LDATS_timelist.rds")

# season
trawl_LDATS_seasons <- 
  lapply(params_list, function(x){ LDA_TS(data = TB_LDATS,
                                          topics = x$topics,
                                          nchangepoints = x$nchangepoints,
                                          nseeds = 10,
                                          formulas = ~sin_year+cos_year,
                                          timename = 'time',
                                          control = list(nit = 1000))})
# saveRDS(trawl_LDATS_seasons, file = "./data/trawl_LDATS_seasons.rds")
trawl_LDATS_seasons <- readRDS(file ="./data/trawl_LDATS_seasons.rds")


#plotting the change in AIC in both LDA and TS with k and nchangepoints
LDATS_summary <- function(x){
  data.frame(
    k = x$'Selected LDA model'$k@k,
    LDA_AIC = min(vapply(x$'LDA models', AIC, 0, USE.NAMES = FALSE)),
    nchangepoints = x$'Selected TS model'$nchangepoints,
    TS_AIC = AIC(x$'Selected TS model'))
}

trawl_LDATS_df <- 
  lapply(trawl_LDATS_null, LDATS_summary) %>% bind_rows %>% mutate(method = "null") %>%
  bind_rows(lapply(trawl_LDATs_time, LDATS_summary) %>% bind_rows %>% mutate(method = "time")) %>%
  bind_rows(lapply(trawl_LDATS_seasons, LDATS_summary) %>% bind_rows %>% mutate(method = "season"))

#plotting the change in AIC in both LDA and TS with k and nchangepoints
ggplot(trawl_LDATS_df, aes(x = nchangepoints, y = TS_AIC, group = factor(k))) + 
  geom_line(aes(color = factor(k))) + facet_grid(~method)

ggplot(trawl_LDATS_df, aes(x = k, y = LDA_AIC))+
  geom_line(aes()) + facet_grid(~method)


