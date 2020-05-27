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
  select(Date, Common_name, Abundance) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() 

seasonal_time <- data.frame(day = 1:365, sin_year = sin(1:365), cos_year = cos(1:365))

# write.csv(TB_trawl_data %>% select(Common_name, species_mod) %>% group_by(Common_name, species_mod) %>% unique, file = "./data/spp_abbrev.csv", quote = FALSE, row.names = FALSE)

#create covariate table
TB_trawl_covariates <- TB_trawl_taxasite %>% select(Date) %>% mutate(Date = as.Date(Date),
                                                                     time = as.integer(Date, origin = min(Date)),
                                                                     day = lubridate::yday(Date)) %>%
  left_join(seasonal_time) %>% data.frame

TB_trawl_taxasite %>% select(-Date) %>% data.frame -> TB_trawl_taxa

# combine taxa data and covariate data to list, name, and conform to LDATS inputs
TB_LDATS <- list(TB_trawl_taxa, TB_trawl_covariates) %>%
  setNames(.,nm = c("document_term_table","document_covariate_table")) %>%
  LDATS::conform_LDA_TS_data(quiet = FALSE)

debugonce(check_LDA_TS_inputs)
check_LDA_TS_inputs(TB_LDATS, topics = 2:5, nseeds = 2, formulas = ~time, nchangepoints = 0:1, timename = "time", 
                    weights = TRUE, control = list(quiet = FALSE))
#should return NULL if everything is okay

lda_model_set <- LDA_set(document_term_table = TB_LDATS$document_term_table,
                         topics = c(2:5),
                         nseeds = 2,
                         control = list(quiet = FALSE))

# saveRDS(lda_model_set, file = "./data/lda_model_set.RDS")
lda_model_set <- readRDS(file = "./data/lda_model_set.rds")
selected_lda_set <- select_LDA(lda_model_set)
plot(selected_lda_set[[1]])

TB_trawl_ts <- TS_on_LDA(LDA_models = selected_lda_set, 
          document_covariate_table = TB_LDATS$document_covariate_table,
          formulas = ~sin_year + cos_year,
          nchangepoints = c(0:5),
          timename = "time",
          weights = document_weights(TB_LDATS$document_term_table),
          control = list(nit = 1000))

# saveRDS(TB_trawl_ts, file = "./data/TB_trawl_ts.rds")
TB_trawl_ts <- readRDS(file = "./data/TB_trawl_ts.rds")

selected_ts <- select_TS(TB_trawl_ts)
plot(selected_ts)


plot(rodents$document_covariate_table$newmoon, rodents$document_covariate_table$sin_year, type = "l")
lines(rodents$document_covariate_table$newmoon, rodents$document_covariate_table$cos_year, type = "l", col = "red")


trawl_LDATS <- LDA_TS(TB_LDATS, topics = 8:12, nseeds = 16, formulas = ~1, nchangepoints = 0:2, timename = "time",
                      weights = TRUE)
