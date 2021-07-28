# load data

## trawl data
load_bioData <- function(){TB_trawl_data <<- readRDS(file = here::here("data/TB_trawl_data.rds"))}
ifelse(file.exists(here::here("data/TB_trawl_data.rds")), suppressWarnings(load_bioData()), source(here::here("analyses/trawl/R/import_bioData.R")));rm(load_bioData)
# TB_trawl_data <<- readRDS(file = "./data/TB_trawl_data.rds")