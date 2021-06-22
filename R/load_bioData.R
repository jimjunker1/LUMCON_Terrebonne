# load data

## trawl data
load_bioData <- function(){TB_trawl_data <<- readRDS(file = "./data/TB_trawl_data.rds")}
ifelse(file.exists("./data/TB_trawl_data.rds"), suppressWarnings(load_bioData()), source("./R/import_bioData.R"));rm(load_bioData)
# TB_trawl_data <<- readRDS(file = "./data/TB_trawl_data.rds")