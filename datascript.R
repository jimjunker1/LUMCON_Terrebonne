source(here::here("packages.R"))
here::i_am("datascript.R")
# call different data import scripts depending on changes to files.

# import biology data if it has been greater than 1 month since download
load(here::here("data/biology_sysDate.rda"))
ifelse((Sys.Date()-biology_sysDate) >= 30, source(here::here("sub-projects/trawl/R/import_bioData.R")),source(here::here("sub-projects/trawl/R/load_bioData.R")))
rm(biology_sysDate)

# important model variables for referencing
iter = 7000
warmup = 5000
thin = 10
adapt_delta = 0.99


# import environmental data
# source(here::here("R/load_enviroData.R"))
