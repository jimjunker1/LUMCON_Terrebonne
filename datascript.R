source(here::here("packages.R"))
# call different data import scripts depending on changes to files.

# import biology data if it has been greater than 1 month since download
load(here::here("data/biology_sysDate.rda"))
ifelse((Sys.Date()-biology_sysDate) >= 30, source(here::here("R/import_bioData.R")),source(here::here("R/load_bioData.R")))
rm(biology_sysDate)

# import environmental data
source(here::here("R/load_enviroData.R"))
