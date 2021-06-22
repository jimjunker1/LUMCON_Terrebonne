source("packages.R")
# call different data import scripts depending on changes to files.

# import biology data if it has been greater than 1 month since download
load("./data/biology_sysDate.rda")
ifelse((Sys.Date()-biology_sysDate) >= 30, source("./R/import_bioData.R"),source("./R/load_bioData.R"))
rm(biology_sysDate)

# import environmental data
source("./R/load_enviroData.R")
