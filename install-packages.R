# Install pacman if it isn't already installed
if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")

# Install packages required for analysis
# Add any packages required for your data cleaning and manipulation to the
# `p_load` function call below
# Do not remove the packages already listed here
# they are important for running the livedat repository
# OG package list: pacman::p_load(git2r, httr, semver, testthat, yaml)

pacman::p_load(git2r, httr, semver, testthat, yaml, data.table, 
               RCurl, plyr, tidyverse, furrr, googlesheets4,
               tictoc, chron, lubridate, httr, TTR, grid,
               gridExtra, ggridges, iNEXT, vegan,
               viridis, broom, bbmle, ggthemes, ggeffects)

pacman::p_load_gh("jimjunker1/junkR")

theme_mod <<- theme_bw() %+replace% theme(panel.grid = element_blank())
theme_set(theme_mod)
yday_summary <<- function(data,...){
  if("yday" %in% colnames(data)){
    data_summ = data %>%
      group_by(yday) %>%
      summarise_if(is.numeric, list(~na.rm_mean(.,na.rm = TRUE), ~max(.,na.rm = TRUE), ~min(., na.rm = TRUE))) %>%
      ungroup()
  } ifelse{
    
  }
    
}
