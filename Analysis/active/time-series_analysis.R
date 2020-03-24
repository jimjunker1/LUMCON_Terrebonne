#time series analysis
source("install-packages.R")#installs all necessary packages
source("datascript.R")#imports data with some minor cleanup


hydro_summ <- hydro_df %>%
  na.omit() %>%
  group_by(yday) %>% 
  summarise_at(vars(temp_C, temp_F, spCon_uS_cm, sal_psu), list(~na.rm_mean(., na.rm = TRUE), ~max(., na.rm = TRUE), ~min(.,na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(DATE = as.Date(yday-1,  origin = "2020-01-01"))

do_summ <- do_df %>%
  na.omit() %>%
  group_by(yday) %>%
  summarise_at(vars(do_mg_L), list(~na.rm_mean(.,na.rm = TRUE), ~max(.,na.rm = TRUE), ~min(., na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(DATE = as.Date(yday-1, origin = "2020-01-01"))