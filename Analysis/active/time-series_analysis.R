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

ggplot(hydro_df, aes(x = time, y = sal_psu)) + geom_line()

ggplot(hydro_df, aes(x = time, y = temp_F)) + geom_line()

env_long <- hydro_df %>% select(-yday, -temp_F) %>% bind_rows(temp_df %>% mutate(temp_C = tempF_to_C(temp_F)) %>% select(-yday, -temp_F)) %>%
  full_join(do_df %>% select(-yday)) %>%
  # full_join(temp_df %>% mutate(temp_C = tempF_to_C(temp_F))%>% select(-yday, - temp_F), all = TRUE) %>%
  pivot_longer(cols = temp_C:do_mg_L, names_to = "env_var", values_to = "value")

source("./figures/env_facet_plot.R");env_plot()
