source(here::here("datascript.R"))
source(here::here("sub-projects/trawl/R/Lorenz.R"))
source(here::here("sub-projects/trawl/R/Evenness.R"))

library(mgcv)
library(brms)
library(bayesplot)
# which years should we remove
year_removals = c("2007","2020","2021")

# create the base data set
TB_trawl_data = TB_trawl_data %>%
  dplyr::filter(year %ni% year_removals)

TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::filter(year(Date) %ni% year_removals) %>%
  dplyr::mutate(Date = as.Date(paste0(as.character(year),"-",as.character(month),"-01"))) %>%
  dplyr::select(Date, date_id, Common_name, Abundance) %>%
  group_by(Date, date_id, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>% group_by(Date, date_id, Common_name) %>%
  # # remove trawls with less than 10 individuals 
  # dplyr::filter(sum(Abundance) >=10) %>%
  # group_by(Date, date_id, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::select(Date, date_id, everything())


hill_time = function(x,...){
  hillR::hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date,-date_id), q = x) %>%
    data.frame %>% dplyr::mutate(q = paste0("q",as.character(x))) %>% setNames(.,c('value','q')) %>%
    bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
}
# debugonce(hill_time)
hill_time_df = furrr::future_map(seq(0,3, by = 0.1), ~hill_time(.x)) %>%
  bind_rows %>%
  group_by(Date, date_id) %>% 
  pivot_wider(names_from = 'q', values_from = 'value') %>%
  dplyr::mutate(across(contains('q'), list(~(.x/q0)*100), .names = "{.col}_eff")) %>%
  dplyr::mutate(month = month(Date)) %>%
  dplyr::select(date_id, month, contains('eff'))

# form_string = bf(q0.1_eff ~ date_id + s(month, bs = 'cc', k = 12))
# mod = brm(form_string,
#           data = hill_time_df,
#           family = gaussian(),
#           cores = 4, seed = 42,
#           iter = 30000,
#           warmup = 25000, thin = 10,
#           control = list(adapt_delta = 0.99),
#           save_pars = save_pars(all = TRUE),
#           file = here::here("sub-projects/trawl/derived-data/q0_eff_lin.rds"))
# summary(mod)[['fixed']]['date_id',]
# conditional_effects(mod, 'date_id')
# posterior_interval(mod, pars = 'date_id')
eff_cols = names(hill_time_df)[grepl("*eff", names(hill_time_df))]
  
model_qD = function(x,...){
  form_string = paste(x," ~ date_id + s(month, bs = 'cc', k = 12)")
  bayes_form = bf(formula = formula(form_string))
  mod = brm(bayes_form,
            data = hill_time_df,
            family = gaussian(),
            cores = 4, seed = 42,
            iter = 30000, warmup = 25000, thin = 10,
            control = list(adapt_delta = 0.99),
            save_pars = save_pars(all = TRUE),
            file = paste0(here::here("sub-projects/trawl/derived-data"),"/",as.character(x),"_full.rds"),
            file_refit = 'on_change')
  
  return(posterior_summary(mod, pars = 'date_id', probs = c(0.025,0.25,0.5,0.75,0.975)))
}

# debugonce(model_qD)
# qD_coefs = lapply(eff_cols, model_qD) %>%
# setNames(., nm = eff_cols)
# saveRDS(qD_coefs, file = here::here("sub-projects/trawl/derived-data/qD_coefs.rds"))
qD_coefs = readRDS( file = here::here("sub-projects/trawl/derived-data/qD_coefs.rds"))
qD_coefs_df = qD_coefs %>% lapply(., data.frame) %>% bind_rows(.id = "model") %>%
  dplyr::mutate(q = as.numeric(sub("^q(\\d{1}|\\d{1}.\\d{1}).*","\\1", model)))

# get_posteriors = function(mod, pars, probs,...){
  # x = readRDS(file = mod)
  # post_int = posterior_interval(x, pars = pars, prob = probs)
  # return(as.numeric(post_int))
# } 
# pars = 'date_id'
# # debugonce(get_posteriors)
# qD_post_int0.5 = list.files(path = here::here("sub-projects/trawl/derived-data/"), pattern = "q.*_eff_full.rds", full.names = TRUE) %>%
#   purrr::map2(., 0.5, ~get_posteriors(mod = .x, pars = 'date_id', prob = .y)) %>%
#   setNames(., nm = eff_cols) %>% bind_rows %>% t %>% data.frame %>%
#   rownames_to_column("model") %>%
#   setNames(., c("model","25%","75%")) %>%
#   dplyr::mutate(q = as.numeric(sub("^q(\\d{1}|\\d{1}.\\d{1}).*","\\1", model)))
# 
# qD_post_int0.9 = list.files(path = here::here("sub-projects/trawl/derived-data/"), pattern = "q.*_eff_full.rds", full.names = TRUE) %>%
#   purrr::map2(., 0.90, ~get_posteriors(mod = .x, pars = 'date_id', prob = .y)) %>%
#   setNames(., nm = eff_cols) %>% bind_rows %>% t %>% data.frame %>%
#   rownames_to_column("model") %>%
#   setNames(., c("model","5%","95%")) %>%
#   dplyr::mutate(q = as.numeric(sub("^q(\\d{1}|\\d{1}.\\d{1}).*","\\1", model)))

qD_coefs_df = qD_coefs_df %>%
  dplyr::mutate(across(c("Estimate",matches("Q\\d{1,2}")), ~.x *365)) 

qD_coefs_df %>% 
  ungroup %>%
  dplyr::filter(q != 0) %>%
  ggplot()+
  geom_ribbon(aes(x = q, ymin = Q2.5, ymax = Q97.5), color = 'lightgrey', alpha = 0.3)+
  geom_ribbon(aes(x = q, ymin = Q25, ymax = Q75), color = 'lightgrey', alpha = 0.5)+
  geom_line(aes(x = q, y = Estimate))+
  # geom_point(aes(x = q, y = Estimate))+
  scale_x_reverse(name = expression(""*q^D), limits = c(3,0), expand = c(0.001,0.001))+
  scale_y_reverse(name = expression('Temporal trend ('*Delta*'%/year)'), limits = c(0,-3), 
                  expand = c(0.001,0.001))+
  coord_flip()
  

ccf(hill_time_df$q0.1_eff, hill_time_df$q1_eff)


hill_time_long = hill_time_df %>%
  dplyr::select(-month) %>%
  group_by(Date, date_id) %>%
  pivot_longer(-Date:-date_id, names_to = 'mode', values_to = 'spp_eff')

hill_time_long %>%
  dplyr::filter(mode != "q0_eff") %>%
  ggplot()+
  geom_line(aes(x = date_id, y = spp_eff, group = mode, color = mode))+
  geom_point(aes(x = date_id, y = spp_eff, group = mode, color = mode))


paste0(here::here("sub-projects/trawl/derived-data/"),"/",as.character(eff_cols[1]),"_full.rds")

hill_time0 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0) %>%
  data.frame %>%  dplyr::mutate(q = "q0") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time0.1 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0.1) %>%
  data.frame %>%  dplyr::mutate(q = "q0.1") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time0.5 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 0.5) %>%
  data.frame %>%  dplyr::mutate(q = "q0.5") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time1 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 1)%>%
  data.frame %>% dplyr::mutate(q = "q1") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)
hill_time2 = hill_taxa(TB_trawl_taxasite %>% dplyr::select(-Date, -date_id), q = 2)%>%
  data.frame %>% dplyr::mutate(q = "q2") %>% setNames(., c('value','q')) %>%
  bind_cols(TB_trawl_taxasite %>% dplyr::select(Date, date_id),.)

hill_time = vctrs::vec_c(hill_time0,hill_time0.1,hill_time0.5,hill_time1,hill_time2) %>%
  group_by(Date, date_id) %>%
  pivot_wider(names_from = 'q', values_from = 'value') %>%
  left_join(mth_group_stats %>%
              dplyr::mutate(month = month(Date)) %>%
              group_by(Date, date_id) %>%
              pivot_wider(-month, names_from = "index", values_from = "median")) %>%
  dplyr::mutate(q0_eff = q0/q0,
                q0.1_eff = (q0.1/q0)*100,
                q0.5_eff = (q0.5/q0)*100,
                q1_eff = (q1/q0)*100,
                q2_eff = (q2/q0)*100,
                q1_sn = q1/S_n,
                q2_sn = q2/S_n)