# run evenness models for GoF measures
even_bayes_full_mod = brm(bf(gini_norm ~ date_id + s(month, bs = 'cc', k = 12)),
                          data = mth_even, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 7000, warmup = 5000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/even_full.rds"),
                          file_refit = 'on_change')
even_bayes_full_mod <- add_criterion(even_bayes_full_mod, c("loo"),
                                     file = here::here("sub-projects/trawl/derived-data/even_full.rds"))

# summary(even_bayes_full_mod)
# plot(even_bayes_full_mod)
#estimate the probability of being zero or below
date_draws = as.data.frame(even_bayes_full_mod)$b_date_id
percentile <- ecdf(date_draws)

# estimate 
# msms <- conditional_smooths(even_bayes_full_mod)
# plot(msms)
ceff <- conditional_effects(even_bayes_full_mod, 'date_id')
# plot(ceff)

even_bayes_s_mod = brm(bf(gini_norm ~ s(month, bs = 'cc', k = 12)),
                          data = mth_even, family = gaussian(),
                          cores = 4, seed = 42,
                          iter = 7000, warmup = 5000, thin = 10,
                          control = list(adapt_delta = 0.99),
                          save_pars = save_pars(all = TRUE),
                          file = here::here("sub-projects/trawl/derived-data/even_s.rds"),
                          file_refit = 'on_change')
even_bayes_s_mod <- add_criterion(even_bayes_s_mod, c("loo"),
                                     file = here::here("sub-projects/trawl/derived-data/even_s.rds"))

even_bayes_null_mod = brm(bf(gini_norm ~ 1),
                       data = mth_even, family = gaussian(),
                       cores = 4, seed = 42,
                       iter = 7000, warmup = 5000, thin = 10,
                       control = list(adapt_delta = 0.99),
                       save_pars = save_pars(all = TRUE),
                       file = here::here("sub-projects/trawl/derived-data/even_null.rds"),
                       file_refit = 'on_change')
even_bayes_null_mod <- add_criterion(even_bayes_null_mod, c("loo"),
                                  file = here::here("sub-projects/trawl/derived-data/even_null.rds"))

loo_compare(even_bayes_full_mod, even_bayes_s_mod, even_bayes_null_mod, criterion = c("loo"), model_names = NULL)
