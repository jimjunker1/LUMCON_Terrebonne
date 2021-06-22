source("datascript.R")

env_series = salinity_series %>% dplyr::select(Date, MC_sal, TB_sal, sss) %>%
  left_join(temperature_series %>% dplyr::select(Date, MC_temp, TB_temp, sst), by = "Date") %>%
  dplyr::mutate(julian = julian(Date, origin = as.Date("2007-01-01")),
                weekly = paste0(year(Date),"-", week(Date))) %>%
  group_by(weekly) %>% dplyr::select(-Date) %>% dplyr::mutate(Date = min(`julian`)) %>%
  group_by(Date, weekly) %>%
  dplyr::summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), NA)))

# estimate univariate state space model

sal_arss = MARSS::MARSS()

set.seed(42);sal_forecast1 = forecast::auto.arima(env_series$TB_sal, xreg = as.matrix(env_series[,c(3)]), stepwise = FALSE, approximation = FALSE)
# sal_forecast1
# set.seed(42);sal_forecast2 = forecast::auto.arima(env_series$TB_sal, xreg = as.matrix(env_series[,c(3,5)]), stepwise = FALSE, approximation = FALSE)
# sal_forecast2
TB_salinity_forecast1 = forecast::forecast(sal_forecast1, xreg = as.matrix(env_series[,3]))
# TB_salinity_forecast2 = forecast::forecast(sal_forecast2, xreg = as.matrix(env_series[,c(3,5)]))
# env_series = env_series %>% bind_cols(forecast1_mean = TB_salinity_forecast$mean) %>%
#   bind_cols(forecast1_up = TB_salinity_forecast$upper) %>% bind_cols(forecast1_dw = TB_salinity_forecast$lower) %>%
kr <- KalmanRun(env_series$TB_sal, sal_forecast1$model)
    
id.na <- which(is.na(env_series$TB_sal))
# for (i in id.na)
#   y[i] <- sal_forecast1$model$Z %*% kr$states[i,]
# alternative to the explicit loop above
y = sapply(id.na, FUN = function(x, Z, alpha) Z %*% alpha[x,], 
       Z = sal_forecast1$model$Z, alpha = kr$states)
y[id.na]



sal_ss <- bsts::AddStaticIntercept(list(), env_series$TB_sal)
sal_ss <- bsts::AddAutoAr(sal_ss, env_series$TB_sal, lags = 10)
# sal_ss <- bsts::AddAut(sal_ss, env_series$TB_sal, lags = 10)
sal_ss <- bsts::AddSeasonal(sal_ss, env_series$TB_sal, nseasons = 6)
sal_ss <- bsts::AddDynamicRegression(sal_ss, env_series$TB_sal ~ env_series$MC_sal, na.action = na.exclude)

model1 = bsts::bsts(TB_sal ~ MC_sal,
                    state.specification = sal_ss,
                    niter = 2.5e3,
                    data = env_series,
                    expected.model.size = 4,
                    seed = 123,
                    na.action = na.exclude)
summary(model1)
plot(model1)

model2 = bsts::bsts(TB_sal ~ sss + MC_sal,
                    family = "gaussian",
                    state.specification = sal_ss,
                    niter = 2.5e3, 
                    data = env_series,
                    seed = 123,
                    na.action = na.exclude)

model3 = bsts::bsts(TB_sal ~ sss, 
                    family="gaussian",
                    state.specification = sal_ss,
                    niter = 2.5e3,
                    data = env_series,
                    seed = 123,
                    na.action = na.exclude)
summary(model3)

MC_pred = predict(model1, newdata = env_series, burn = 200, na.action = na.exclude)

MC_fit = env_series %>% dplyr::select(Date, weekly, TB_sal, MC_sal) %>%
  dplyr::filter(!is.na(MC_sal)) %>% 
  bind_cols(MC_pred_mean = MC_pred$mean,
            MC_pred_lower = MC_pred$interval[1,],
            MC_pred_upper = MC_pred$interval[2,])
lines(MC_fit$MC_pred_mean)
plot(MC_fit$TB_sal, MC_fit$MC_pred_mean)
sss_pred = predict(model2, newdata = env_series, burn = 200, na.action = na.exclude)
  
sss_fit = env_series %>% dplyr::select(Date, weekly, TB_sal, sss) %>%
  dplyr::filter(!is.na(sss)) %>% 
  bind_cols(sss_pred_mean = sss_pred$mean,
            sss_pred_lower = sss_pred$interval[1,],
            sss_pred_upper = sss_pred$interval[2,])  
  

plot(sss_fit$TB_sal, )
salinity_fit = env_series %>% 
  dplyr::select(Date, weekly, TB_sal, MC_sal, sss) %>%
  bind_rows()
  

ggplot(salinity_series, aes(x = Date, y = TB_sal_series)) + geom_line(color = 'blue') +
  geom_line(data = temperature_series, aes(x = Date, y = TB_tempC_series), color = 'red')+
  scale_x_date(date_breaks = "year")

sal_ts = ts(salinity_series$TB_sal_series)
sal_ts = zoo::zoo(as.matrix(salinity_series[,-5:-7]), order.by = salinity_series$Date)

plot(salinity_series$TB_sal, salinity_series$MC_sal)
plot(salinity_series$TB_sal, salinity_series$sss)

plot(sal_ts)
# fit bsts model
model_components <- list()
# add a linear trend component
model_components <- AddAutoAr(model_components, y = salinity_series$TB_sal_series)

# add a seasonal component
model_components <- AddSeasonal(model_components, y = salinity_series$TB_sal_series, nseasons = 12)
summary(model_components)
# fit the model
fit <- bsts(model_components,salinity_series$TB_sal,  niter = 1e3)
summary(fit)
# view the model 
burnin <- 250

tibble(
  `Date` = as.Date(salinity_series$Date),
  trend = colMeans(fit$state.contributions[-(1:burnin), 'Ar1',]),
  seasonality = colMeans(fit$state.contributions[-(1:burnin),'seasonal.12.1',])) %>%
  pivot_longer(-Date, names_to = 'component', values_to = 'value') -> y

ggplot(y, aes(x = Date, y = value)) + geom_line() + facet_grid(component ~ ., scales = 'free')
