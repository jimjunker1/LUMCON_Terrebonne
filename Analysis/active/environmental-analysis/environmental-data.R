## source in environmental data ##

source("datascript.R")

## combine the TB and MC datasets to fill in gaps
full_dates = seq.Date(from = as.Date("2008-01-01"), to = Sys.Date(), by = 1) %>% data.frame %>% setNames('Date')

##temperature
temp_combined = left_join(MC_temperature_df %>% rename(MC_temp = 'temp_C'), TB_temp_df %>% rename(TB_temp = 'temp_C'), by = 'Date') %>%
  right_join(full_dates, by = 'Date', all = TRUE) %>% ungroup() %>%
  join(regional_temp_df)
#bivariate plot of TB and MC temperature series
ggplot(temp_combined, aes(x = MC_temp, y = TB_temp)) + geom_point(size = 1.2, alpha = 0.5)
ggplot(temp_combined, aes(x = sst, y = TB_temp)) + geom_point(size = 1.2, alpha = 0.5)

ggplot(temp_combined, aes(x = Date)) + geom_line(aes(y = MC_temp), colour = 'red') +
  geom_line(aes(y = TB_temp), colour = 'blue') +
  geom_line(aes(y = sst), colour = 'black')

# build a forecasted model with Marine center and regional SST to fill NAs
library(forecast)
TB_temp_ts = ts(temp_combined$TB_temp)
set.seed(42);ar.model = auto.arima(temp_combined$TB_temp, xreg = as.matrix(temp_combined[,c(2,4)]))
ar.model

TB_nas = is.na(temp_combined$TB_temp)
sum(TB_nas)
set.seed(123);TB_forecast = forecast(ar.model, xreg = as.matrix(temp_combined[,c(2,4)]))

temperature_series = temp_combined %>% 
  bind_cols(TB_forecast$mean) %>% 
  rename(forecast = '...5') %>%
  mutate(TB_tempC_series = ifelse(is.na(TB_temp), forecast,TB_temp),
         TB_tempC_series = ifelse(is.na(TB_tempC_series), sst, TB_tempC_series))
    
ggplot(temperature_series, aes( x = Date)) +
  geom_line(aes(y = TB_tempC_series), colour = 'blue') +
  geom_line(aes(y = sst), colour = 'red')+
  geom_line(aes(y = TB_temp), colour = 'black')

## salinity

salinity_combined = left_join(MC_salinity_df %>% rename(MC_sal = 'sal_psu'), TB_salinity_df %>% rename(TB_sal = 'sal_psu')) %>%
  right_join(full_dates,by = 'Date', all = TRUE)

ggplot(salinity_combined, aes(x = MC_sal, y = TB_sal)) + geom_point(size = 1.2, alpha = 0.5) 

ggplot(salinity_combined, aes(Date)) +
  geom_line(aes(y = MC_sal), colour = 'blue')+
  geom_line(aes(y = TB_sal), colour = 'black')
