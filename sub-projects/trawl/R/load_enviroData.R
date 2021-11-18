here::i_am('sub-projects/trawl/R/load_enviroData.R')
## source in environmental data ##

load_environmental_data <- function(){
temp_paths <- list.files(path = here::here("data/environment"), "Terrebonne.*temp_.*2020.*", full.names = TRUE)
temp_files <- lapply(temp_paths, read_csv,
                     col_types = cols('Date/Time' = col_datetime(format = "%m/%d/%Y %H:%M")))

roundUP <<- function(x,to=10){
  to*(x%/%to + as.logical(x%%to))
}

tempF_to_C <<- function(F){(5/9)*(F-32)}

TB_temp_df <<- bind_rows(temp_files) %>%
  dplyr::select(datetime = 'Date/Time', temp_C = 'Water Temp C') %>%
  dplyr::mutate(Date = as.Date(datetime)) %>%
  dplyr::select(-datetime) %>%
  group_by(Date) %>%
  summarise_all(mean, na.rm = TRUE)

# do_paths <- list.files(path = "./data/", "*DO_*", full.names = TRUE)
# do_files <- lapply(do_paths, read_csv, skip = 1,
#                    col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))
# 
# TB_do_df <<- bind_rows(do_files) %>%
#   rename(datetime = 'TS', do_mg_L = 'mg L') %>%
#   dplyr::mutate(yday = yday(datetime))
# 
# hydro_paths <- list.files(path = "./data/", "hydro{1}.*.csv", full.names= TRUE)
# hydro_files <- lapply(hydro_paths, read_csv, col_type = cols(ArrayID = col_skip(), Year = col_skip(), JDay= col_skip(),
#                                                              X4 = col_skip(), X5 = col_skip(), X6 = col_skip(), Time = col_skip(), X8 = col_skip(), `Date/Time` = col_character(),
#                                                              DO = col_skip(), Value4 = col_double(), Depth = col_skip(), `Depth feet` = col_skip()))
# TB_hydro_df <<- bind_rows(hydro_files) %>%
#   rename(datetime = 'Date/Time', temp_C = 'Temp', temp_F = 'Water Temp F', spCon_uS_cm = 'Value4', sal_psu = 'Sal') %>%
#   dplyr::mutate(datetime = parse_date_time(datetime, orders = c("%m/%d/%Y %H:%M", "%m/%d/%y %H:%M")),
#          yday = yday(datetime))

sal_paths = list.files(path = here::here("data/environment"), "Terrebonne.*salinity.*2020.*", full.names = TRUE)
TB_salinity_df <<- lapply(sal_paths, read_csv,
                          col_types = cols('Date/Time' = col_datetime(format = "%m/%d/%Y %H:%M"))) %>%
  bind_rows %>%
  dplyr::select(datetime = 'Date/Time', sal_psu = contains('salinity')) %>%
  dplyr::mutate(Date = as.Date(datetime)) %>%
  dplyr::select(-datetime) %>%
  group_by(Date) %>%
  summarise_all(mean, na.rm = TRUE)

marine_center_paths = list.files(path = here::here("data/environment"), "MarineCenter{1}.*.csv", full.names = TRUE)
marine_center_files = lapply(marine_center_paths, read_csv, col_type = cols(`Date/Time` = col_character())) %>%
  map(~.x %>% dplyr::rename(datetime = 'Date/Time') %>%
        dplyr::mutate(datetime = parse_date_time(datetime, orders = c("%m/%d/%Y %H:%M", "%m/%d/%y %H:%M")),
               Date = as.Date(datetime)) %>%
        dplyr::select(-datetime) %>%
        group_by(Date) %>%
        summarise_all(mean, na.rm = TRUE))

MC_salinity_df <<- marine_center_files[[1]] %>% dplyr::select(Date, sal_psu = "Sal")
MC_temperature_df <<- marine_center_files[[2]] %>% dplyr::select(Date, temp_C = 'Water Temp C')
MC_depth_df <<- marine_center_files[[3]] %>% dplyr::select(Date, depth_m = 'Water Depth M')
## combine the TB, MC, and regional datasets to fill in gaps

regional_temp_paths = list.files(path = here::here("data/environment"), "^X.*.nc", full.names = TRUE)
regional_temp_list = ncdf4::nc_open(regional_temp_paths)
regional_time = ncvar_get(regional_temp_list, "time")
regional_time_units = ncatt_get(regional_temp_list, 'time', 'units')
temp_array = ncvar_get(regional_temp_list, 'sst')
regional_time_date = as.Date(regional_time, origin =  unlist(strsplit(regional_time_units$value," "))[3], format = "%Y-%m-%d")
temp_avg = apply(temp_array,3,mean, na.rm = TRUE)
regional_temp_df <<- data.frame(Date = as.Date(regional_time_date), sst = temp_avg) %>% as_tibble

full_dates = seq.Date(from = as.Date("2009-01-01"), to = as.Date("2019-12-31"), by = 1) %>% data.frame %>% setNames('Date')

##temperature
temp_combined = left_join(MC_temperature_df %>% rename(MC_temp = 'temp_C'), TB_temp_df %>% rename(TB_temp = 'temp_C'), by = 'Date') %>%
  right_join(full_dates, by = 'Date', all = TRUE) %>% ungroup() %>%
  join(regional_temp_df)
set.seed(42);ar.model = auto.arima(temp_combined$TB_temp, xreg = as.matrix(temp_combined[,c(2,4)]))
set.seed(123);TB_forecast = forecast(ar.model, xreg = as.matrix(temp_combined[,c(2,4)]))

temperature_series <<- temp_combined %>% 
  bind_cols(TB_forecast$mean) %>% 
  rename(forecast = '...5') %>%
  dplyr::mutate(TB_tempC_series = ifelse(is.na(TB_temp), forecast,TB_temp),
         TB_tempC_series = ifelse(is.na(TB_tempC_series), sst, TB_tempC_series))
## salinity
regional_salinity_paths = list.files(path = here::here("data/environment"),"^aggregate.*.nc", full.names = TRUE)
regional_salinity_list = ncdf4::nc_open(regional_salinity_paths[[1]])
regional_salinity_time = ncdf4::ncvar_get(regional_salinity_list, 'time')
regional_salinity_time_units = ncdf4::ncatt_get(regional_salinity_list,'time', 'units')
regional_salinity_time_date = as.Date(regional_salinity_time, origin = unlist(strsplit(regional_salinity_time_units$value," "))[3], format = "%Y-%m-%d")
regional_salinity_array = ncdf4::ncvar_get(regional_salinity_list, 'l3m_data')
# regional_salinity_array[[2]] = ncdf4::ncvar_get(regional_salinity_list[[2]], 'sss_smap_40km')
regional_salinity_array = apply(regional_salinity_array,2,mean, na.rm = TRUE)
regional_sss_df <<- data.frame(Date = as.Date(regional_salinity_time_date), sss = regional_salinity_array)
salinity_combined = left_join(MC_salinity_df %>% rename(MC_sal = 'sal_psu'), TB_salinity_df %>% rename(TB_sal = 'sal_psu')) %>%
  right_join(full_dates,by = 'Date', all = TRUE) %>% join(regional_sss_df)
set.seed(42);sal.ar.model1 = auto.arima(salinity_combined$TB_sal, xreg = as.matrix(salinity_combined[,2]))
sal.ar.model1
set.seed(42);sal.ar.model2 = auto.arima(salinity_combined$TB_sal, xreg = as.matrix(salinity_combined[,c(2,4)]))
sal.ar.model2
TB_sal_forecast1 = forecast(sal.ar.model1, xreg = as.matrix(salinity_combined[,2]))
TB_sal_forecast2 = forecast(sal.ar.model2, xreg = as.matrix(salinity_combined[,c(2,4)]))

salinity_series <<- salinity_combined %>%
  bind_cols(forecast1 = TB_sal_forecast1$mean) %>%
  bind_cols(forecast2 =TB_sal_forecast2$mean) %>%
  dplyr::mutate(TB_sal_series = ifelse(is.na(TB_sal), forecast1, TB_sal),
         TB_sal_series = ifelse(is.na(TB_sal_series), forecast2, TB_sal_series))

 }
 load_environmental_data()
 rm(load_environmental_data)
