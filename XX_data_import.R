# data manipulation

data_import <<- function(){
temp_paths <- list.files(path = "./data/", "*temp_*", full.names = TRUE)
temp_files <- lapply(temp_paths, read_csv, skip = 1,
                     col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))

temp_df <<- bind_rows(temp_files)

do_paths <- list.files(path = "./data/", "*DO_*", full.names = TRUE)
do_files <- lapply(do_paths, read_csv, skip = 1,
                   col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))

do_df <<- bind_rows(do_files)

hydro_paths <- list.files(path = "./data/", "hydro{1}.*.csv", full.names= TRUE)
hydro_files <- lapply(hydro_paths, read_csv, col_type = cols(ArrayID = col_skip(), Year = col_skip(), JDay= col_skip(),
X4 = col_skip(), X5 = col_skip(), X6 = col_skip(), Time = col_skip(), X8 = col_skip(), `Date/Time` = col_character(),
Value4 = col_skip(), Depth = col_skip(), `Depth feet` = col_skip()))
hydro_df <<- bind_rows(hydro_files) %>%
  rename(time = 'Date/Time', temp_C = 'Temp', temp_F = 'Water Temp F', sal_psu = 'Sal', do_mg_L = 'DO') %>%
  mutate(time = parse_date_time(time, orders = c("%m/%d/%Y %H:%M", "%m/%d/%y %H:%M")))
}
data_import()
