print("Have this script run whatever data cleaning you do")

# data manipulation

data_import <<- function(){
  temp_paths <- list.files(path = "./data/", "*temp_*", full.names = TRUE)
  temp_files <- lapply(temp_paths, read_csv, skip = 1,
                       col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))
  
  temp_df <<- bind_rows(temp_files) %>%
    rename(time = 'TS', temp_F = 'F') %>%
    mutate(yday = yday(time))
  
  do_paths <- list.files(path = "./data/", "*DO_*", full.names = TRUE)
  do_files <- lapply(do_paths, read_csv, skip = 1,
                     col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))
  
  do_df <<- bind_rows(do_files) %>%
    rename(time = 'TS', do_mg_L = 'mg L') %>%
    mutate(yday = yday(time))
  
  hydro_paths <- list.files(path = "./data/", "hydro{1}.*.csv", full.names= TRUE)
  hydro_files <- lapply(hydro_paths, read_csv, col_type = cols(ArrayID = col_skip(), Year = col_skip(), JDay= col_skip(),
                                                               X4 = col_skip(), X5 = col_skip(), X6 = col_skip(), Time = col_skip(), X8 = col_skip(), `Date/Time` = col_character(),
                                                               DO = col_skip(), Value4 = col_double(), Depth = col_skip(), `Depth feet` = col_skip()))
  hydro_df <<- bind_rows(hydro_files) %>%
    rename(time = 'Date/Time', temp_C = 'Temp', temp_F = 'Water Temp F', spCon_uS_cm = 'Value4', sal_psu = 'Sal') %>%
    mutate(time = parse_date_time(time, orders = c("%m/%d/%Y %H:%M", "%m/%d/%y %H:%M")),
           yday = yday(time))
  
  shallow<<- read_csv(file = "https://www.dropbox.com/s/bw1klp28wx5i7lx/TB2_Compiled_Data.csv?dl=1") %>% as.data.frame() %>%
    column_to_rownames('X1')
  deep<<- read_csv(file = "https://www.dropbox.com/s/581r05guoyj6aoc/Patch%20Mosaic%20Master2.csv?dl=1") %>% as.data.frame() %>%
    column_to_rownames('X1')
  site_list<<-read_csv(file ="https://www.dropbox.com/s/fmahl3dc34ss33x/terrebonne_site_master.csv?dl=1") %>% as.data.frame()
  TB_2018_core_meta<<-read_csv(file = "./data/metadata/TB_2018_Core-Metadata.csv", 
                               col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
                                                `Van Veen Grab` = "i", Core = "i"))
  TB_2019_core_meta<<- read_csv(file = "./data/metadata/TB_2019_Core-Metadata.csv",
                                col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
                                                 `Van Veen Grab` = "i", Core = "i"))
 
  TB_trawl_meta <- sheets_get("1mUG6oYBGGONBzgXPsvS6XZkNRUIojxtI7BV8WUTh-3I")#sheetID extracted with as_sheets_ID(*sheet URL*)
  TB_trawl_tabs <- nrow(TB_trawl_meta[['sheets']])
  TB_trawl_list <- vector('list', TB_trawl_tabs)
  for(tab in 1:TB_trawl_tabs){
    TB_trawl_list[[tab]] <- sheets_read(TB_trawl_meta[['spreadsheet_id']], sheet = tab, col_types = "T?d")
  }
  TB_trawl_data <<- bind_rows(TB_trawl_list)
  
  
}
data_import()

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

# working script for updating metadata and sample counts for each sampling date.

TB_2018_core_meta = TB_2018_core_meta %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))
TB_2019_core_meta = TB_2019_core_meta %>% rename(Date = "Date Taken") %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"))

TB_2018_core_summ = TB_2018_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())
TB_2019_core_summ = TB_2019_core_meta %>% group_by(Date, Station) %>% summarise(core_count = n())

TB_core_summ = bind_rows(TB_2018_core_summ, TB_2019_core_summ)

## working with trawl data
unique_spp = TB_trawl_data %>%
  mutate(Species = str_replace_all(Species, "(?i)-|mixed", ""),
         Species = str_replace(Species, "-", ","),
         Species2 = str_remove(Species, ".+;"),
         Species2 = ifelse(Species2 == Species, NA, Species2),
         Species = str_remove(Species,";.+"))
  unlist()%>%
  str_replace("(?i)-|Mixed", "") 
  strsplit(gsub("(.+);(.+)","\\1\1\\2", .),"\1")

