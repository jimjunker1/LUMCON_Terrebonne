source("install-packages.R")
# data manipulation

data_import <<- function(){
  temp_paths <- list.files(path = "./data/", "*temp_*", full.names = TRUE)
  temp_files <- lapply(temp_paths, read_csv, skip = 1,
                       col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))
  
  roundUP <<- function(x,to=10){
    to*(x%/%to + as.logical(x%%to))
  }
  
  tempF_to_C <<- function(F){(5/9)*(F-32)}
  
  temp_df <<- bind_rows(temp_files) %>%
    rename(datetime = 'TS', temp_F = 'F') %>%
    mutate(yday = yday(datetime))
  
  do_paths <- list.files(path = "./data/", "*DO_*", full.names = TRUE)
  do_files <- lapply(do_paths, read_csv, skip = 1,
                     col_types = cols(TS = col_datetime(format = "%m/%d/%Y %H:%M")))
  
  do_df <<- bind_rows(do_files) %>%
    rename(datetime = 'TS', do_mg_L = 'mg L') %>%
    mutate(yday = yday(datetime))
  
  hydro_paths <- list.files(path = "./data/", "hydro{1}.*.csv", full.names= TRUE)
  hydro_files <- lapply(hydro_paths, read_csv, col_type = cols(ArrayID = col_skip(), Year = col_skip(), JDay= col_skip(),
                                                               X4 = col_skip(), X5 = col_skip(), X6 = col_skip(), Time = col_skip(), X8 = col_skip(), `Date/Time` = col_character(),
                                                               DO = col_skip(), Value4 = col_double(), Depth = col_skip(), `Depth feet` = col_skip()))
  hydro_df <<- bind_rows(hydro_files) %>%
    rename(datetime = 'Date/Time', temp_C = 'Temp', temp_F = 'Water Temp F', spCon_uS_cm = 'Value4', sal_psu = 'Sal') %>%
    mutate(datetime = parse_date_time(datetime, orders = c("%m/%d/%Y %H:%M", "%m/%d/%y %H:%M")),
           yday = yday(datetime))
  # 2020-05-12 
  # these are throwing errors and not loading due to network timeout. May need to reset/update DNS servers if this continues.
  
  # shallow<<- read_csv(file = "https://www.dropbox.com/s/bw1klp28wx5i7lx/TB2_Compiled_Data.csv?dl=1") %>% as.data.frame() %>%
  #   column_to_rownames('X1')
  # deep<<- read_csv(file = "https://www.dropbox.com/s/581r05guoyj6aoc/Patch%20Mosaic%20Master2.csv?dl=1") %>% as.data.frame() %>%
  #   column_to_rownames('X1')
  # site_list<<-read_csv(file ="https://www.dropbox.com/s/fmahl3dc34ss33x/terrebonne_site_master.csv?dl=1") %>% as.data.frame()
  # TB_2018_core_meta<<-read_csv(file = "./data/metadata/TB_2018_Core-Metadata.csv", 
  #                              col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
  #                                               `Van Veen Grab` = "i", Core = "i"))
  # TB_2019_core_meta<<- read_csv(file = "./data/metadata/TB_2019_Core-Metadata.csv",
  #                               col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
  #                                                `Van Veen Grab` = "i", Core = "i"))
 #    # options(gargle_oauth_email = "jjunker@lumcon.edu")
 #    gs4_auth(email = "jjunker@lumcon.edu")
 #    # gs4_deauth()
 #    TB_trawl_meta <- gs4_get("1mUG6oYBGGONBzgXPsvS6XZkNRUIojxtI7BV8WUTh-3I")#sheetID extracted with as_sheets_ID(*sheet URL*)
 #    TB_trawl_tabs <- nrow(TB_trawl_meta[['sheets']])
 #    TB_trawl_list <- vector('list', TB_trawl_tabs)
 #    for(tab in 1:TB_trawl_tabs){
 #      TB_trawl_list[[tab]] <- sheets_read(TB_trawl_meta[['spreadsheet_id']], sheet = tab, col_types = "T?d", range = cell_cols("A:C"))
 #    }
 #    TB_trawl_data <- bind_rows(TB_trawl_list)
 # # 
 # ##### Data entry fixes ####
 #    TB_trawl_data <- TB_trawl_data %>%
 #     filter(!is.na(Species)) %>%
 #     mutate(Species = str_replace(Species,"spp(?!\\.)","spp."),
 #            species_mod = case_when(grepl("Litopanaeus|Litopenaeus|Farfantepenaeus|Farfantepanaeus|Farfantapenaeus", Species, ignore.case = TRUE) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
 #                                    grepl("Ariopsis felis|Ariopsisis|Ariosis|Bagre", Species) ~ "Ariopsis felis/Bagre marinus, Hardhead Catfish/Gafftopsail Catfish",
 #                                    grepl("Family Naticidae", Species) ~ "Family: Naticidae/Uropsalpinx cinerea, Moon Snail/Oyster drill snail",
 #                                    grepl("midshipmen", Species) ~ "Porichthys plectrodon, Atlantic Midshipman",
 #                                    grepl("shimp", Species) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
 #                                    grepl("Puffer", Species, ignore.case = TRUE) ~ "Sphoeroides parvus, Least Pufferfish",
 #                                    grepl("hermit crab", Species) ~ "Superfamily: Paguroidea, Hermit Crab",
 #                                    grepl("jellyish", Species, ignore.case = TRUE) ~ "Class: Scyphozoa",
 #                                    grepl("mantis", Species, ignore.case = TRUE) ~ "Order: Stomatopoda, Mantis Shrimp",
 #                                    grepl("Menticcirhus spp.", Species) ~ "Menticcirhus spp., Kingfish",
 #                                    grepl("oyster drill", Species, ignore.case = TRUE) ~ "Uropsalpinx cinerea, Oyster drill snail",
 #                                    grepl("Prionutus spp. Sea Robin", Species, ignore.case = TRUE) ~ "Prionutus spp., Sea Robin",
 #                                    grepl("silver kingfish", Species, ignore.case = TRUE) ~ "Megalops atlanticus, Atlantic Tarpon/Silver Kingfish",
 #                                    grepl("southern flounder", Species, ignore.case = TRUE) ~ "Paralichthys lethostigma, Southern Flounder",
 #                                    grepl("stonecrab", Species, ignore.case = TRUE) ~ "Menippe mercenaria, Stone Crab",
 #                                    grepl("Bristleworm|Polychae", Species) ~ "Class: Polychaeta",
 #                                    grepl("Isopod", Species) ~ "Order: Isopoda",
 #                                    grepl("Mojarra spp.", Species, ignore.case = TRUE) ~ "Family: Gerreidae, Mojarra spp.",
 #                                    grepl("Croaker*", Species, ignore.case = TRUE) ~ "Micropogonius undulatus, Atlantic Croaker",
 #                                    grepl("butter*", Species, ignore.case = TRUE) ~"Peprilus spp., American Butterfish",
 #                                    grepl("brief squid", Species, ignore.case = TRUE) ~ "Lolliguncula brevis, Atlantic Brief Squid",
 #                                    grepl("Anchovy|Ancovy|anchovie", Species, ignore.case = TRUE) ~ "Anchoa spp., Anchovy",
 #                                    grepl("robin", Species, ignore.case = TRUE) ~ "Prionotus spp., Sea Robin",
 #                                    grepl("ophiuroid*", Species, ignore.case = TRUE) ~ "Class: Ophiuroidae, Brittle star",
 #                                    grepl("F. Sciaenidae|F. Scianidae", Species, ignore.case = TRUE) ~ "F. Sciaenidae, Drum",
 #                                    grepl("Family: Monocanthidae|Monocanthidae", Species, ignore.case = TRUE) ~ "Family: Monacanthidae, Filefish",
 #                                    grepl("Clibanarius vittatus|Clibenarius vittatum", Species, ignore.case = TRUE) ~ "Clibanarius vittatus, Hermit Crab",
 #                                    grepl("Leiostomus xanthurus",Species, ignore.case = TRUE) ~ "Leiostomus xanthurus, Spot",
 #                                    grepl("Selene setapinnis", Species, ignore.case = TRUE) ~ "Selene setapinnis, Atlantic Moonfish",
 #                                    grepl("Chaetodipterus faber", Species, ignore.case = TRUE) ~ "Chaetodipterus faber, Atlantic Spadefish",
 #                                    grepl("Clibanarius spp.", Species, ignore.case = TRUE) ~ "Clibanarius spp., Hermit Crab",
 #                                    grepl("Larimus fasciatus", Species, ignore.case = TRUE) ~ "Larimus fasciatus, Banded Drum",
 #                                    TRUE ~ NA_character_),
 #            species_mod = if_else(is.na(species_mod), Species, species_mod),
 #            species_mod = str_replace(species_mod, "F. ", "Family: "),
 #            species_mod = str_replace(species_mod, "Family ","Family: "),
 #            species_mod = str_replace_all(species_mod, c("C. ", "O. "), c("Class: ","Order: ")),
 #            species_mod = str_replace_all(species_mod, c("Porychthys","Preprilus"), c("Porichthys","Peprilus")),
 #            species_mod = str_replace_all(species_mod, "Porselain", "Porcelain"),
 #            species_mod = str_replace(species_mod, "crabs","crab"),
 #            species_mod = str_replace(species_mod, "Mackeral","Mackerel"),
 #            Common_name = gsub('.*,\\s*(.*)','\\1',species_mod),
 #            Common_name = str_to_title(Common_name),
 #            # Common_name = str_replace_all(Common_name, c(".drum", "blue"), c(" Drum", "Blue")),
 #            # Common_name = str_replace_all(Common_name, c(".crab","flounder"),c(" Crab","Flounder")),
 #            species_mod = str_remove(species_mod, ",.+")) %>%
 #     mutate(year = year(Date),
 #            month = month(Date))
 #    # unique(levels(as.factor(TB_trawl_data$species_mod)))
 # 
 # ##### Specific scientific and common name fixes #####
 #    ## replace bad species names with good species names ##
 #   bad_names = list(c("Ariosis felis","Bagre marinus"),
 #                    c("bay anchovie"),
 #                    c("burfish", "Burfish"),
 #                    c("Chloroscombrus chysurus","Chloroscombrus chysurus","Chloroscomborus chysurus"),
 #                    c("Citharychthys spilopterus"),
 #                    c("Clibenarius vittatum"),
 #                    c("crevalle jack"),
 #                    c("Ctenophora"),
 #                    c("Cynoscin arenarius"),
 #                    c("Dasyatus sabina"),
 #                    c("Family: Rajiidae"),
 #                    c("Family: Scianidae"),
 #                    c("Farfantapenaeus aztecus"),
 #                    c("Menidia berylinia"),
 #                    c("Menticirrhus americanus"),
 #                    c("Menticirrhus saxatilis"),
 #                    c("Menticirrhus spp."),
 #                    c("Mojarra spp.,"),
 #                    c("Obsanus beta"),
 #                    c("Peprilus berti"),
 #                    c("Prionutus spp."),
 #                    c("Sciaenops ocellata"),
 #                    c("Stellifer lanceolata"),
 #                    c("Sygnathus scovelli"),
 #                    c("Urophysis floridanus")
 # 
 #   )
 #    good_names = list("Ariopsis felis/Bagre marinus",
 #                     "Anchoa mitchilli",
 #                     "Chilomycterus schoepfi",
 #                     "Chloroscombrus chrysurus",
 #                     "Citharichthys spilopterus",
 #                     "Clibanarius vittatus",
 #                     "Caranx hippos",
 #                     "Phylum: Ctenophora",
 #                     "Cynoscion arenarius",
 #                     "Dasyatis sabina",
 #                     "Family: Rajidae",
 #                     "Family: Sciaenidae",
 #                     "Litopenaeus setiferus/Farfantepenaeus aztecus",
 #                     "Menidia beryllina",
 #                     "Menticcirhus americanus",
 #                     "Menticcirhus saxatilis",
 #                     "Menticcirhus spp.",
 #                     "Mojarra spp.",
 #                     "Opsanus beta",
 #                     "Peprilus burti",
 #                     "Prionotus spp.",
 #                     "Sciaenops ocellatus",
 #                     "Stellifer lanceolatus",
 #                     "Syngnathus scovelli",
 #                     "Urophycis floridanus"
 #   )
 #    keyval <- setNames(rep(good_names, lengths(bad_names)), unlist(bad_names))
 #    # unique(levels(as.factor(TB_trawl_data$species_mod)))
 # 
 #    # coordinate common names for species either: w/ no common name in data or common name errors
 #   species_names = list(c("Class: Ophiuroidea"),
 #                        c("Class: Scyphozoa"),
 #                        c("Class: Polychaeta"),
 #                        c("Family: Brachyura", "Crab"),
 #                        c("Family: Blenniidae"),
 #                        c("Family: Gobiidae"),
 #                        c("Family: Malacostraca"),
 #                        c("Family: Ophiuroidae"),
 #                        c("Family: Rajidae"),
 #                        c("Family: Sciaenidae"),
 #                        c("Family: Monacanthidae"),
 #                        c("Family: Syngnathidae"),
 #                        c("Megalops atlanticus"),
 #                        c("Paralichthys lethostigma", "Paralichthys Lethostigma"),
 #                        c("Spotted Trout"),
 #                        c("Litopenaeus setiferus/Farfantepenaeus aztecus"),
 #                        c("Menippe mercenaria", "Menippe Mercenaria"),
 #                        c("Micropogonius undulatus"),
 #                        c("Order: Isopoda"),
 #                        c("Selene setapinnis"),
 #                        c("Stellifer lanceolatus","Stellifer Lanceolatus"),
 #                        c("Syngnathus scovelli"),
 #                        c("Larimus fasciatus"),
 #                        c("Anchoa mitchilli"),
 #                        c("Baywhiff"),
 #                        c("atlantic bumper","Chloroscombrus chrysurus", "Chloroscombrus Chrysurus"),
 #                        c("burfish"),
 #                        c("crevalle jack"),
 #                        c("Hogchoaker"),
 #                        c("Sciaenops ocellatus", "Sciaenops Ocellatus"),
 #                        c("Carcharhinus limbatus","Carcharhinus Limbatus")
 #   )
 #    common_names = list("Brittle star",#Class: Ophiuroidea
 #                       "Jellyfish",#Class: Scyphozoa
 #                       "Polychaetes",#Class: Polychaeta
 #                       "Crab spp.",#Family: Brachyura
 #                       "Blennies",#Family: Blenniidae
 #                       "Goby",#Family: Gobiidae"
 #                       "Crustacean indet.",#Family: Malacostraca
 #                       "Brittle star",#Family: Ophiuroidae
 #                       "Skate",#Family: Rajidae
 #                       "Drum",#"Family: Sciaenidae"
 #                       "Filefish",#Family: Monacanthidae
 #                       "Pipefish",#Family: Syngnathidae
 #                       "Silver kingfish",#Megalops atlanticus
 #                       "Southern flounder",#Paralichthys lethostigma
 #                       "Spotted Seatrout",#Spotted Trout
 #                       "White Shrimp/Brown Shrimp",
 #                       "Florida Stone Crab",
 #                       "Atlantic Croaker",
 #                       "Isopod",
 #                       "Atlantic Moonfish",
 #                       "Star Drum",
 #                       "Gulf Pipefish",
 #                       "Banded Drum",
 #                       "Bay Anchovy",
 #                       "Bay Whiff",
 #                       "Atlantic Bumper",
 #                       "Burrfish",
 #                       "Crevalle Jack",
 #                       "Hogchoker",
 #                       "Red Drum",
 #                       "Blacktip Shark"
 #   )
 #    commonkey <- setNames(rep(common_names, lengths(species_names)), unlist(species_names))
 #    TB_trawl_data <<- TB_trawl_data %>%
 #     mutate(species_mod = recode(species_mod, !!!keyval),
 #            Common_name = recode(Common_name, !!!commonkey),
 #             date_id = paste0(as.character(year),"-",as.character(month)))
 #   
#### End scientific and common name fixes #### 
   #### End Data entry fixes ####
  # saveRDS(TB_trawl_data, file = "./data/TB_trawl_data.rds")
  TB_trawl_data <<- readRDS(file = "./data/TB_trawl_data.rds")
  # load(ocecolors)
  
  # working with bird data
 # ebird_data <<- read_ebd("./data/ebd_US-LA_relFeb-2020/ebd_US-LA_relFeb-2020.txt") %>%
 #   auk_bbox(bbox = c(-90.727, 29.122, -90.58, 29.43))
  }
data_import()
print("Have this script run whatever data cleaning you do")