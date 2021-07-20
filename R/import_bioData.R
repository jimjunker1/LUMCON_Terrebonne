import_bioData <- function(){
   source(here::here("R/clean_trawl_spp.R"))
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
   # options(gargle_oauth_email = "jjunker@lumcon.edu")
   gs4_auth(email = "jjunker@lumcon.edu")
   # gs4_deauth()
   TB_trawl_meta <- gs4_get("1mUG6oYBGGONBzgXPsvS6XZkNRUIojxtI7BV8WUTh-3I")#sheetID extracted with as_sheets_ID(*sheet URL*)
   TB_trawl_tabs <- nrow(TB_trawl_meta[['sheets']])
   TB_trawl_list <- vector('list', TB_trawl_tabs)
   for(tab in 1:TB_trawl_tabs){
     TB_trawl_list[[tab]] <- sheets_read(TB_trawl_meta[['spreadsheet_id']], sheet = tab, col_types = "T?d", range = cell_cols("A:C"))
   }
   TB_trawl_data <- bind_rows(TB_trawl_list) %>% dplyr::select(Date, Species, Abundance)
#
##### Data entry fixes ####
   TB_trawl_data <- TB_trawl_data %>%
    filter(!is.na(Species)) %>%
    mutate(Species = str_replace(Species,"spp(?!\\.)","spp."),
           species_mod = case_when(grepl("Litopanaeus|Litopenaeus|Farfantepenaeus|Farfantepanaeus|Farfantapenaeus", Species, ignore.case = TRUE) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
                                   grepl("Ariopsis felis|Ariopsisis|Ariosis|Bagre", Species) ~ "Ariopsis felis/Bagre marinus, Hardhead Catfish/Gafftopsail Catfish",
                                   grepl("Family Naticidae", Species) ~ "Family: Naticidae/Uropsalpinx cinerea, Moon Snail/Oyster drill snail",
                                   grepl("midshipmen", Species) ~ "Porichthys plectrodon, Atlantic Midshipman",
                                   grepl("shimp", Species) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
                                   grepl("Puffer", Species, ignore.case = TRUE) ~ "Sphoeroides parvus, Least Pufferfish",
                                   grepl("hermit crab", Species) ~ "Superfamily: Paguroidea, Hermit Crab",
                                   grepl("jellyish", Species, ignore.case = TRUE) ~ "Class: Scyphozoa",
                                   grepl("mantis", Species, ignore.case = TRUE) ~ "Order: Stomatopoda, Mantis Shrimp",
                                   grepl("Menticcirhus spp.", Species) ~ "Menticcirhus spp., Kingfish",
                                   grepl("oyster drill", Species, ignore.case = TRUE) ~ "Uropsalpinx cinerea, Oyster drill snail",
                                   grepl("Prionutus spp. Sea Robin", Species, ignore.case = TRUE) ~ "Prionutus spp., Sea Robin",
                                   grepl("silver kingfish", Species, ignore.case = TRUE) ~ "Megalops atlanticus, Atlantic Tarpon/Silver Kingfish",
                                   grepl("southern flounder", Species, ignore.case = TRUE) ~ "Paralichthys lethostigma, Southern Flounder",
                                   grepl("stonecrab", Species, ignore.case = TRUE) ~ "Menippe mercenaria, Stone Crab",
                                   grepl("Bristleworm|Polychae", Species) ~ "Class: Polychaeta",
                                   grepl("Isopod", Species) ~ "Order: Isopoda",
                                   grepl("Mojarra spp.", Species, ignore.case = TRUE) ~ "Family: Gerreidae, Mojarra spp.",
                                   grepl("Croaker*", Species, ignore.case = TRUE) ~ "Micropogonius undulatus, Atlantic Croaker",
                                   grepl("butter*", Species, ignore.case = TRUE) ~"Peprilus spp., American Butterfish",
                                   grepl("brief squid", Species, ignore.case = TRUE) ~ "Lolliguncula brevis, Atlantic Brief Squid",
                                   grepl("Anchovy|Ancovy|anchovie", Species, ignore.case = TRUE) ~ "Anchoa spp., Anchovy",
                                   grepl("robin", Species, ignore.case = TRUE) ~ "Prionotus spp., Sea Robin",
                                   grepl("ophiuroid*", Species, ignore.case = TRUE) ~ "Class: Ophiuroidae, Brittle star",
                                   grepl("F. Sciaenidae|F. Scianidae", Species, ignore.case = TRUE) ~ "F. Sciaenidae, Drum",
                                   grepl("Family: Monocanthidae|Monocanthidae", Species, ignore.case = TRUE) ~ "Family: Monacanthidae, Filefish",
                                   grepl("Clibanarius vittatus|Clibenarius vittatum", Species, ignore.case = TRUE) ~ "Clibanarius vittatus, Hermit Crab",
                                   grepl("Leiostomus xanthurus",Species, ignore.case = TRUE) ~ "Leiostomus xanthurus, Spot",
                                   grepl("Selene setapinnis", Species, ignore.case = TRUE) ~ "Selene setapinnis, Atlantic Moonfish",
                                   grepl("Chaetodipterus faber", Species, ignore.case = TRUE) ~ "Chaetodipterus faber, Atlantic Spadefish",
                                   grepl("Clibanarius spp.", Species, ignore.case = TRUE) ~ "Clibanarius spp., Hermit Crab",
                                   grepl("Larimus fasciatus", Species, ignore.case = TRUE) ~ "Larimus fasciatus, Banded Drum",
                                   TRUE ~ NA_character_),
           species_mod = if_else(is.na(species_mod), Species, species_mod),
           species_mod = str_replace(species_mod, "F. ", "Family: "),
           species_mod = str_replace(species_mod, "Family ","Family: "),
           species_mod = str_replace_all(species_mod, c("C. ", "O. "), c("Class: ","Order: ")),
           species_mod = str_replace_all(species_mod, c("Porychthys","Preprilus"), c("Porichthys","Peprilus")),
           species_mod = str_replace_all(species_mod, "Porselain", "Porcelain"),
           species_mod = str_replace(species_mod, "crabs","crab"),
           species_mod = str_replace(species_mod, "Mackeral","Mackerel"),
           Common_name = gsub('.*,\\s*(.*)','\\1',species_mod),
           Common_name = str_to_title(Common_name),
           # Common_name = str_replace_all(Common_name, c(".drum", "blue"), c(" Drum", "Blue")),
           # Common_name = str_replace_all(Common_name, c(".crab","flounder"),c(" Crab","Flounder")),
           species_mod = str_remove(species_mod, ",.+")) %>%
    mutate(year = year(Date),
           month = month(Date))
   # unique(levels(as.factor(TB_trawl_data$species_mod)))

##### Specific scientific and common name fixes #####
   ## replace bad species names with good species names ##
  bad_names = list(c("Ariosis felis","Bagre marinus"),
                   c("bay anchovie"),
                   c("burfish", "Burfish"),
                   c("Chloroscombrus chysurus","Chloroscombrus chysurus","Chloroscomborus chysurus"),
                   c("Citharychthys spilopterus"),
                   c("Clibenarius vittatum"),
                   c("crevalle jack"),
                   c("Ctenophora"),
                   c("Cynoscin arenarius"),
                   c("Dasyatus sabina"),
                   c("Family: Rajiidae"),
                   c("Family: Scianidae"),
                   c("Farfantapenaeus aztecus"),
                   c("Menidia berylinia"),
                   c("Menticirrhus americanus"),
                   c("Menticirrhus saxatilis"),
                   c("Menticirrhus spp."),
                   c("Mojarra spp.,"),
                   c("Obsanus beta"),
                   c("Peprilus berti"),
                   c("Prionutus spp."),
                   c("Sciaenops ocellata"),
                   c("Stellifer lanceolata"),
                   c("Sygnathus scovelli"),
                   c("Urophysis floridanus")

  )
   good_names = list("Ariopsis felis/Bagre marinus",
                    "Anchoa mitchilli",
                    "Chilomycterus schoepfi",
                    "Chloroscombrus chrysurus",
                    "Citharichthys spilopterus",
                    "Clibanarius vittatus",
                    "Caranx hippos",
                    "Phylum: Ctenophora",
                    "Cynoscion arenarius",
                    "Dasyatis sabina",
                    "Family: Rajidae",
                    "Family: Sciaenidae",
                    "Litopenaeus setiferus/Farfantepenaeus aztecus",
                    "Menidia beryllina",
                    "Menticcirhus americanus",
                    "Menticcirhus saxatilis",
                    "Menticcirhus spp.",
                    "Mojarra spp.",
                    "Opsanus beta",
                    "Peprilus burti",
                    "Prionotus spp.",
                    "Sciaenops ocellatus",
                    "Stellifer lanceolatus",
                    "Syngnathus scovelli",
                    "Urophycis floridanus"
  )
   keyval <- setNames(rep(good_names, lengths(bad_names)), unlist(bad_names))
   # unique(levels(as.factor(TB_trawl_data$species_mod)))

   # coordinate common names for species either: w/ no common name in data or common name errors
  species_names = list(c("Class: Ophiuroidea"),
                       c("Class: Scyphozoa"),
                       c("Class: Polychaeta"),
                       c("Family: Brachyura", "Crab"),
                       c("Family: Blenniidae"),
                       c("Family: Gobiidae"),
                       c("Family: Malacostraca"),
                       c("Family: Ophiuroidae"),
                       c("Family: Rajidae"),
                       c("Family: Sciaenidae"),
                       c("Family: Monacanthidae"),
                       c("Family: Syngnathidae"),
                       c("Megalops atlanticus"),
                       c("Paralichthys lethostigma", "Paralichthys Lethostigma"),
                       c("Spotted Trout"),
                       c("Litopenaeus setiferus/Farfantepenaeus aztecus"),
                       c("Menippe mercenaria", "Menippe Mercenaria"),
                       c("Micropogonius undulatus"),
                       c("Order: Isopoda"),
                       c("Selene setapinnis"),
                       c("Stellifer lanceolatus","Stellifer Lanceolatus"),
                       c("Syngnathus scovelli"),
                       c("Larimus fasciatus"),
                       c("Anchoa mitchilli"),
                       c("Baywhiff"),
                       c("atlantic bumper","Chloroscombrus chrysurus", "Chloroscombrus Chrysurus"),
                       c("burfish"),
                       c("crevalle jack"),
                       c("Hogchoaker"),
                       c("Sciaenops ocellatus", "Sciaenops Ocellatus"),
                       c("Carcharhinus limbatus","Carcharhinus Limbatus")
  )
   common_names = list("Brittle star",#Class: Ophiuroidea
                      "Jellyfish",#Class: Scyphozoa
                      "Polychaetes",#Class: Polychaeta
                      "Crab spp.",#Family: Brachyura
                      "Blennies",#Family: Blenniidae
                      "Goby",#Family: Gobiidae"
                      "Crustacean indet.",#Family: Malacostraca
                      "Brittle star",#Family: Ophiuroidae
                      "Skate",#Family: Rajidae
                      "Drum",#"Family: Sciaenidae"
                      "Filefish",#Family: Monacanthidae
                      "Pipefish",#Family: Syngnathidae
                      "Silver kingfish",#Megalops atlanticus
                      "Southern flounder",#Paralichthys lethostigma
                      "Spotted Seatrout",#Spotted Trout
                      "White Shrimp/Brown Shrimp",
                      "Florida Stone Crab",
                      "Atlantic Croaker",
                      "Isopod",
                      "Atlantic Moonfish",
                      "Star Drum",
                      "Gulf Pipefish",
                      "Banded Drum",
                      "Bay Anchovy",
                      "Bay Whiff",
                      "Atlantic Bumper",
                      "Burrfish",
                      "Crevalle Jack",
                      "Hogchoker",
                      "Red Drum",
                      "Blacktip Shark"
  )
   commonkey <- setNames(rep(common_names, lengths(species_names)), unlist(species_names))
   TB_trawl_data <<- TB_trawl_data %>%
    mutate(species_mod = recode(species_mod, !!!keyval),
           Common_name = recode(Common_name, !!!commonkey),
            date_id = paste0(as.character(year),"-",as.character(month)))

#### End scientific and common name fixes #### 
#### Confirm Species and Common names from taxonomic databases
   debugonce(clean_trawl_spp)
   new_spp =  TB_trawl_data %>% clean_trawl_spp(spp_col = "species_mod", common_col = "Common_name",
                                     taxonomy_list = here::here("data/taxonomy_valid_common.rds")) %>% unique
   
#### End Data entry fixes ####
 saveRDS(TB_trawl_data, file = "./data/TB_trawl_data.rds")
 biology_sysDate = Sys.Date()
 save(biology_sysDate, file = "./data/biology_sysDate.rda")
}