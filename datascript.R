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
  sheets_auth(email = "jjunker@lumcon.edu")
  TB_trawl_meta <- sheets_get("1mUG6oYBGGONBzgXPsvS6XZkNRUIojxtI7BV8WUTh-3I")#sheetID extracted with as_sheets_ID(*sheet URL*)
  TB_trawl_tabs <- nrow(TB_trawl_meta[['sheets']])
  TB_trawl_list <- vector('list', TB_trawl_tabs)
  for(tab in 1:TB_trawl_tabs){
    TB_trawl_list[[tab]] <- sheets_read(TB_trawl_meta[['spreadsheet_id']], sheet = tab, col_types = "T?d", range = cell_cols("A:C"))
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

# unique_spp  = TB_trawl_data %>%
#   pull(Species) %>%
#   as.factor() %>%
#   unique()
# unique_spp

TB_trawl_data = TB_trawl_data %>%
  filter(!is.na(Species)) %>%
  mutate(Species = str_replace(Species,"spp(?!\\.)","spp."),
    species_mod = case_when(grepl("Litopanaeus|Litopenaeus|Farfantepenaeus|Farfantepanaeus", Species) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
                                 grepl("Ariopsis felis|Ariopsisis|Ariosis", Species) ~ "Ariopsis felis/Bagre marinus, Hardhead Catfish/Gafftopsail Catfish",
                                 grepl("Family Naticidae", Species) ~ "Family: Naticidae/Uropsalpinx cinerea, Moon Snail/Oyster drill snail",
                                 grepl("midshipmen", Species) ~ "Porichthys plectrodon",
                                 grepl("shimp", Species) ~ "Litopenaeus setiferus/Farfantepenaeus aztecus,White Shrimp/Brown Shrimp",
                                 grepl("Puffer", Species, ignore.case = TRUE) ~ "Sphoeroides parvus, Least Pufferfish",
                                 grepl("hermit crab", Species) ~ "Superfamily: Paguroidea, Hermit crabs",
                                 grepl("jellyish", Species, ignore.case = TRUE) ~ "Class: Scyphozoa",
                                 grepl("mantis", Species, ignore.case = TRUE) ~ "Order: Stomatopoda, Mantis Shrimp",
                                 grepl("Menticcirhus spp.", Species) ~ "Menticcirhus spp., Kingfish",
                                 grepl("oyster drill", Species, ignore.case = TRUE) ~ "Uropsalpinx cinerea, Oyster drill snail",
                                 grepl("Prionutus spp. Sea Robin", Species, ignore.case = TRUE) ~ "Prionutus spp., Sea Robin",
                                 grepl("silver kingfish", Species, ignore.case = TRUE) ~ "Megalops atlanticus, Atlantic Tarpon/Silver Kingfish",
                                 grepl("southern flounder", Species, ignore.case = TRUE) ~ "Paralichthys lethostigma, Southern Flounder",
                                 grepl("stonecrab", Species, ignore.case = TRUE) ~ "Menippe mercenaria, Stonecrab",
                                 grepl("Bristleworm|Polychae", Species) ~ "Class: Polychaeta",
                                 grepl("Isopod", Species) ~ "Order: Isopoda",
                                 grepl("Mojarra spp.", Species, ignore.case = TRUE) ~ "Family: Gerreidae, Mojarra spp.",
                                 TRUE ~ NA_character_),
         species_mod = if_else(is.na(species_mod), Species, species_mod),
         species_mod = str_replace(species_mod, "F. ", "Family: "),
         species_mod = str_replace(species_mod, "Family ","Family: "),
         species_mod = str_replace_all(species_mod, c("C. ", "O. "), c("Class: ","Order: ")),
         species_mod = str_replace_all(species_mod, c("Porychthys","Preprilus"), c("Porichthys","Peprilus")),
         Common_name = gsub('.*,\\s*(.*)','\\1',species_mod),
         species_mod = str_remove(species_mod, ",.+")) %>%
  mutate(year = year(Date),
         month = month(Date))

unique_spp  = TB_trawl_data %>%
  pull(species_mod) %>%
  as.factor() %>%
  unique()

unique(levels(as.factor(TB_trawl_data$species_mod)))
## replace bad names with good names ##
bad_names = list(c("Ariosis felis","Bagre marinus"),
                 c("bay anchovie"),
                 c("burfish"),
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
                 c("Urophysis floridanus"),
                 c("Porichthys plectrodon")
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
                  "Urophycis floridanus",
                  "Atlantic Midshipman"
                  )

keyval <- setNames(rep(good_names, lengths(bad_names)), unlist(bad_names))

unique(levels(as.factor(TB_trawl_data$species_mod)))

species_names = list(c("Class: Ophiuroidea"),
                     c("Class: Scyphozoa"),
                     c("Class: Polychaeta"),
                     c("Family: Brachyura"),
                     c("Family: Blenniidae"),
                     c("Family: Gobiidae"),
                     c("Family: Malacostraca"),
                     c("Family: Ophiuroidae"),
                     c("Family: Rajidae"),
                     c("Family: Sciaenidae"),
                     c("Megalops atlanticus"),
                     c("Paralichthys lethostigma"),
                     c("Spotted Trout"),
                     c("Litopenaeus setiferus/Farfantepenaeus aztecus"),
                     c("Menippe mercenaria"),
                     c("Micropogonius undulatus"),
                     c("Order: Isopoda"),
                     c("Selene setapinnis"),
                     c("Stellifer lanceolatus"),
                     c("Syngnathus scovelli")
                     )
common_names = list("Brittle star",#Class: Ophiuroidea
                    "Jellyfish",#Class: Scyphozoa
                    "Polychaetes",#Class: Polychaeta
                    "Crab",#Family: Brachyura
                    "Blennies",#Family: Blenniidae
                    "Goby",#Family: Gobiidae"
                    "Crustacean indet.",#Family: Malacostraca
                    "Brittle star",#Family: Ophiuroidae
                    "Skate",#Family: Rajidae
                    "Drum",#"Family: Sciaenidae"
                    "Silver kingfish",#Megalops atlanticus
                    "Southern flounder",#Paralichthys lethostigma
                    "Spotted Seatrout",#Spotted Trout
                    "White Shrimp/Brown Shrimp",
                    "Florida Stone Crab",
                    "Atlantic Croaker",
                    "Isopod",
                    "Atlantic Moonfish",
                    "Star Drum",
                    "Gulf Pipefish"
                    )

commonkey <- setNames(rep(common_names, lengths(species_names)), unlist(species_names))

TB_trawl_data = TB_trawl_data %>%
  mutate(species_mod = recode(species_mod, !!!keyval),
         Common_name = recode(Common_name, !!!commonkey),
         date_id = paste0(as.character(year),"-",as.character(month)))

TB_trawl_taxasite<- TB_trawl_data %>%
  select(year, month, date_id, species_mod) %>%
  group_by(date_id, species_mod) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = species_mod, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

set.seed(123)
TB_NMDS <- vegan::metaMDS(TB_trawl_taxasite %>% select(-c(1:2)) %>% column_to_rownames("date_id"), distance = "jaccard", trymax = 1000, autotransform = FALSE)
TB_NMDS
plot(TB_NMDS)

TB_NMDS.scrs <- as.data.frame(scores(TB_NMDS, display = 'sites'))
TB_NMDS.scrs <- TB_NMDS.scrs %>% bind_cols(TB_trawl_taxasite %>% select(1:2)) %>% mutate(year = factor(year))
TB_NMDS.spp <- as.data.frame(scores(TB_NMDS, display = "species"))

set.seed(123)
vec.spp <- envfit(TB_NMDS, TB_trawl_taxasite %>% select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 1000, na.rm = TRUE)
vec.spp

sig_spp = unname(which(vec.spp[['vectors']]$pvals < 0.001))
spp.scrs = vec.spp[['vectors']]$arrows[sig_spp,] %>% data.frame() %>% rownames_to_column("species")

hulls <- ddply(TB_NMDS.scrs, "year", find_hull)

# NMDS_plot <- 
ggplot(TB_NMDS.scrs) +
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
  # arrow = arrow(length(unit(0.5, "cm")))) +
  ggrepel::geom_text_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic")+
  scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS.scrs$year)))])+
  scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS.scrs$year)))])


set.seed(123)
TB_NMDS_common <- vegan::metaMDS(TB_trawl_commonsite %>% select(-c(1:2)) %>% column_to_rownames("date_id"), distance = 'jaccard', trymax = 1000, autotransform = FALSE)
TB_NMDS_common
plot(TB_NMDS_common)

TB_NMDS_common.scrs <- as.data.frame(scores(TB_NMDS_common, display = 'sites'))
TB_NMDS_common.scrs <- TB_NMDS_common.scrs %>% bind_cols(TB_trawl_commonsite %>% select(1:2)) %>% mutate(year = factor(year))
TB_NMDS_common.spp <- as.data.frame(scores(TB_NMDS_common, display = "species"))

set.seed(123)
vec.common_spp <- envfit(TB_NMDS_common, TB_trawl_commonsite %>% select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 1000, na.rm = TRUE)
vec.common_spp

set.seed(123)
vec.common_env <- envfit(TB_NMDS_common ~ year, TB_trawl_commonsite %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                          month = factor(month)), permutations = 1000, na.rm = TRUE)
vec.common_env

sig.common_spp = unname(which(vec.common_spp[['vectors']]$pvals < 0.001))
common_spp.scrs = vec.common_spp[['vectors']]$arrows[sig.common_spp,] %>% data.frame() %>% rownames_to_column("species")

year_centroid = vec.common_env[['factors']]$centroids %>% data.frame() %>% rownames_to_column('year')
## create an environmental surface of year over the 
year_surf = ordisurf(TB_NMDS_common ~ as.numeric(year), data  = TB_trawl_commonsite %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                     month = factor(month)), na.rm = TRUE)
year_surf.grid = year_surf$grid
year_surf.mite = expand.grid(NMDS1 = year_surf.grid$x, NMDS2 = year_surf.grid$y)
year_surf.mite$z = as.vector(year_surf.grid$z)
year_surf.na = as.data.frame(na.omit(year_surf.mite))


hulls <- ddply(TB_NMDS_common.scrs, "year", find_hull)

common_NMDS_plot <-
  ggplot(TB_NMDS_common.scrs) +
  stat_contour(data = year_surf.na, aes(x = NMDS1, y = NMDS2, z = z), colour = 'grey', binwidth = 1)+
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  # geom_point(data = year_centroid, aes(x = NMDS1, y = NMDS2), size = 3, colour= 'black') +
  # geom_path(data = year_centroid, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
  geom_point(data = common_spp.scrs, aes(x = NMDS1, y = NMDS2), size = 1.5, shape =5, colour = 'black', fill = 'grey')+
  
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
               # arrow = arrow(length(unit(0.5, "cm")))) +
    ggrepel::geom_text_repel(data = common_spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic")+
  scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common.scrs$year)))])+
  scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common.scrs$year)))])

tiff("./figures/TB_ord_comm.tif", res = 600, height = 10, width = 19.1,units = "in", compression ="lzw")
common_NMDS_plot
dev.off()


##addressing outlier points
TB_trawl_outlier_rm <- TB_trawl_commonsite %>%
  filter(date_id != "2018-1" | "2007-7")

set.seed(100)
TB_NMDS_common_out <- vegan::metaMDS(TB_trawl_outlier_rm %>% select(-c(1:2)) %>% column_to_rownames("date_id"), distance = 'jaccard', trymax = 5000, autotransform = FALSE)
TB_NMDS_common_out
plot(TB_NMDS_common_out)

TB_NMDS_common_out.scrs <- as.data.frame(scores(TB_NMDS_common_out, display = 'sites'))
TB_NMDS_common_out.scrs <- TB_NMDS_common_out.scrs %>% bind_cols(TB_trawl_outlier_rm %>% select(1:2)) %>% mutate(year = factor(year))
TB_NMDS_common_out.spp <- as.data.frame(scores(TB_NMDS_common_out, display = "species"))

set.seed(123)
vec.common_out_spp <- envfit(TB_NMDS_common_out, TB_trawl_outlier_rm %>% select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 1000, na.rm = TRUE)
vec.common_out_spp

set.seed(123)
vec.common_out_env <- envfit(TB_NMDS_common_out ~ year, TB_trawl_outlier_rm %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                 month = factor(month)), permutations = 1000, na.rm = TRUE)
vec.common_out_env

sig.common_out_spp = unname(which(vec.common_out_spp[['vectors']]$pvals < 0.001))
common_out_spp.scrs = vec.common_out_spp[['vectors']]$arrows[sig.common_out_spp,] %>% data.frame() %>% rownames_to_column("species")

year_centroid_out = vec.common_out_env[['factors']]$centroids %>% data.frame() %>% rownames_to_column('year')
## create an environmental surface of year over the 
year_surf_out = ordisurf(TB_NMDS_common_out ~ as.numeric(year), data  = TB_trawl_outlier_rm %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                                 month = factor(month)), na.rm = TRUE)
year_surf_out.grid = year_surf_out$grid
year_surf_out.mite = expand.grid(NMDS1 = year_surf_out.grid$x, NMDS2 = year_surf_out.grid$y)
year_surf_out.mite$z = as.vector(year_surf_out.grid$z)
year_surf_out.na = as.data.frame(na.omit(year_surf_out.mite))


hulls <- ddply(TB_NMDS_common_out.scrs, "year", find_hull)

common_NMDS_plot <-
  ggplot(TB_NMDS_common.scrs) +
  stat_contour(data = year_surf.na, aes(x = NMDS1, y = NMDS2, z = z), colour = 'grey', binwidth = 1)+
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  # geom_point(data = year_centroid, aes(x = NMDS1, y = NMDS2), size = 3, colour= 'black') +
  # geom_path(data = year_centroid, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
  geom_point(data = common_spp.scrs, aes(x = NMDS1, y = NMDS2), size = 1.5, shape =5, colour = 'black', fill = 'grey')+
  
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
  # arrow = arrow(length(unit(0.5, "cm")))) +
  ggrepel::geom_text_repel(data = common_spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic")+
  scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common.scrs$year)))])+
  scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common.scrs$year)))])

tiff("./figures/TB_ord_comm.tif", res = 600, height = 10, width = 19.1,units = "in", compression ="lzw")
common_NMDS_plot
dev.off()



# unique_spp = TB_trawl_data %>%
#   mutate(Species = str_remove_all(Species, "(?i)-|mix.*"),
#          Species = str_replace(Species, "-", ","),
#          Species2 = str_remove(Species, ".+;"),
#          Species2 = ifelse(Species2 == Species, NA, Species2),
#          Species = str_remove(Species,";.+"))
#   unlist()%>%
#   str_replace("(?i)-|Mixed", "") 
#   strsplit(gsub("(.+);(.+)","\\1\1\\2", .),"\1")