source(here::here('datascript.R'))

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::mutate(julian_date = round(julian(Date, origin = "2007-01-05")),
                biweekly = ceiling(julian_date/30),
                month = month(Date),
                year = year(Date)) %>%
  dplyr::select(julian_date, Date,year, month, biweekly, everything())



# trawl_rarefaction analysis

trawl_mob_in = mobr::make_mob_in(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly),
                           TB_trawl_taxasite %>% dplyr::select(julian_date:biweekly),
                           coord_names = "")
trawl_mod_stats = mobr::get_mob_stats(trawl_mob_in)

mobr::plot_rarefaction(trawl_mob_in,'year', method = "IBR", lwd = 3, leg_loc = 'none' )
