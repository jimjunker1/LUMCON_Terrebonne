

# Calculated dissimilarity for all temporal distances based on pres-absence

TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(year, month, date_id, species_mod) %>%
  dplyr::filter(year != "2021") %>%
  group_by(date_id, species_mod) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = species_mod, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

