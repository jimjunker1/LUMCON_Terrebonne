
# Full species accumulation curve (SAC)

##create taxa by site matrix
TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

sp1 <- vegan::specaccum(TB_trawl_commonsite %>% select(-year, -month, -date_id))
spp_mod.df <- data.frame(trawls = sp1[["sites"]], richness = sp1[["richness"]], richness_sd = sp1[["sd"]])
# plot(sp1, ci.type = "poly", col = "nlue", lwd = 2, ci.lty = 0)
sp2 <- vegan::specaccum(TB_trawl_taxasite %>% select(-year, -month, -date_id), 'random')
boxplot(sp2, col = "yellow", add = TRUE, pch = "+")

SAC_full <- ggplot(spp_mod.df, aes(x = trawls, y = richness)) + geom_path(size = 2, colour = 'red') +
  geom_ribbon(aes(x = trawls, ymin = (richness - 1.96*richness_sd), ymax = (richness + 1.96*richness_sd)), colour = "red", fill = "red", alpha = 0.5) +
  scale_y_continuous(name = expression(Sigma~"species"), breaks = seq(0,roundUP(max(spp_mod.df$richness), to = 10), by = 10), labels = seq(0,roundUP(max(spp_mod.df$richness), to = 10), by = 10)) +
  scale_x_continuous(name = "Trawls", breaks = seq(0,roundUP(max(spp_mod.df$trawls), to = 10), by = 20));SAC_full

tiff(filename = "./figures/trawl_specaccum_full.tiff", height = 5, width = 5, units = "in",res = 450, compression = "lzw")
SAC_full
dev.off()

## SAC for each year
TB_trawl_commondate <- TB_trawl_data %>%
  select(year, Date, Common_name) %>%
  group_by(Date, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()


TB_trawl_taxasite_yr <- TB_trawl_commondate %>%
  group_split(year) %>%
  map(., ~.x %>%  select(-year, -Date) %>% vegan::specaccum(.) %>%
       rlist::list.remove(names(.) %ni% c("sites","richness","sd")) %>%
        bind_cols) %>% 
  map2(., TB_trawl_commondate %>% group_split(year), ~..1 %>% 
         mutate(year = ..2 %>% select(year) %>% unlist %>% as.factor)) %>%
  bind_rows

SAC_yr = ggplot(TB_trawl_taxasite_yr, aes(x = sites, y = richness, group = year, colour = year)) +
  geom_path(size = 1.5) +# geom_ribbon(aes(x = sites, ymin = (richness - 1.96*sd), ymax = (richness + 1.96*sd), fill = year), alpha = 0.2) +
  scale_y_continuous(name = expression(Sigma~"species"),breaks = seq(0, roundUP(max(TB_trawl_taxasite_yr$richness), to = 10), by = 10), labels = seq(0, roundUP(max(TB_trawl_taxasite_yr$richness), to = 10), by = 10)) +
  scale_x_continuous(name = "Trawls", breaks = seq(0, roundUP(max(TB_trawl_taxasite_yr$sites), to = 10), by = 5)) +
  scale_colour_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE);SAC_yr


### need to break this out by individual trawls not by month ##
tiff(file = "./figures/trawl_specaccum_yr.tiff", height = 5, width = 5, units = "in", res = 450, compression = "lzw")
SAC_yr
dev.off()
