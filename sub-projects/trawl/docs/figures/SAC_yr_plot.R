SAC_yr_plot <- function(){
SAC_yr = ggplot(TB_trawl_taxasite_yr, aes(x = sites, y = richness, group = year, colour = year)) +
  geom_path(size = 1.5) +# geom_ribbon(aes(x = sites, ymin = (richness - 1.96*sd), ymax = (richness + 1.96*sd), fill = year), alpha = 0.2) +
  scale_y_continuous(name = expression(Sigma~"species"),breaks = seq(0, roundUP(max(TB_trawl_taxasite_yr$richness), to = 10), by = 10), labels = seq(0, roundUP(max(TB_trawl_taxasite_yr$richness), to = 10), by = 10)) +
  scale_x_continuous(name = "Trawls", breaks = seq(0, roundUP(max(TB_trawl_taxasite_yr$sites), to = 10), by = 5)) +
  scale_colour_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE);SAC_yr

### need to break this out by individual trawls not by month ##
tiff(file = "./figures/trawl_specaccum_yr.tiff", height = 5, width = 5, units = "in", res = 450, compression = "lzw")
SAC_yr
dev.off()
}