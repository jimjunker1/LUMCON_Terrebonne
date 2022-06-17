SAC_full_plot <- function() {
SAC_full <- ggplot(spp_mod.df, aes(x = trawls, y = richness)) + geom_path(size = 2, colour = 'red') +
  geom_ribbon(aes(x = trawls, ymin = (richness - 1.96*richness_sd), ymax = (richness + 1.96*richness_sd)), colour = "red", fill = "red", alpha = 0.5) +
  scale_y_continuous(name = expression(Sigma~"species"), breaks = seq(0,roundUP(max(spp_mod.df$richness), to = 10), by = 10), labels = seq(0,roundUP(max(spp_mod.df$richness), to = 10), by = 10)) +
  scale_x_continuous(name = "Trawls", breaks = seq(0,roundUP(max(spp_mod.df$trawls), to = 10), by = 20));SAC_full

tiff(filename = "./figures/trawl_specaccum_full.tiff", height = 5, width = 5, units = "in",res = 450, compression = "lzw")
SAC_full
dev.off()
}