env_plot <- function(){
env_facet_plot <- ggplotGrob(ggplot(env_long, aes(x = datetime, y = value, group = env_var, colour = env_var)) +
                               geom_line() + scale_colour_viridis_d() +
                               facet_grid(rows = vars(env_var), scales = "free_y")+
  theme(axis.title = element_blank(), legend.position = "none"))

png(file = "./figures/env_facet_plot.png", res = 400, height = 10, width = 5, units = "in")
grid.draw(env_facet_plot)
dev.off()
}
