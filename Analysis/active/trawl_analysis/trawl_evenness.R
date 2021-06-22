# trawl_evenness Chao and Ricotta 
source("./R/Evenness.R")
## ++++ Helper functions ++++ ##
#'
#'
#'
#'@param x the species abundance/frequency/productivity vector
tax_q_profile <- function(x, q, name1){
  x <- as.data.frame(x)
  x$q.order <-as.character(q)
  x1 <- x %>% pivot_longer(-q.order, names_to = 'Date', values_to = 'evenness') %>%
    dplyr::mutate(Date = as.Date(Date, format = "%Y-%m-%d"),
                  year = year(Date))
  ggplot(x1, aes(q.order, evenness))+
    geom_line(aes(color = year, group = Date), size = 1.1)+
    theme_tufte(ticks = TRUE)+
    geom_rangeframe(sides = "lb")+
    theme(legend.position = c(1,1),
          legend.justification = c("right", "top"), 
          legend.title = element_blank(),
          legend.spacing.y = unit(0.0003,'cm'),
          legend.key.size = unit(0.4,'cm'),
          strip.text = element_blank())+
    scale_x_discrete(breaks=seq(0, 3, 0.5))+
    scale_colour_viridis()+
    # theme(legend.key.width = unit(2,"cm"))+
    # theme(plot.title = element_text(size=20, face="bold.italic",hjust = 0.5))+
    ggtitle(name1)+
    xlab("Diversity order q")
  #ylim(c(0, 1))
}

evenness_profile_function = function(x, q.seq, evenness.type = c("E1", "E2", "E3", "E4", "E5", "E6"),...){
  name1 = evenness.type
  ind_evenness1 <- apply(x, 2, function(x) new_fun(x = x, q.order = q.seq, evenness.type = evenness.type))
  ind_evenness2 <- array(unlist(ind_evenness1), c(length(q.seq), dim(x)[2], length(evenness.type)))
  ind_evenness3 <- aperm(ind_evenness2, c(1,3,2))
  dimnames(ind_evenness3)[[1]] <- paste0("q = ", q.seq)
  dimnames(ind_evenness3)[[3]] <- colnames(x)
  dimnames(ind_evenness3)[[2]] <- name1
  return(ind_evenness3 = ind_evenness3)
}
## ++++ End Helper functions ++++ ##
##create taxa by site matrix
TB_trawl_taxasite_t <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  column_to_rownames("Date") %>%
  as.matrix %>% t

debugonce(evenness_profile_function)
debugonce(new_fun)
x = evenness_profile_function(TB_trawl_taxasite_t, q.seq = seq(0, 2, 0.05))
debugonce(tax_q_profile)
evenness_plot = tax_q_profile(x[, 3, ],seq(0, 2, 0.05), "E5")+ ylim(c(0, 1)) +facet_wrap(~year) + theme(plot.title = element_blank())
# 
