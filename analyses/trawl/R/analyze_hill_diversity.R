##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

here::i_am("analyses/trawl/R/analyze_hill_diversity.R")
library(hillR)

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(Date, Common_name, Abundance) %>%
  dplyr::filter(year(Date) %ni% c("2007","2020","2021")) %>%
  group_by(Date, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  na.omit %>% group_by(Date) %>%
  # remove trawls with less than 10 individuals 
  dplyr::filter(sum(Abundance) >=10) %>%
  group_by(Date, Common_name) %>%
  pivot_wider(names_from = Common_name, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  ungroup() %>%
  dplyr::mutate(julian_date = julian(Date, origin = as.Date("2009-01-01")),
                biweekly = ceiling(julian_date/30),
                month = month(Date),
                year = year(Date)) %>%
  dplyr::select(julian_date, Date,year, month, biweekly, everything())

## Estimated qD = 2 for all communities
biweekly_sampling = TB_trawl_taxasite %>%
  group_by(biweekly) %>%
  dplyr::summarise(n = n()) %>%
  right_join(TB_trawl_taxasite %>% dplyr::select(Date, ))

trawl_qD = hillR::hill_taxa(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly), q = 2) %>% 
  data.frame %>% dplyr::mutate(Date = TB_trawl_taxasite$Date,
                               qD= 2) %>% 
  bind_rows(
hillR::hill_taxa(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly), q= 1) %>%
  data.frame %>% dplyr::mutate(Date = TB_trawl_taxasite$Date,
                               qD = 1)) %>%
  bind_rows(hillR::hill_taxa(TB_trawl_taxasite %>% dplyr::select(-julian_date:-biweekly), q= 0) %>%
              data.frame %>% dplyr::mutate(Date = TB_trawl_taxasite$Date,
                                           qD = 0))

ggplot(trawl_qD, aes(x =Date, y =., group = qD, color = qD)) +
  geom_point(size = 2, alpha = 0.5)+
  geom_line(alpha = 0.5)+
  geom_smooth( se = FALSE, span = 0.1)+
  scale_y_continuous(limits = c(0,10))





plot(site_hill2)
  ## ++++ Helper functions ++++ ##
  #'
  #'
  #'
  #'@param x the species abundance/frequency/productivity vector
  tax_q_profile <- function(x, q, name1){
    x <- as.data.frame(x)
    x$q.order <-as.character(q)
    x1 <- x %>% pivot_longer(-q.order, names_to = 'site', values_to = 'evenness') %>%
      dplyr::mutate(site = factor(site, levels = names(stream_order_list)))
    ggplot(x1, aes(q.order, evenness))+
      geom_line(aes(color = site, group = site), size = 1.1)+
      theme_tufte(ticks = TRUE)+
      geom_rangeframe(sides = "lb")+
      theme(legend.position = c(1,1),
            legend.justification = c("right", "top"), 
            legend.title = element_blank(),
            legend.spacing.y = unit(0.0003,'cm'),
            legend.key.size = unit(0.4,'cm'),
            strip.text = element_blank())+
      scale_x_discrete(breaks=seq(0, 3, 0.5))+
      scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels)+
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
    dimnames(ind_evenness3)[[2]] <- colnames(x)
    dimnames(ind_evenness3)[[3]] <- name1
    return(ind_evenness3 = ind_evenness3)
  }
  ## ++++ End Helper functions ++++ ##

  # ann_spp_flux_list = ann_spp_flux_boots %>% named_group_split(boot_id) %>%
  #   map(~.x %>%
  #         dplyr::select(-boot_id) %>%
  #         group_by(site) %>%
  #         pivot_wider( names_from = taxon, values_from = flux_mg_m_y, values_fill = 0) %>%
  #         ungroup %>%
  #         dplyr::mutate(site = factor(site, levels = names(stream_order_list))) %>%
  #         column_to_rownames('site') %>%
  #         .[names(stream_order_list),] %>%
  #         as.matrix) #%>%
    # map2(., list(seq(0:5)), ~.x %>% hillR::hill_taxa(.x, q = .y))
  
  ann_spp_flux = ann_spp_flux_summary %>% data.frame %>%
    dplyr::select(site, taxon, flux_mg_m_y_mean) %>%
    pivot_wider(names_from = taxon, values_from = flux_mg_m_y_mean, values_fill = 0) %>%
    ungroup %>%
    dplyr::mutate(site = factor(site, levels=names(stream_order_list))) %>%
    column_to_rownames('site') %>%
    .[names(stream_order_list), ] %>%
    as.matrix %>% t()
  
  
  # debugonce(evenness_profile_function)#;debugonce(new_fun)
  x = evenness_profile_function(ann_spp_flux, q.seq = seq(0, 2, 0.05))
  
  evenness_plot = tax_q_profile(x[, , 3],seq(0, 2, 0.05), "E5")+ ylim(c(0, 1)) + theme(plot.title = element_blank())
# 
#   hill_taxa_0 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 0) %>% .[names(stream_order_list)]
#   hill_taxa_1 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 1)
#   hill_taxa_2 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 2)
#   hill_taxa_3 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 3)
#   hill_taxa_4 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 4)
#   hill_taxa_5 = hillR::hill_taxa(ann_spp_flux_list[[1]], q = 5)
#   
#   
#   
#   data.frame(site = names(hill_taxa_0),
#              hill0 = hill_taxa_0, 
#              hill1 = hill_taxa_1,
#              hill2 = hill_taxa_2,
#              hill3 = hill_taxa_3,
#              hill4 = hill_taxa_4,
#              hill5 = hill_taxa_5) %>%
#     pivot_longer(-site, names_to = 'q', values_to = 'value') %>%
#     dplyr::mutate(site = factor(site, levels = names(stream_order_list))) %>%
#     ggplot(aes(x = q, y = log(value), group = site, colour = site)) + geom_point() + geom_line()
  
  ####################################################################################
  #
  # Example for (1). Alpine species example (See Figure 1 for data and Figure 2 for output)
  #
  ####################################################################################

  ##########caculate evenness##########
  # data = ann_spp_flux_list[[1]] %>% t
  # q <- seq(0, 3, 0.05)
  # name1 <- c("E1", "E2", "E3", "E4", "E5", "E6")
  # ind_evenness1 <- lapply(list(data[, 1], data[, 2], data[, 3], data[,4], data[,5], data[,6]), function(x) new_fun(x = x, q.order = q))
  # ind_evenness2 <- array(unlist(ind_evenness1), c(length(q), 6, 6))
  # ind_evenness3 <- aperm(ind_evenness2, c(1, 3, 2))
  # dimnames(ind_evenness3)[[1]] <- paste0("q = ", q)
  # dimnames(ind_evenness3)[[2]] <- colnames(data)
  # dimnames(ind_evenness3)[[3]] <- name1
  # Gini_indices <- apply(data,2,gini_even)# compute the Gini evenness indices for each assemblage.
  # 
  # ##########plot ##############
  # # debugonce(tax_q_profile)
  # ggarrange(tax_q_profile(ind_evenness3[, , 1],q, name1[1])+ylim(c(0, 1)),
  #           tax_q_profile(ind_evenness3[, , 2],q, name1[2])+ylim(c(0, 1)),
  #           tax_q_profile(ind_evenness3[, , 3],q, name1[3])+ylim(c(0, 1)),
  #           tax_q_profile(ind_evenness3[, , 4],q, name1[4])+ylim(c(0, 1)),
  #           tax_q_profile(ind_evenness3[, , 5],q, name1[5])+ylim(c(0, 1)),
  #           tax_q_profile(ind_evenness3[, , 6],q, name1[6])+ylim(c(0, 1)),
  #           ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
  return(evenness_plot)
}
