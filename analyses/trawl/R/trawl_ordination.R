## terrebonne bay trawl working analysis
## works with TB_trawl_data
# source("install-packages.R")#installs all necessary packages
source(here::here("datascript.R"))#imports data with some minor cleanup


## summarise the monthly sampling effort
## look at this over time
## 

## talk to murt about getting trawl locations
##

## look into Ernest community time analysis from 

## automate a second state analysis to look at the correlation of distance metrics
## then run an mds on the new correlation

##create taxa by site matrix
TB_trawl_taxasite <- TB_trawl_data %>%
  dplyr::select(year, month, date_id, species_mod) %>%
  group_by(date_id, species_mod) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = species_mod, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

set.seed(123)
TB_NMDS <- vegan::metaMDS(TB_trawl_taxasite %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), distance = "jaccard",
                          try = 1000, trymax = 3000, maxit = 9999, autotransform = FALSE, parallel = 3)
TB_NMDS
vegan::stressplot(TB_NMDS)
plot(TB_NMDS)

TB_NMDS.scrs <- as.data.frame(scores(TB_NMDS, display = 'sites'))
TB_NMDS.scrs <- TB_NMDS.scrs %>% bind_cols(TB_trawl_taxasite %>% dplyr::select(1,3:5)) %>% mutate(year = factor(year))
TB_NMDS.spp <- as.data.frame(scores(TB_NMDS, display = "species"))

set.seed(123)
vec.spp <- envfit(TB_NMDS, TB_trawl_taxasite %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 1000, na.rm = TRUE)
vec.spp

sig_spp = unname(which(vec.spp[['vectors']]$pvals < 0.001))
spp.scrs = vec.spp[['vectors']]$arrows[sig_spp,] %>% data.frame() %>% rownames_to_column("species")

hulls <- ddply(TB_NMDS.scrs, "year", find_hull)

# NMDS_plot <- 
ggplot(TB_NMDS.scrs) +
  # geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
  # arrow = arrow(length(unit(0.5, "cm")))) +
  ggrepel::geom_text_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic")+
  scale_fill_viridis(name = "Year", discrete = TRUE)+
  scale_color_viridis(name = "Year", discrete = TRUE)
  
# create taxa by site matrix with common names
TB_trawl_commonsite <- TB_trawl_data %>%
  dplyr::select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

set.seed(123)
TB_NMDS_common <- vegan::metaMDS(TB_trawl_commonsite %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), distance = 'jaccard', trymax = 1000, autotransform = FALSE)
TB_NMDS_common
plot(TB_NMDS_common)

TB_NMDS_common.scrs <- as.data.frame(scores(TB_NMDS_common, display = 'sites'))
TB_NMDS_common.scrs <- TB_NMDS_common.scrs %>% bind_cols(TB_trawl_commonsite %>% dplyr::select(1:2)) %>% mutate(year = factor(year))
TB_NMDS_common.spp <- as.data.frame(scores(TB_NMDS_common, display = "species"))

set.seed(123)
vec.common_spp <- envfit(TB_NMDS_common, TB_trawl_commonsite %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 1000, na.rm = TRUE)
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
  filter(date_id %ni% c("2018-1","2007-7"))

set.seed(100)
TB_NMDS_common_out <- vegan::metaMDS(TB_trawl_outlier_rm %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), distance = 'jaccard', trymax = 5000, autotransform = FALSE)
TB_NMDS_common_out
plot(TB_NMDS_common_out)

TB_NMDS_common_out.scrs <- as.data.frame(scores(TB_NMDS_common_out, display = 'sites'))
TB_NMDS_common_out.scrs <- TB_NMDS_common_out.scrs %>% bind_cols(TB_trawl_outlier_rm %>% dplyr::select(1:2)) %>% mutate(year = factor(year), month = factor(month))
TB_NMDS_common_out.spp <- as.data.frame(scores(TB_NMDS_common_out, display = "species"))

set.seed(123)
vec.common_out_spp <- envfit(TB_NMDS_common_out, TB_trawl_outlier_rm %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 5000, na.rm = TRUE)
vec.common_out_spp

set.seed(123)
vec.common_out_env <- envfit(TB_NMDS_common_out ~ year+month, TB_trawl_outlier_rm %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                         month = factor(month)), permutations = 5000, na.rm = TRUE)
vec.common_out_env

sig.common_out_spp = unname(which(vec.common_out_spp[['vectors']]$pvals < 0.0002))
common_out_spp.scrs = vec.common_out_spp[['vectors']]$arrows[sig.common_out_spp,] %>% data.frame() %>% rownames_to_column("species")

year_centroid_out = vec.common_out_env[['factors']]$centroids %>% data.frame %>% 
  rownames_to_column('date_id') %>% filter(grepl("year", date_id, ignore.case = TRUE)) %>%
                                             mutate(year = str_remove(date_id, 'year'))

month_centroid_out = vec.common_out_env[['factors']]$centroids %>% data.frame %>% 
  rownames_to_column("date_id") %>% filter(grepl("month", date_id, ignore.case = TRUE)) %>% 
  mutate( month = str_remove(date_id, "month"))

## create an environmental surface of year over the 
year_surf_out = ordisurf(TB_NMDS_common_out ~ as.numeric(year), data  = TB_trawl_outlier_rm %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                                         month = factor(month)), na.rm = TRUE)
year_surf_out.grid = year_surf_out$grid
year_surf_out.mite = expand.grid(NMDS1 = year_surf_out.grid$x, NMDS2 = year_surf_out.grid$y)
year_surf_out.mite$z = as.vector(year_surf_out.grid$z)
year_surf_out.na = as.data.frame(na.omit(year_surf_out.mite))


yr_hulls <- ddply(TB_NMDS_common_out.scrs, "year", find_hull)
mth_hulls <- ddply(TB_NMDS_common_out.scrs, "month", find_hull)
# yr_common_out_NMDS_plot <-
  ggplot(TB_NMDS_common_out.scrs %>% dplyr::select(-month))+
  geom_point(aes(x = NMDS1, y = NMDS2, colour = year), size = 2) +
  geom_point(data = year_centroid_out, aes(x = NMDS1, y = NMDS2, fill = year), shape = 21, size = 3, colour = "black") +
    geom_path(data = year_centroid_out, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
  scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])+
  scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])
  
    
#mnth_common_out_NMDS_plot <-
  ggplot(TB_NMDS_common_out.scrs %>% dplyr::select(-year))+
    geom_point(aes( x= NMDS1, y = NMDS2, colour = month), size = 2) +
    geom_point(data = month_centroid_out, aes(x = NMDS1, y = NMDS2, fill = month), shape = 21, size = 3, colour = "black") +
    geom_point(data = month_centroid_out, aes(x = NMDS1, y = NMDS2), size = 1, colour = "darkgrey") #+
    scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])+
    scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])
    
  
  common_out_NMDS_plot <-
  ggplot(TB_NMDS_common_out.scrs %>% dplyr::select(-month)) +
  # stat_contour(data = year_surf_out.na, aes(x = NMDS1, y = NMDS2, z = z), colour = 'grey', binwidth = 1)+
  # geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = year), size = 2)+
  geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  geom_point(data = year_centroid_out, aes(x = NMDS1, y = NMDS2, fill = year), shape = 21, size = 3, colour= 'black') +
  geom_path(data = year_centroid_out, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
  geom_point(data = common_out_spp.scrs, aes(x = NMDS1, y = NMDS2), size = 1.5, shape =5, colour = 'black', fill = 'grey')+
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
  # arrow = arrow(length(unit(0.5, "cm")))) +
  ggrepel::geom_text_repel(data = common_out_spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic", alpha = 0.5)+
  annotate('text', x = -1, y = 1, label = "stress = 0.26", hjust = 0)+
  annotate('text', x = -1, y = 0.9, label = "distance: jaccard", hjust = 0)#+
  # scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])+
  # scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_NMDS_common_out.scrs$year)))])

tiff("./figures/TB_ord_comm_outrm.tif", res = 600, height = 4.5, width = 6,units = "in", compression ="lzw")
common_out_NMDS_plot
dev.off()

#### working with abundance data ####
TB_trawl_taxasite_N <- TB_trawl_data %>%
  filter(date_id %ni% c("2018-1","2007-7"))%>%
  dplyr::select(year, month, date_id, Common_name, Abundance) %>%
  na.omit() %>%
  group_by(year, month, date_id, Common_name) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  group_by(date_id) %>%
  mutate(rel_abun = Abundance/sum(Abundance)) %>%
  dplyr::select(year, month, date_id, Common_name, rel_abun) %>%
  pivot_wider(names_from = Common_name, values_from = rel_abun, values_fill = list(rel_abun = 0)) %>%
  ungroup() 

#Hellinger pre-transformation of the species matrix
TB_trawl_taxasite_Ntrans <- decostand (TB_trawl_taxasite_N %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), "hellinger")

set.seed(123)
TB_relN_NMDS <- vegan::metaMDS(TB_trawl_taxasite_Ntrans, distance = "bray", trymax = 1000, autotransform = FALSE)
TB_relN_NMDS 
plot(TB_relN_NMDS )

TB_relN_NMDS.scrs <- as.data.frame(scores(TB_relN_NMDS, display = 'sites'))
TB_relN_NMDS.scrs <- TB_relN_NMDS.scrs %>% bind_cols(TB_trawl_taxasite_N  %>% dplyr::select(1:2)) %>% mutate(year = factor(year))
TB_relN_NMDS.spp <- as.data.frame(scores(TB_relN_NMDS, display = "species"))

set.seed(123)
vec.spp <- envfit(TB_relN_NMDS, TB_trawl_taxasite_N  %>% dplyr::select(-c(1:2)) %>% column_to_rownames("date_id"), permutations = 5000, na.rm = TRUE)
vec.spp

set.seed(123)
vec.common_env <- envfit(TB_relN_NMDS ~ year+month, TB_trawl_taxasite_N %>% column_to_rownames("date_id") %>% mutate(year = factor(year),
                                                                                                                 month = factor(month)), permutations = 5000, na.rm = TRUE)
vec.common_env

sig_spp = unname(which(vec.spp[['vectors']]$pvals < 0.001))
spp.scrs = vec.spp[['vectors']]$arrows[sig_spp,] %>% data.frame() %>% rownames_to_column("species")

y_hulls <- ddply(TB_relN_NMDS.scrs, "year", find_hull)
m_hulls <- ddply(TB_relN_NMDS.scrs, 'month', find_hull)

year_centroid = vec.common_env[['factors']]$centroids %>% data.frame() %>% rownames_to_column('date_id') %>% filter(grepl("year", date_id)) %>% mutate(year = str_remove(date_id, 'year'))
month_centroid = vec.common_env[['factors']]$centroids %>% data.frame() %>% rownames_to_column('date_id') %>% filter(grepl("month", date_id)) %>% mutate(year = str_remove(date_id, 'month'))

relN_NMDS_plot <- 
ggplot(TB_relN_NMDS.scrs) +
  annotate('text', x = 0, y = 0, label = "draft:months only", angle = 45, alpha = 0.1, size = 10)+
  # geom_point(data = year_centroid, aes(x = NMDS1, y = NMDS2, fill = year), shape = 21, size = 3, colour= 'black') +
  # geom_path(data = year_centroid, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
  geom_point(data = month_centroid, aes(x = NMDS1, y = NMDS2, fill = year), shape = 22, size = 3, colour= 'black') +
  geom_path(data = month_centroid, aes(x = NMDS1, y = NMDS2), size =1 , colour = 'darkgrey')+
    # stat_contour(data = year_surf_out.na, aes(x = NMDS1, y = NMDS2, z = z), colour = 'grey', binwidth = 1)+
  # geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = year, colour = year), alpha = 0.3) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = year), size = 2, alpha = 0.5)+
  # geom_point(data = hulls, aes(x = NMDS1, y = NMDS2,  fill = year, colour = year), size = 2, shape = 21) +
  geom_point(data = spp.scrs, aes(x = NMDS1, y = NMDS2), size = 1.5, shape =5, colour = 'black', fill = 'grey')+
  # geom_segment(data = spp.scrs, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), size = 1, 
  # arrow = arrow(length(unit(0.5, "cm")))) +
  ggrepel::geom_text_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = species), fontface = "italic", alpha = 0.5)+
  annotate('text', x = -1, y = 1, label = "stress = 0.26", hjust = 0)+
  annotate('text', x = -1, y = 0.9, label = "distance: bray", hjust = 0)#+
  # scale_fill_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_relN_NMDS.scrs$year)))])+
  # scale_colour_manual(name = "Year", values = ocecolors[['temperature']][seq(1,256, length.out = length(levels(TB_relN_NMDS.scrs$year)))])

tiff("./figures/TB_ord_comm_relN.tif", res = 600, height = 4.5, width = 6,units = "in", compression ="lzw")
relN_NMDS_plot
dev.off()
