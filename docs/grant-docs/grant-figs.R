# NSF-BIO-OCE grant figures
source("install-packages.R")
TB_trawl_data <- readRDS("./data/TB_trawl_data.rds")

temp_TS = read_csv(file = "./data/TerrebonneBay_temp_2008-2020.csv", col_type = cols(ArrayID = col_skip(), yeardata = col_skip(), jday= col_skip(),
                                                                               X4 = col_skip(), X5 = col_skip(), X6 = col_skip(), timedata = col_skip(),
                                                                               X8 = col_skip(), `Date/Time` = col_datetime(format = "%m/%d/%Y %H:%M"),
                                                                               `Water Temp C` = col_number(), `Water Temp <U+00B0>F` = col_skip())) %>%
  rename(datetime = "Date/Time", temp_C = "Water Temp C") %>%
  mutate(Date = as.Date(datetime,origin = "1970-01-01")) %>%
  group_by(Date) %>%
  summarise_at(vars("temp_C"), na.rm_mean)

sal_TS = read_csv(file = "./data/TerrebonneBay_salinity_2008-2020.csv", col_type = cols(ArrayID = col_skip(), yeardata = col_skip(), jday = col_skip(), X4 = col_skip(),
                                                                                        X5 = col_skip(), X6 = col_skip(), timedata = col_skip(), X8 = col_skip(), 
                                                                                        X9 = col_datetime(format = "%m/%d/%Y %H:%M"), salinity1_avg = col_number())) %>%
                    rename(datetime = 'X9',sal_psu = 'salinity1_avg') %>%
  na.omit %>%
  mutate(Date = as.Date(datetime,origin = "1970-01-01")) %>%
  group_by(Date) %>%
  summarise_at(vars("sal_psu"), na.rm_mean)

date_TS = seq.Date(as.Date("2000-01-01"), as.Date(Sys.Date()), by = 'days') %>% data.frame %>% setNames(.,c("Date"))
env_summ <- date_TS %>% left_join(temp_TS, all = TRUE) %>% left_join(sal_TS, all = TRUE)

week_TS = seq.Date(as.Date("2000-01-01"), as.Date(Sys.Date()), by = 'weeks') %>% data.frame %>% setNames(.,c("Date"))
env_week = week_TS %>% left_join(temp_TS) %>% left_join(sal_TS)


TB_trawl_samplings <- TB_trawl_data %>% select(Date) %>%
  group_by(Date) %>% distinct %>% ungroup %>%
  mutate(Date = as.Date(Date), 
         Date_num = as.numeric(Date)) %>%
  left_join(env_summ, all = TRUE)

date_breaks = as.numeric(as.Date(paste0(seq(2007,2020, by = 1),"-01-01")))
date_labels = year(as.Date(date_breaks, origin = "1970-01-01"))

ggplot(TB_trawl_samplings, aes(x = Date, y = 0))+ geom_vline(aes(xintercept = as.numeric(ymd(Date)))) +
  scale_x_continuous(breaks = date_breaks, labels = date_labels) + theme_minimal() +
  theme(axis.title = element_blank(), panel.grid = element_blank())

sub1 = subset(env_week, Date < as.Date("2012-01-01"))
sub2 = subset(env_summ, Date > as.Date("2012-01-02"))
sal_means = c(na.rm_mean(sub1$sal_psu, na.rm = TRUE), na.rm_mean(sub2$sal_psu, na.rm = TRUE)) 

sal_plot = ggplot(env_week, aes(x = Date, y = sal_psu)) +
  geom_path(aes(x= as.numeric(ymd(Date)), y = sal_psu), color = "grey", size = 2) +
  annotate("segment", xend = as.numeric(ymd("2007-01-01")), x = as.numeric(ymd("2012-01-01")),
           y = sal_means[1], yend = sal_means[1], size = 1.5, color = "black")+
  annotate('segment', xend = as.numeric(ymd("2016-01-02")), x = as.numeric(ymd("2020-08-12")),
           y = sal_means[2], yend = sal_means[2], size = 1.5, color = "black")+
  # geom_path(aes(x = as.numeric(ymd(Date)), y = temp_C), color = 'red', size = 2) +
  labs(y = "Salinity (PSU)") +
  scale_x_continuous(limits = c(as.Date("2007-01-01"), as.Date(Sys.Date())), breaks = date_breaks, labels = date_labels) + theme_minimal() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16), axis.title.y = element_text(size = 20),
        panel.grid = element_blank())
  
tiff("./figures/sal_plot.tiff", res = 600, height = 7, width = 10, units = "in", compression = 'lzw')
sal_plot
dev.off()

### barcode plot
TB_trawl_samplings <- TB_trawl_data %>% select(Date) %>%
  group_by(Date) %>% distinct %>% ungroup %>%
  mutate(Date = as.Date(Date), 
         Date_num = as.numeric(Date))

core_samplings = data.frame(Date = as.Date(c("2018-06-01", "2019-01-15","2019-04-015","2019-07-01","2019-07-30","2019-08-15","2019-09-15",
                                             "2019-10-15","2020-01-15","2020-04-15","2020-06-15", "2020-08-14")))
date_breaks = as.numeric(as.Date(paste0(seq(2007,2020, by = 1),"-01-01")))
date_labels = year(as.Date(date_breaks, origin = "1970-01-01"))

trawl_sampling_plot = ggplot(TB_trawl_samplings, aes(x = Date, y = 0))+ 
  geom_vline(aes(xintercept = as.numeric(ymd(Date)))) +
  # geom_vline(data = core_samplings, aes(xintercept = as.numeric(ymd(Date))), colour = 'grey', linetype = "dashed", size = 1.2) +
  scale_x_continuous(breaks = date_breaks, labels = date_labels) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_text(size = 20), panel.grid = element_blank()) 

core_sampling_plot = ggplot(TB_trawl_samplings, aes(x = Date, y = 0))+ 
  # geom_vline(aes(xintercept = as.numeric(ymd(Date)))) +
  geom_vline(data = core_samplings, aes(xintercept = as.numeric(ymd(Date))), colour = 'darkgrey', linetype = "dashed", size = 1.2) +
  scale_x_continuous(limits = c(as.numeric(ymd("2007-01-01")),NA), breaks = date_breaks, labels = date_labels) + 
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) 

tiff("./figures/samplings.tiff",res = 450, width = 12, height = 4, units = "in", compression = "lzw")
grid.arrange(core_sampling_plot, trawl_sampling_plot, ncol = 1)
dev.off()
# figure 4
#create time-accumulation df from trawl data
TB_trawl_data %>%
  filter(!duplicated(Common_name)) %>%
  mutate(richness = 1:n(), 
         Date = as.Date(Date),
         richness = sort(richness))-> TB_time_accum

TB_annual_accum <- TB_trawl_data %>%
  mutate(pres = 1) %>% group_by(year) %>%
  filter(!duplicated(Common_name)) %>%
  summarise(richness = sum(pres)) %>%
  mutate(Date = as.Date(as.character(year), format = "%Y", origin = "2000-01-01"))

TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()

# partition beta diversity into nestedness and turnover #
# uses `betapart` package #
TB_time_df <- TB_time_accum %>%
  select(Date, year, richness) %>%
  mutate(data_type = 'continuous') %>%
  bind_rows(TB_annual_accum %>% mutate(data_type = "annual"))

#create moving jaccard calc by year-month
# create a dataframe for time series of distance compared to first sampling 
unique_combinations <- expand.grid(date_A = unique(levels(as.factor(TB_trawl_commonsite$date_id))), 
                                   Date = unique(levels(as.factor(TB_trawl_commonsite$date_id)))) %>%
  filter(date_A == as.character(min(TB_trawl_commonsite$date_id)))

TB_date_jaccard = vegdist((TB_trawl_commonsite %>% select(-year, -month) %>% 
                             column_to_rownames("date_id")), method = "jaccard") %>%
  as.matrix %>% as.data.frame %>% rownames_to_column("date_A") %>% 
  pivot_longer(cols = 2:dim(.)[2], names_to = "Date", values_to = "jaccard") %>%
  inner_join(unique_combinations) %>% filter(date_A != Date) %>% 
  mutate(Date = gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', Date),
         Date = as.Date(paste0(Date,"-01"), format = "%Y-%m-%d"))

# partitioning beta diversity metrics from trawl dataset
# create named list of each sampling date
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>%
    group_split() %>%
    rlang::set_names(names)
}

date_lists = TB_trawl_commonsite %>%
  named_group_split(date_id) %>%
  map(., ~.x %>% select(-year, -month,-date_id))

beta_lists <- lapply(date_lists, function(x,r){ betapart::beta.temp(x = x, y = r, index.family = 'jaccard')}, x = date_lists[[1]]) 
beta_dates <- gsub('^([0-9]{4}-)([0-9]{1})$', '\\10\\2', names(beta_lists))
beta_df <- bind_rows(beta_lists) %>% mutate(Date = as.Date(paste0(beta_dates,"-01"), "%Y-%m-%d")) %>% #select(Date, everything()) %>%
  pivot_longer(-Date, names_to = "partition", values_to = "value")

rich_time <- 
  ggplot(TB_time_df %>% filter(data_type == "annual", year != "2020"), aes(x = Date, y = richness, group = data_type)) +
  geom_line(colour = "black",size = 1.2) +
  geom_point(colour = "grey",size = 1.5) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + 
  scale_y_continuous(name = "Species Richness", limits = c(0,80))+
  scale_colour_manual(name = "", values = viridis(4)[c(2,3)]) +
  theme_minimal() +
  theme(legend.position = "none", legend.justification = c(0, 1), 
        legend.background = element_rect(fill = "transparent", colour = NA),
        axis.title.x = element_blank(),
        legend.spacing.y = unit(-0.1, 'cm'), panel.grid = element_blank()) 

tiff("./figures/annual_richness.tiff", res = 600, height = 5, width = 5, units = "in", compression = 'lzw')
rich_time
dev.off()
# create taxa by site matrix with common names
TB_trawl_commonsite <- TB_trawl_data %>%
  select(year, month, date_id, Common_name) %>%
  group_by(date_id, Common_name) %>%
  distinct() %>% 
  mutate(pres = 1) %>%
  pivot_wider(names_from = Common_name, values_from = pres, values_fill = list(pres = 0)) %>%
  ungroup()
##addressing outlier points
TB_trawl_outlier_rm <- TB_trawl_commonsite %>%
  filter(date_id %ni% c("2018-1","2007-7"))

### Legendre PC with Hellinger
# Hellinger pre-transformation of the species matrix
tb.h <- decostand (TB_trawl_outlier_rm %>% select(-c(1:2)) %>% column_to_rownames("date_id"), "hellinger")
tb.h.pca <- rda(tb.h)
tb.h.pca.summ <- summary(tb.h.pca)

#pull coordinates and loadings
TB <- data.frame(matrix(ncol=2,nrow=105, dimnames=list(NULL, c("PC1", "PC2"))))

df1  <- data.frame(tb.h.pca.summ$sites[,1:2])
TB$PC1 <- df1[,1]
TB$PC2 <- df1[,2]
TB$date_id <- row.names(tb.h.pca.summ$sites)
TB <- left_join(TB, TB_trawl_outlier_rm, by=c("date_id"="date_id"))

Var = factor(TB$year, levels =c("2007","2009","2010","2011",
                                "2012","2013","2014","2015",
                                "2016","2017","2018","2019",
                                "2020"))
#pull loadings
TB.loadings  <- data.frame(tb.h.pca.summ$species[,1:2])
TB.loadings$Species <- row.names(tb.h.pca.summ$species)
#filtering out top species contributions
TB.loadings.top <- TB.loadings %>%
  filter(PC1>quantile(PC1, prob=.99) | 
           PC1<quantile(PC1, prob=.01)| 
           PC2>quantile(PC2, prob=.99)  | 
           PC2<quantile(PC2, prob=.01)) 

#analyses
# tb.h.dl <- dist(tb.h)
# pc_test<- betadisper(tb.h.dl, TB$year)
# anova(pc_test)
# permutest(pc_test, pairwise = TRUE)

TB.centroids <- TB %>%
  group_by(year) %>%
  summarise(
    meanPC1 = mean(PC1),
    meanPC2 = mean(PC2)
  )
TB$year <- as.integer(TB$year)

#figure
pca_grant = ggplot(data = TB, aes(x=PC1, y=PC2, color=as.factor(year))) + 
  geom_point(aes(color=as.factor(year)), size =3, alpha = 0.6)+
  geom_path(data=TB.centroids, aes(x=meanPC1, y=meanPC2), arrow=arrow(), size=1, color="grey20")+
  geom_label(data=TB.centroids, aes(x=meanPC1, y=meanPC2, label=year), color="grey20")+
  # ggrepel::geom_text_repel(data=TB.loadings.top, aes(x=PC1, y=PC2, label=Species), color="grey30")+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  labs(color = "Year") +
  theme_minimal() +
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = c(1,0), legend.justification = c(1,0),
        legend.background = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank())+
  scale_colour_viridis(discrete = TRUE)

png("./figures/trawl_PCA_grant.png", res = 450, height = 5, width = 5, units = 'in')
pca_grant
dev.off()
