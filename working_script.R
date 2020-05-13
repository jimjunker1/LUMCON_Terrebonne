TB_trawl_samplings <- TB_trawl_data %>% select(Date) %>%
  group_by(Date) %>% distinct %>% ungroup %>%
  mutate(Date = as.Date(Date), 
         Date_num = as.numeric(Date))

date_breaks = as.numeric(as.Date(paste0(seq(2007,2020, by = 1),"-01-01")))
date_labels = year(as.Date(date_breaks, origin = "1970-01-01"))

ggplot(TB_trawl_samplings, aes(x = Date, y = 0))+ geom_vline(aes(xintercept = as.numeric(ymd(Date)))) +
  scale_x_continuous(breaks = date_breaks, labels = date_labels) +
  theme(axis.title = element_blank()) +
  ggtitle('Sampling events')
                                      

TB_month_hist <- TB_trawl_samplings %>%
  mutate(month_num = lubridate::month(Date), 
         month_name = lubridate::month(Date, label = TRUE)) 

  ggplot(TB_month_hist, aes(month_num)) + geom_bar(width = 1) +
    scale_x_continuous(breaks = seq(1,12, by = 1), labels = lubridate::month(seq(1,12,by = 1), label = TRUE)) +
    theme_minimal() +
    theme(axis.title = element_blank())
  
  
#### Combine relative abundance by year and plot how this ####
  
  TB_ann_relN = TB_trawl_data %>% select(year, Common_name, Abundance) %>% na.omit %>%
    group_by(year, Common_name) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
    mutate(rel_N = Abundance/sum(Abundance, na.rm = TRUE)) %>%
    left_join(TB_trawl_data %>% 
                # Turn filter back on to use 2007 as reference year
                # filter(year == 2007) %>%
                na.omit %>%
                group_by(Common_name) %>%
                summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
                mutate(rel_N = Abundance/sum(Abundance, na.rm = TRUE)) %>%
                mutate(rank = dense_rank(desc(rel_N))) %>%
                select(Common_name, rank))
  
  TB_trawl_data %>% filter(year == 2007) %>%
    group_by(Common_name) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE))
    
  TB_relN_ann = TB_trawl_data %>% select(year, Common_name, Abundance) %>% na.omit %>%
    group_by(year, Common_name) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
    mutate(rel_N = Abundance/sum(Abundance, na.rm = TRUE),
           rank = dense_rank(desc(rel_N))) %>%
    ungroup() %>%  mutate(year = as.factor(year))
# widen out data
  
  TB_relN_ann_wide <- TB_relN_ann %>%
    pivot_wider(id_cols = c(year,Common_name), names_from = year, values_from = rank) %>%
    group_by(Common_name)
  