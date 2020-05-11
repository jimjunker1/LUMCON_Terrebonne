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
  