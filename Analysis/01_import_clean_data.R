# Initiation code for loading packages and raw data manipulation

load_packages = function(){
  if(!require("pacman")) install.packages("pacman")
  library(pacman)
  package.list <- c("twitteR","data.table", "RCurl","plyr","tidyverse","furrr",
                    "tictoc","chron","lubridate","httr","TTR", 
                    "grid","gridExtra", "ggridges",
                    "viridis", "broom","bbmle","ggthemes", "ggeffects")
  p_load(char = package.list, install = T)
  rm("package.list")
  #  len_freq_path = getURL("https://raw.githubusercontent.com/jimjunker1/secprod_workflow/master/len_freq/len_freq_function.R",ssl.verifypeer = FALSE)
  #  len_freq <<- eval(parse(text = len_freq_path))
  cbbPalette <<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  load("./object_files/ocecolors.rda")
  theme_mod <<- theme_bw() %+replace% theme(panel.grid = element_blank())
  theme_black <<- function() {theme_bw() %+replace% theme(panel.background = element_rect(fill = 'transparent', colour = NA),panel.grid = element_blank(), axis.ticks = element_line(color = 'white'),
                                                          axis.title = element_text(color = 'white'), axis.text = element_text(color = 'white'), plot.background =element_rect(fill = 'transparent', colour = NA),
                                                          panel.border = element_rect(fill = NA, colour = 'white'))}
  C_to_overkt <<- function(a){1/(8.61733*10^-5*(a+273.15))}#overkt function
  overkt_to_C <<- function(a){1/(a*(8.61733*10^-5)) - 273.15}
  C_to_overkt_stand15 <<- function(a){(1/(8.61733e-5*(15+273.15)) - (1/(8.61733e-5*(a+273.15))))}
  # overkt_to_C <<- Vectorize(function(x) {
  #   VGAM::reciprocal(x*(8.61733*10^-5))-273.15
  # })
  Mode <<- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  na.rm_mean <<- function(...,na.rm=TRUE){mean(c(...),na.rm=na.rm)}
  options(stringsAsFactors = F)

  mass_corrPB <<- function(a,b){a*(b^-0.25)}
  mass_corrBIO <<- function(a,b){a*(b^0.25)}
  
  temp_corrPB <<- function(a,b){a*(b^0.65)}
  temp_corrBIO <<- function(a,b){a*(b^-0.65)}
  
  '%ni%' <<- Negate('%in%')
  
  myspread <<- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
      unite(temp, !!keyq, variable) %>%
      spread(temp, value)
  }
  get_legend<<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
}
load_packages()

# raw data manipulation
# extracts data from raw data sheets to populate analysis datasets
# outputs .RDS dataset objects of currated data to './data' folder for further analysis
# can be run each time or re-run just when raw data is added

data_manipulation = function() {
  
  
}