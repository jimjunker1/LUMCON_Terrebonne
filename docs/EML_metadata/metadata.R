#  metadata creation file  #
##  Hydro monitoring data  ##
attributes <-
  tibble::tribble(
    ~attributeName, ~attributeDefinition,                                                 ~formatString, ~definition,        ~unit,   ~numberType,
    "datetime",    "datetime stamp of collection. 2020-01-01 00:00:01",                 "YYYY-MM-DD HH:MM:SS",            NA, NA,       NA,
    "year",       "year, 2012",                                                         "YYYY",        NA,                 NA,       NA,
    "day",        "Julian day. Range: 1 - 209.",                                      "DDD",         NA,                 NA,       NA,
    "i.flag",     "is variable Real, Interpolated or Bad (character/factor)",           NA,            NA,                 NA,       NA,
    "variable",   "what variable being measured (character/factor)", NA,            NA,                 NA,       NA,
    "value.i",    "value of measured variable on year/day/hour.min. Real",       NA,            NA,                 NA,       NA)

i.flag <- c(R = "real",
            I = "interpolated",
            B = "bad")
variable <- c(
  temp_F  = "Temperature in Fahrenheit",
  temp_C = "Temperature in Celcius",
  do_mg_L = "Dissolved Oxygen Concentration mg/L",
  spCon_uS_cm  = "Specific conductivity microSeimens/cm",
  sal_psu = "Salinity in Practical Salinity Units (PSU; g/kg)"
  )

value.i <- c(
  temp_F  = paste0(expression(degree*"F")),
  temp_C = paste0(expression(degree*"C")),
  do_mg_L = "[DO]; mg/L",
  spCon_us_cm = paste0(expression(mu*"Siemens/cm")),
  sal_psu = "Practical Salinity Units"
)

## Write these into the data.frame format
factors <- rbind(
  data.frame(
    attributeName = "i.flag",
    code = names(i.flag),
    definition = unname(i.flag)
  ),
  data.frame(
    attributeName = "variable",
    code = names(variable),
    definition = unname(variable)
  ),
  data.frame(
    attributeName = "value.i",
    code = names(value.i),
    definition = unname(value.i)
  )
)


attributeList <- EML::set_attributes(attributes, factors, col_classes = c("Date", "Date", "Date", "factor", "factor", "numeric"))

physical <- set_physical("hydro.csv")

dataTable <- list(
  entityName = "hydro.csv",
  entityDescription = "Historical hydrologic monitoring data",
  physical = physical,
  attributeList = attributeList)

# # Geographic coverage of trawl data  ##
geographicDescription <- "Terrebonne Bay, Chauvin Louisiana"


coverage <-
  set_coverage(begin = '2000-10-11', end = '2020-01-01',
               geographicDescription = geographicDescription,
               west = -90.60, east = -90.40,
               north = 29.30, south = 29.08,
               altitudeMin = 0, altitudeMaximum = 0,
               altitudeUnits = "meter")

list.files(path = "./data/", "*temp_*", full.names = TRUE)
methods_file <- list.files(path = "./text/EML_metadata/","hydro.*\\.md", full.names = TRUE)
methods <- EML::set_methods(methods_file)
