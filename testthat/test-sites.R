library(testthat)
context("checks that input sites are valid")

TB_2018_core_meta<<-read_csv(file = "./data/metadata/TB_2018_Core-Metadata.csv", 
                             col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
                                              `Van Veen Grab` = "i", Core = "i"))
TB_2019_core_meta<<- read_csv(file = "./data/metadata/TB_2019_Core-Metadata.csv",
                              col_types = cols(.default = "?", Date = col_date(format = "%m/%d/%Y"), Time = "i",
                                               `Van Veen Grab` = "i", Core = "i"))
core_data = bind_rows(TB_2018_core_meta, TB_2019_core_meta)

site_list <<- c("TB2A","TB2B","TB2C",
                "TB3A","TB3V","TB3C",
                "TB4A","TB4B","TB4C")

test_that("sites are valid", {
  
  expect_true(all(core_data$Station %in% site_list))
  
})