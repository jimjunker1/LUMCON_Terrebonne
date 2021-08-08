here::here()
here::i_am("packages.R")
# Install pacman if it isn't already installed
if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")

# Install packages required for analysis
# Add any packages required for your data cleaning and manipulation to the
# `p_load` function call below
# Do not remove the packages already listed here
# they are important for running the livedat repository
# OG package list: pacman::p_load(git2r, httr, semver, testthat, yaml)

pacman::p_load(
  git2r, httr, semver, testthat, yaml, here, EML, data.table,
  RCurl, plyr, tidyverse, furrr, googlesheets4,
  tictoc, chron, lubridate, httr, TTR, grid, egg,
  gridExtra, ggridges, iNEXT, vegan, rlist, pipeR,
  viridis, broom, bbmle, ggthemes, ggeffects, bsts,
  rfishbase, auk, betapart, plot3D, ncdf4, forecast,
  taxize, hillR, codyn, gridExtra, fuzzySim
)

devtools::install_github(c("jimjunker1/junkR")) # , "ropensci/rfishbase"))
# library(junkR)
# devtools::install_github('MoBiodiv/mobr')
library(mobr)


theme_mod <<- function() {
  theme_bw() %+replace% theme(panel.grid = element_blank())
}
theme_set(theme_mod())
rm(theme_mod)
find_hull <<- function(df) df[chull(df$NMDS1, df$NMDS2), ]

"%ni%" <<- Negate("%in%")
# bootstrapped sampling for estimating variability withing groups
boot_sample_groups <- function(abund_mat, group_var, qD, n_perm, ...) {
  boot_list <- vector("list", n_perm)

  # sample rows and calculate abundance vector
  sample_hill <- function(boot_n, ...) {
    sample_dat <- by(abund_mat,
      INDICES = abund_mat[, group_var], FUN = sample_frac,
      replace = TRUE
    )
    class(sample_dat) <- "list"
    sample_dat <- bind_rows(sample_dat)

    agg_groups <- stats::aggregate(sample_dat[, -1],
      by = sample_dat[, 1],
      FUN = "sum"
    )

    group_col <- agg_groups[, group_var]
    qD_col <- paste0("q", parse(text = qD))
    boot_col <- paste0("boot", parse(text = boot_n))

    dat_groups <- agg_groups %>%
      apply(., 1, function(x) hillR::hill_taxa(comm = x, q = qD)) %>%
      bind_cols(group_col, .) %>%
      setNames(., c(group_var, qD_col)) %>%
      dplyr::mutate(boot_id = boot_col)

    return(dat_groups)
  }
  # debugonce(sample_hill)
  boot_list <- purrr::map(1:n_perm, ~ sample_hill(boot_n = .x)) %>% bind_rows()
  return(boot_list)
}

# bootstrapped sampling for estimating variability withing groups
boot_gini_groups <- function(abund_mat, group_var, n_perm, ...) {
  boot_list <- vector("list", n_perm)

  # sample rows and calculate abundance vector
  sample_even <- function(boot_n, ...) {
    sample_dat <- by(abund_mat,
      INDICES = abund_mat[, group_var], FUN = sample_frac,
      replace = TRUE
    )
    class(sample_dat) <- "list"
    sample_dat <- bind_rows(sample_dat)

    agg_groups <- stats::aggregate(sample_dat[, -1],
      by = sample_dat[, 1],
      FUN = "sum"
    )

    group_col <- agg_groups[, group_var]
    boot_col <- paste0("boot", parse(text = boot_n))
    agg_dat <- agg_groups[, -1]

    dat_groups <- agg_dat %>%
      apply(., 1, gini_even) %>%
      t() %>%
      data.frame() %>%
      bind_cols(group_col, .) %>%
      setNames(., c(group_var, "gini_raw", "gini_norm")) %>%
      dplyr::mutate(boot_id = boot_col)

    return(dat_groups)
  }
  # debugonce(sample_even)
  boot_list <- purrr::map(1:n_perm, ~ sample_even(boot_n = .x)) %>% bind_rows()
  return(boot_list)
}

pct_rare <- function(x, rare_thres) {
  S <- sum(x > 0)
  if (S > 0) {
    N <- sum(x)
    if (rare_thres == "N/S") {
      rare_thres <- N / S
      100 * (sum(x[x > 0] <= rare_thres) / S)
    } else {
      100 * (sum(x[x > 0] <= (rare_thres * N)) / S)
    }
  } else {
    0
  }
}

# yday_summary <<- function(data,...){
#   if("yday" %in% colnames(data)){
#     data_summ = data %>%
#       group_by(yday) %>%
#       summarise_if(is.numeric, list(~na.rm_mean(.,na.rm = TRUE), ~max(.,na.rm = TRUE), ~min(., na.rm = TRUE))) %>%
#       ungroup()
#   } ifelse{
#
#   }
#
# }
