####################################

# helper functions

#####################################

library(tidyverse)
library(magrittr)   # used for the %$% pipe
library(lubridate)  # datetime formatting
library(here)       # file management
library(ncdf4)
library(tidync)
library(patchwork)  # plotting
library(rstan)
library(tidybayes)

rstan::rstan_options(auto_write = TRUE)



#' Function to transform climate data (real or climate model output)
#' Does things like transform to anomaly space, conversion to water year, and optional aggregation to some yearly average
#' @param monthly_climate dataframe of monthly climate in long format, with columns `year`, `month`, and `clim`
#' @param clim_name string, name of climate variable we wish to use
#' @param w_agg vector of length 12 which sums to 1. Used if aggregating to yearly values for baseline models. defaults to NULL, which returns monthly values.
#' @returns a dataframe of log-transformed, anomalized climate data. If w_agg is null, then dimensions will match `clim_df` with columns `water_year`, `month`, and `clim`. If `w_agg` is not NULL, then the dataframe will be of yearly values, with columns `water_year` and `clim`
#' @example process_climate(sim_precip_df, name = "log_precip")
process_climate <- function(monthly_climate, clim_name, w_agg = NULL, anom = T, standardize_names = T, water_year = T){

  if(standardize_names){
    #using `all_of()` is a funny trick to force R to use clim_name in the local env
    # see https://github.com/tidyverse/dplyr/issues/6975
    monthly_df <- monthly_climate %>%
      rename(clim = all_of(clim_name))
  } else{
    monthly_df <- monthly_climate
  }


  ## convert to water year --------
  if(water_year){

    clim_df <- monthly_df %>%
      mutate(water_year = ifelse(month %in% c('Oct', 'Nov', 'Dec'), year + 1, year)
             # water_year = ifelse(month %in% c('Sep', 'Oct', 'Nov', 'Dec'), year + 1, year)
      ) %>%
      relocate(water_year) %>%
      # delete first and last years (incomplete water year data)
      group_by(water_year) %>%
      mutate(num_months = n()) %>%
      filter(num_months == 12) %>%
      ungroup() %>%
      select(-year) %>%
      rename(year = water_year)
  } else{
    clim_df <- monthly_df %>%
      group_by(year) %>%
      mutate(num_months = n()) %>%
      filter(num_months == 12)
  }


  # model expects the order to be calendar year
  month_lookup <- tibble(
    month = month.abb,
    month_no = 1:12
  )

  clim_df <- clim_df %>%
    left_join(month_lookup, by = 'month') %>%
    group_by(year) %>%
    arrange(month_no, .by_group = T) %>%
    ungroup()


  # turn weights into lookup table for convenience
  w_lookup <- tibble(
    w_agg = w_agg,
    month = month.abb
  )




  ## aggregate monthly data to yearly ------

  if(!is.null(w_agg)){
    yearly_climate <- clim_df %>%
      left_join(w_lookup, by = 'month') %>%
      group_by(year) %>%
      summarize(
        # q: should we use clim instead of clim_anom here.
        clim = sum(w_agg * clim)
      )

    return(yearly_climate)
  }

  ## transform to standardized anomaly space (z scores) -----

  if(anom){

    clim_means <- clim_df %>%
      group_by(month) %>%
      summarize(
        clim_mean = mean(clim),
        clim_var = sd(clim)
      )

    clim_df <- clim_df %>%
      left_join(clim_means, by = 'month') %>%
      mutate(clim = (clim - clim_mean)/clim_var)

  }


  if(standardize_names){
    # future functions expect `clim` to be the column name.
    # future functions expet `year` to be the column name
    clim_df <- clim_df %>%
      select(year,
             month,
             clim)
  }

  clim_df

}

#' a set of "trees" as a linear function of PRECIP
#' @param sim_df is a dataset which has columns ,,,,
#' @param w_true are the 12 monthly weights that the tree will be sensitive to
#' @return a list containing 2 dataframes: the first with monthly growth, the second with simulated tree rings (yearly)
make_tree_clim_linear <- function(sim_df, n_trees = NULL,
                                  w_true = rep(1/12, 12), sigma2_w = 1/2,
                                  beta0 = rnorm(1), beta1 = rnorm(1), beta_age = runif(1), take_log = F){

  if(take_log){
    sim1_trw <- sim_df %>%
      group_by(year) %>%
      summarize(
        clim_yearw = log(sum(w_true * clim) + 0.001), ## add log here!!!
        trw_sim = beta0 + beta1 * clim_yearw,
        .groups = "drop")

  } else{
    sim1_trw <- sim_df %>%
      group_by(year) %>%
      summarize(
        clim_yearw = sum(w_true * clim),
        trw_sim = beta0 + beta1 * clim_yearw,
        .groups = "drop")
  }



  sim1_df <- sim1_trw


  n_trees <- n_trees
  sim_trws <- map_dfr(1:n_trees, ~ tibble(
    tree_id = .x,
    year = sim1_df$year,
    age = 1:nrow(sim1_df),
    clim = sim1_df$clim_yearw
  ) %>%
    mutate(age_std = (age - mean(age)) / sd(age)) %>%
    mutate(trw_mean = sim1_df$trw_sim + beta_age * age_std) %>%
    mutate(
      trw_scale = sim1_df %$%
        scale(trw_mean)[,1] * sqrt(1 - sigma2_w) + rnorm(nrow(sim1_df), 0, sqrt(sigma2_w))
    ) %>%
    mutate(trw = trw_scale + abs(min(trw_scale))) # force positive trw. note this changes expected intercept in fit.

  )


  list(sim_growth = sim_df,
       sim_trw = sim_trws
  )
}


#'
#' @param clim_df output of `process_climate()`
#' @param n_trees positive integer, number of pseudoproxy trees to create
#' @param snr positive number (typically less than 2). Signal to noise ratio of the pseudoproxies
#' @param w_true vector of length 12 which sums to 1. Defines the growing season of the pseudoproxies.
#' @param beta0 intercept of the trw ~ climate relationship
#' @param beta1 slope of the trw ~ climate relationship
#' @returns: list of dataframes (length n_trees) with columns `tree_id`, `year`, `trw`, and `age`.
#' @example format_trw_linear(sim_precip_df, 10, 1, 1500, rep(1/12, 12))
format_trw_linear <- function(clim_df, n_trees, snr, start_year_calib,
                               w_true, beta0 = 1, beta1 = 2, beta_age = 0, include_growth = F, take_log = F){

  ### generative part: pick a first year of data for each tree
  first_year_obs <- sample(min(clim_df$year):start_year_calib, n_trees, replace = T)

  ## we can also pick out a last year of obs for the trees, if we don't want them all to live to current day.
  # last_year_obs <- sample((min(sim_df$year) + 2):max(sim_df$year))

  start_year_reconstruction <- min(first_year_obs)


  noise <- 1 / (1 + snr^2)
  fake_trees_full <- make_tree_clim_linear(clim_df, w_true = w_true,
                                   n_trees = n_trees, sigma2_w = noise,
                                   beta0 = beta0, beta1 = beta1, beta_age = beta_age,
                                   take_log = take_log)

  sim_trw_df <- fake_trees_full$sim_trw

  # remove data so that the pseudoproxies don't all start/end at the same time.
  fake_trees <- sim_trw_df %>%
    group_by(tree_id) %>%
    group_split() %>%
    imap(
      ~.x %>%
        filter(year >= first_year_obs[.y]) %>%
        mutate(age = 1:nrow(.)) %>%
        select(tree_id, year, trw, age)
    )

  if(!include_growth){
    return(fake_trees)
  }

  list(
    trw = fake_trees,
    monthly_growth = fake_trees_full$monthly_growth
  )



}


#'
#' @param clim_df output of `process_climate()`
#' @param n_trees positive integer, number of pseudoproxy trees to create
#' @param snr positive number (typically less than 2). Signal to noise ratio of the pseudoproxies
#' @param month_nums a vector containing the months we want the tree to grow in.
#' @returns: list of dataframes (length n_trees) with columns `tree_id`, `year`, `trw`, and `age`.
#' @example format_trw_linear(sim_precip_df, 10, 1, 1500, rep(1/12, 12))
format_trw_vslite <- function(clim_df, n_trees, snr, start_year_calib, month_nums = c(12,1,2),
                              sm_limited = T,
                              include_extras = F, ...){

  ### generative part: pick a first year of data for each tree
  first_year_obs <- sample(min(clim_df$year):start_year_calib, n_trees, replace = T)

  ## we can also pick out a last year of obs for the trees, if we don't want them all to live to current day.
  # last_year_obs <- sample((min(sim_df$year) + 2):max(sim_df$year))

  start_year_reconstruction <- min(first_year_obs)


  noise <- 1 / (1 + snr^2)

  fake_trees_full <- make_tree_vslite(clim_df, n_trees = n_trees,
                                      month_nums = month_nums,
                                      sigma2_w = noise, sm_limited = sm_limited, ...)

  sim_trw_df <- fake_trees_full$sim_trw

  # remove data so that the pseudoproxies don't all start/end at the same time.

  fake_trees <- sim_trw_df %>%
    group_by(tree_id) %>%
    group_split() %>%
    imap(
      ~.x %>%
        filter(year >= first_year_obs[.y]) %>%
        mutate(age = 1:nrow(.)) %>%
        select(tree_id, year, trw, age)
    )

  if(!include_extras){
    return(fake_trees)
  }

  list(
    trw = fake_trees,
    monthly_growth = fake_trees_full$monthly_growth,
    cor_monthly = fake_trees_full$cor_monthly,
    cor_by_month = fake_trees_full$cor_by_month,
    prop_moisture_limited = fake_trees_full$prop_moisture_limited
  )



}


#' Wrapper around VSLiteR::leakybucket.monthly(), to get leaky bucket soil moisture from the VSLiteR Package
#' Default parameters are taken from the VSLiteR package
#' @param sim_df monthly dataframe with columns temp (in deg C), precip (in mm), and soil_m (in v/v)
#'
get_leaky_soilm <- function(sim_df,
                            Mmax = 0.76, Mmin = 0.01, alph = 0.093, m.th = 4.886,
                            mu.th = 5.8, rootd = 1000, M0 = 0.2){
  require(VSLiteR)

  ## getting things in the right formatting for VSLiteR
  temp_formatted <- sim_df %>%
    select(year, month, temp_c) %>%
    pivot_wider(names_from = year, values_from = temp_c) %>%
    select(-month) %>%
    as.matrix()

  precip_formatted <- sim_df %>%
    select(year, month, precip_mmmonth) %>%
    pivot_wider(names_from = year, values_from = precip_mmmonth) %>%
    select(-month) %>%
    as.matrix()

  phi <- unique(sim_df$lat)

  syear <- as.numeric(min(sim_df$year))
  eyear <- as.numeric(max(sim_df$year))

  ## --------

  ## getting soil moisture from leaky bucket
  soil_m_raw <- VSLiteR::leakybucket.monthly(syear, eyear, phi, temp_formatted, precip_formatted,
                                             Mmax, Mmin, alph, m.th, mu.th, rootd, M0)

  soil_m <- soil_m_raw %>%
    as_tibble() %>%
    set_colnames(syear:eyear) %>%
    rownames_to_column('month_no') %>%
    pivot_longer(-month_no,
                 names_to = 'year',
                 values_to = 'M')  %>%
    mutate(month_no = as.numeric(month_no),
           year = as.numeric(year)) %>%
    arrange(year, month_no) %>%
    mutate(month = month.abb[month_no]) %>%
    select(
      year,
      # month_no,
      month,
      M
    )

  soil_m
}

#' a set of trees with the same vs-lite parameters
#' uses the `VSLiteR` package (defaults also taken from here)
#' @param phi Latitude (deg N)
#' @param sim_df monthly dataframe with columns temp (in deg C), precip (in mm), and soil_m (in v/v)
make_tree_vslite <- function(sim_df, n_trees = NULL,
                             month_nums = 1:12, sigma2_w = 1/2, sm_limited = T,
                             T1 = 8, T2 = 23, M1 = 0.01, M2 = 0.05,
                             hydroclim = "P", Mmax = 0.76, Mmin = 0.01, alph = 0.093, m.th = 4.886,
                             mu.th = 5.8, rootd = 1000, M0 = 0.2
){
  require(VSLiteR)

  ## getting things in the right formatting for VSLiteR
  temp_formatted <- sim_df %>%
    select(year, month, temp_c) %>%
    pivot_wider(names_from = year, values_from = temp_c) %>%
    select(-month) %>%
    as.matrix()

  precip_formatted <- sim_df %>%
    select(year, month, precip_mmmonth) %>%
    pivot_wider(names_from = year, values_from = precip_mmmonth) %>%
    select(-month) %>%
    as.matrix()

  phi <- unique(sim_df$lat)

  syear <- as.numeric(min(sim_df$year))
  eyear <- as.numeric(max(sim_df$year))

  ## --------

  ## getting soil moisture from leaky bucket
  soil_m_raw <- VSLiteR::leakybucket.monthly(syear, eyear, phi, temp_formatted, precip_formatted,
                           Mmax, Mmin, alph, m.th, mu.th, rootd, M0)

  soil_m <- soil_m_raw %>%
    as_tibble() %>%
    set_colnames(syear:eyear) %>%
    rownames_to_column('month_no') %>%
    pivot_longer(-month_no,
                 names_to = 'year',
                 values_to = 'M')



  # try to make location-agnostic "moisture sensitive trees" by using quantiles
  q_M <- quantile(soil_m$M, c(0, .1, .5, .75, .95))
  q_T <- quantile(sim_df$temp_c, c(0, .1, .25, .5, .75, .9))


  # just one ad-hoc way of making moisture limited trees
  if(sm_limited){
    M1 = q_M["10%"]
    M2 = q_M["95%"]
    T1 = q_T["0%"]
    T2 = q_T["10%"]
  }

  ## Running VSLite
  # I_0 = 0 uses all 12 months. it's just useful for the sanity check
  #
  sim_vslite <- VSLiteR::VSLite(syear, eyear, phi, temp_formatted, precip_formatted, I_0 = 0,
                                M1 = M1, M2 = M2,
                                T1 = M2, T2 = T2
                                )


  ## -----------
  ### formatting VSLite output ------------
  gt <- sim_vslite$gT %>%
    as_tibble() %>%
    rownames_to_column('month_no') %>%
    pivot_longer(-month_no,
                 names_to = 'year',
                 values_to = 'gT')

  gm <- sim_vslite$gM %>%
    as_tibble() %>%
    set_colnames(colnames(sim_vslite$gT)) %>%
    rownames_to_column('month_no') %>%
    pivot_longer(-month_no,
                 names_to = 'year',
                 values_to = 'gM')


  ge <- sim_vslite$gE %>%
    as_tibble() %>%
    rownames_to_column('month_no') %>%
    rename('gE' = V1)

  # combine all the growth components
  # calculate total growth
  g_df <- gt %>%
    left_join(gm, by = c('year', 'month_no')) %>%
    left_join(ge, by = 'month_no') %>%
    mutate(g_tot = gE * pmin(gT, gM)) %>%
    mutate(month = as.integer(month_no),
           year = as.integer(year))

  # quick test ------------
  # manually calculate trw_sim and compare to built in function output
  trw_norm_manual <- g_df %>%
    group_by(year) %>%
    summarise(trw = sum(g_tot)) %>%
    mutate(trw_norm = (trw - mean(trw)) / sd(trw))

  # VSLite-calculated.
  trw_norm_fnc <- sim_vslite$trw %>%
    as_tibble() %>%
    pivot_longer(everything(),
                 names_to = 'year',
                 values_to = 'trw')

  if(!identical(trw_norm_manual$trw_norm, trw_norm_fnc$trw)){
    stop('something went wrong with the TRW calculation!')
  }


  #-------
  ## further restrict the months to make the correct answer more obvious
  g_final_df <- g_df %>%
    filter(month %in% month_nums)


  ## get some diagnostics about how well we created moisture-limited trees

  # what proportion of months were limited by moisture?
  prop_moisture_limited <- g_final_df %$%
    mean(gM < gT)

  # what's the correlation between monthly growth and precip?
  # in the linear ppes, this is 1 (because this is before adding noise)
  cor_precip <- g_final_df %>%
    left_join(sim_df, by = c('year', 'month' = 'month_no')) %$%
    cor(g_tot, precip_mmmonth)

  cor_temp <- g_final_df %>%
    left_join(sim_df, by = c('year', 'month' = 'month_no')) %$%
    cor(g_tot, temp_c)

  # ---------
  # calculate trw in restricted months
  trw_norm_selectedmonths <- g_df %>%
    filter(month %in% month_nums) %>%
    group_by(year) %>%
    summarise(trw = sum(g_tot)) %>%
    mutate(trw_norm = (trw - mean(trw)) / sd(trw))


  # ---------
  # get cor(trw, month precip) for each month
  # as a proxy for the "optimal reconstruction season"

  cor_months_df <- sim_df %>%
    left_join(trw_norm_selectedmonths, by = 'year') %>%
    group_by(month_no) %>%
    summarize(cor_mon = cor(precip_mmmonth, trw),
              cor_mon_temp = cor(temp_c, trw))


  #--------

  sim_trw_df <- map_dfr(1:n_trees, ~ tibble(
    tree_id = .x,
    year = trw_norm_selectedmonths$year,
    trw = trw_norm_selectedmonths %$%
      trw_norm * sqrt(1 - sigma2_w) + rnorm(nrow(trw_norm_selectedmonths), 0, sqrt(sigma2_w)))
  )

  list(sim_trw = sim_trw_df,
       monthly_growth = g_df,
       cor_monthly = cor_precip,
       cor_monthly_t = cor_temp,
       cor_by_month = cor_months_df,
       prop_moisture_limited = prop_moisture_limited
  )


}




#' get trees in a standardized format to accomodate ppe and real data
#' @description reads in the blue oaks data, formatted as the output of `trw_to_long()`
#' @param start_year_calib integer between 1896 and 2004. if not null, selects only trees with overlap with the calibration period (start_year_calib - 2004)
#' @returns: list of dataframes (length of n_trees) with columns `tree_id`, `year`, `trw`, and `age`
format_trw_lobos <- function(start_year_calib = NULL){
  lobos_trees <- read_csv(here('data', 'blueoaks_long.csv')) %>%
    filter(site %in% c('loslobos', 'figeroa')) %>%
    dplyr::select(tree_id, year, age, trw)

  # filter out trees that have data processing issues
  lobos_trees <- lobos_trees %>%
    filter(!str_detect(tree_id, 'FIG103')) %>%
    filter(!str_detect(tree_id, 'FIG105'))

  if(!is.null(start_year_calib)){
    lobos_trees <- lobos_trees %>%
      as_tibble() %>%
      group_by(tree_id) %>%
      mutate(min_year = min(year), max_year = max(year)) %>%
      filter(max_year > start_year_calib)%>%# make sure there's overlap
      ungroup()
  }

  lobos_trees %>%
    group_by(tree_id) %>%
    group_split()
}

#' do any tree transformations we're thinking of (e.g. log transform, filtering for length)
#' @param trw_list output of `format_trw_*()`
#' @returns list of dataframes, matching dimension of `trw_list`, with added column `trw_mod`
process_trw <- function(trw_list){

  # test that no years are missing
  test_seq <- trw_list %>%
    bind_rows() %>%
    group_by(tree_id) %>%
    summarize(test = all.equal(year, min(year, na.rm = T):max(year, na.rm = T))) %>%
    pull(test) %>%
    all()

  if(!test_seq){
    stop('some trees have gaps in years')
  }

  # workflow wants the column to be called "trw"
  trw_list %>%
    map(~.x %>%
          rename(trw_original = trw) %>%
          # mutate(trw = sqrt(trw_original)) %>%
          # mutate(trw = trw_original + 0.01)
          # mutate(trw = log(trw_original + 100)) %>%
          mutate(trw = (trw_original - mean(trw_original)) / sd(trw_original))
          # mutate(trw = (trw_original - mean(trw_original)) / sd(trw_original))
        )


}


#' set up all modelling parameters common to both yearly and monthly models.
#' @param clim_df output of process_climate()
#' @param trw_dfs output of process_trw_*()
#' @param start_year_reconstruction year to start reconstructing. If null, defaults to first year of tree data. Useful to change for testing, to make computation quicker. if null, defaults to first year of tree data
#' @param start_year_calib year to start using observed climate data. Years after start_year_cailb define the "training period", and years of data in clim_df prior to start_year_calib define the validation period.
#'
stan_setup <- function(clim_df, trw_dfs,
                       start_year_calib, start_year_reconstruct = NULL,
                       w_param = rep(1, 12), month_model = T, month_means = NULL, month_sds = NULL
                       ){

  n_trees <- length(trw_dfs)

  # get first year of data for each tree
  first_years_obs <- trw_dfs %>%
    map_dbl(
      ~.x %>%
        filter(!is.na(trw)) %>%
        pull(year) %>%
        min()
    )

  # get last year of data for each tree
  last_years_obs <- trw_dfs %>%
    map_dbl(
      ~.x %>%
        filter(!is.na(trw)) %>%
        pull(year) %>%
        max()
    )

  # check when we want to start reconstruction
  if(is.null(start_year_reconstruct)){
    start_year_reconstruct <- min(first_years_obs)
  }

  last_year_reconstruct <- max(last_years_obs)

  # pull out monthly calibration climate
  clim_obs <- clim_df %>%
    filter(between(year,
                   start_year_calib,
                   last_year_reconstruct)
           ) %>%
    pull(clim)

  # for evaluation, get monthly validation climate
  # e.g. vector of earliest tree year -> first instrumental observed year
  clim_unobs <- clim_df %>%
    filter(year < start_year_calib,
           year >= start_year_reconstruct) %>%
    pull(clim)

  # concatenate the trees into one large vector
  trw_list <- trw_dfs %>%
    map(~.x %>%
          pull(trw))

  trw <- unlist(trw_list)
  trw_lengths <- trw_list %>%
    map_dbl(~length(.x))


  age_list <- trw_dfs %>%
    map(~.x %>%
          pull(age))

  age <- unlist(age_list)


  if(month_model){
    df_fit <- list(
      N_years = max(last_years_obs) - min(first_years_obs) + 1,
      N_years_obs = length(clim_obs) / 12,
      N_trees = n_trees,
      trw = trw,
      age = age,
      clim_obs = clim_obs,
      s = trw_lengths,
      start_years = first_years_obs,
      w_param = w_param,
      w_djf = c(1/3, 1/3, rep(0, 9), 1/3),
      month_means_mod = month_means,
      month_sds_mod = month_sds
    )
  } else{
    df_fit <- list(
      N_years = max(last_years_obs) - min(first_years_obs) + 1,
      N_years_obs = length(clim_obs),
      N_trees = n_trees,
      trw = trw,
      age = age,
      clim_obs = clim_obs,
      s = trw_lengths,
      start_years = first_years_obs
    )
  }



  df_fit
}



#' set up all modelling parameters common to both yearly and monthly models.
#' but only in calibration time.
#' @param clim_df output of process_climate()
#' @param trw_dfs output of process_trw_*()
#' @param start_year_reconstruction year to start reconstructing. If null, defaults to first year of tree data. Useful to change for testing, to make computation quicker. if null, defaults to first year of tree data
#' @param start_year_calib year to start using observed climate data. Years after start_year_cailb define the "training period", and years of data in clim_df prior to start_year_calib define the validation period.
#'
stan_setup_calib <- function(clim_df, trw_dfs,
                       start_year_calib, w_param = rep(1, 12)
){

  # filter just calibration time period
  trw_dfs <- trw_dfs %>%
    map(
      ~.x %>%
        filter(year >= start_year_calib)
    ) %>%
    discard(function(x) nrow(x) == 0)

  n_trees <- length(trw_dfs)

  # get first year of data for each tree
  first_years_obs <- trw_dfs %>%
    map_dbl(
      ~.x %>%
        # filter(!is.na(trw)) %>%
        pull(year) %>%
        min()
    )

  # get last year of data for each tree
  last_years_obs <- trw_dfs %>%
    map_dbl(
      ~.x %>%
        # filter(!is.na(trw)) %>%
        pull(year) %>%
        max()
    )


  last_year_reconstruct <- max(last_years_obs)

  # pull out monthly calibration climate
  clim_obs <- clim_df %>%
    filter(between(year,
                   start_year_calib,
                   last_year_reconstruct)
    ) %>%
    pull(clim)

  # concatenate the trees into one large vector
  trw_list <- trw_dfs %>%
    map(~.x %>%
          # filter(!is.na(trw)) %>%
          pull(trw))

  trw <- unlist(trw_list)
  trw_lengths <- trw_list %>%
    map_dbl(~length(.x))


  age_list <- trw_dfs %>%
    map(~.x %>%
          # filter(!is.na(trw)) %>%
          pull(age))

  age <- unlist(age_list)

  # halfyear_prior <- c(10,10,rep(1, 6), rep(10, 4))
  # w_param <- rep(1, 12)

  df_fit <- list(
    N_years = max(last_years_obs) - min(first_years_obs) + 1,
    N_years_obs = length(clim_obs) / 12,
    N_trees = n_trees,
    trw = trw,
    age = age,
    clim_obs = clim_obs,
    s = trw_lengths,
    start_years = first_years_obs,
    w_param = w_param
  )


  df_fit
}

#' copy of VSLiteR::VSLite, so that I can play with months
#' Almost all of this code is copied directly from VSLiteR::VSLite!!
VSLite_mod <- function (syear, eyear, phi, T, P, T1 = 8, T2 = 23, M1 = 0.01,
          M2 = 0.05, Mmax = 0.76, Mmin = 0.01, alph = 0.093, m.th = 4.886,
          mu.th = 5.8, rootd = 1000, M0 = 0.2, substep = 0, I_0 = 1,
          I_f = 12, hydroclim = "P")
{
  nyrs <- length(syear:eyear)
  Gr <- gT <- gM <- M <- potEv <- matrix(NA, 12, nyrs)
  if (hydroclim == "M") {
    M = P
  }
  else {
    if (substep == 1) {
      M <- leakybucket.submonthly(syear, eyear, phi, T,
                                  P, Mmax, Mmin, alph, m.th, mu.th, rootd, M0)
    }
    else {
      M <- leakybucket.monthly(syear, eyear, phi, T, P,
                               Mmax, Mmin, alph, m.th, mu.th, rootd, M0)
    }
    if (substep != 1 && substep != 0) {
      cat("'substep' param must either be set to 1 or 0.")
      return
    }
  }
  gE <- compute.gE(phi)
  gT <- std.ramp(T, T1, T2)
  gM <- std.ramp(M, M1, M2)
  Gr <- kronecker(matrix(1, 1, nyrs), gE) * pmin(gT, gM)
  width <- matrix(NA, nyrs, 1)
  if (phi > 0) {
    if (I_0 < 0) {
      startmo <- 13 + I_0
      endmo <- I_f
      width[1] <- sum(Gr[1:endmo, 1]) + sum(rowMeans(Gr[startmo:12,
      ]))
      for (cyear in 2:nyrs) {
        width[cyear] <- colSums(Gr[startmo:12, cyear -
                                     1]) + colSums(Gr[1:endmo, cyear])
      }
    }
    else {
      startmo <- I_0 + 1
      endmo <- I_f
      width <- colSums(Gr[startmo:endmo, ])
    }
  }
  if (phi < 0) {
    startmo <- 7 + I_0
    endmo <- I_f - 6
    for (cyear in 1:(nyrs - 1)) {
      width(cyear) <- sum(Gr[startmo:12, cyear]) + sum(Gr[1:endmo,
                                                          cyear + 1])
    }
    width[nyrs] <- sum(Gr[startmo:12, nyrs]) + sum(rowMeans(Gr[1:endmo,
    ]))
  }
  trw <- t((width - mean(width))/sd(width))
  out <- list(trw = trw, gT = gT, gM = gM, gE = gE, M = M,
              potEv = potEv, sample.mean.width = mean(width), sample.std.width = sd(width))
  return(out)
}



# source: https://github.com/betanalpha/knitr_case_studies/blob/master/fitting_the_cauchy/stan_utility.R
check_energy <- function(fit, quiet=FALSE) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      if (!quiet) print(sprintf('Chain %s: E-FMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('E-FMI indicated no pathological behavior')
    if (quiet) return(TRUE)
  } else {
    if (!quiet) print('  E-FMI below 0.2 indicates you may need to reparameterize your model')
    if (quiet) return(FALSE)
  }
}

# source: https://github.com/claudiofronterre/onco/blob/master/R/pairs_stan.R
pairs_stan <- function(chain, stan_model, pars) {
  energy <- as.matrix(sapply(get_sampler_params(stan_model, inc_warmup = F),
                             function(x) x[,"energy__"]))
  pars <- extract(stan_model, pars = pars, permuted = F)
  df <- data.frame(energy[,chain], pars[,chain,])
  names(df)[1] <- "energy"
  GGally::ggpairs(df, title = paste0("Chain", chain),
                  lower = list(continuous = GGally::wrap("points", alpha = 0.2)))
}
