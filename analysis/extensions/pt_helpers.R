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


#' a set of "trees" as a linear function of PRECIP
#' @param sim_df is a dataset which has columns temp, precip, year, month
#' @param w_true are the 12 monthly weights that the tree will be sensitive to
#' @return a list containing 2 dataframes: the first with monthly growth, the second with simulated tree rings (yearly)
make_tree_p_t_linear <- function(sim_df, n_trees = NULL,
                                  w_p_true = rep(1/12, 12), w_t_true = rep(1/12, 12), sigma2_w = 1/2,
                                  beta0 = rnorm(1), betat = rnorm(1), betap = rnorm(1), beta_age = runif(1), take_log_p = F){

  if(take_log_p){
    sim1_trw <- sim_df %>%
      group_by(year) %>%
      summarize(
        precip_yearw = log(sum(w_p_true * precip) + 0.001),
        temp_yearw = sum(w_t_true * temp),
        trw_sim = beta0 + betap * precip_yearw + betat * temp_yearw,
        .groups = "drop")

  } else{
    sim1_trw <- sim_df %>%
      group_by(year) %>%
      summarize(
        precip_yearw = sum(w_p_true * precip),
        temp_yearw = sum(w_t_true * temp),
        trw_sim = beta0 + betap * precip_yearw + betat * temp_yearw,
        .groups = "drop")
  }


  sim1_df <- sim1_trw


  n_trees <- n_trees
  sim_trws <- map_dfr(1:n_trees, ~ tibble(
    tree_id = .x,
    year = sim1_df$year,
    age = 1:nrow(sim1_df),
    precip = sim1_df$precip_yearw,
    temp = sim1_df$temp_yearw
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
#' @param clim_df output of `process_climate()`, joined temp and precip
#' @param n_trees positive integer, number of pseudoproxy trees to create
#' @param snr positive number (typically less than 2). Signal to noise ratio of the pseudoproxies
#' @param w_true vector of length 12 which sums to 1. Defines the growing season of the pseudoproxies.
#' @param beta0 intercept of the trw ~ climate relationship
#' @param beta1 slope of the trw ~ climate relationship
#' @returns: list of dataframes (length n_trees) with columns `tree_id`, `year`, `trw`, and `age`.
#' @example format_trw_linear(sim_precip_df, 10, 1, 1500, rep(1/12, 12))
format_trw_p_t_linear <- function(clim_df, n_trees, snr, start_year_calib,
                              w_p_true,  w_t_true, beta0 = 1, betat = 2, betap = 2, beta_age = 0, include_growth = F, take_log = F){

  ### generative part: pick a first year of data for each tree
  first_year_obs <- sample(min(clim_df$year):start_year_calib, n_trees, replace = T)

  ## we can also pick out a last year of obs for the trees, if we don't want them all to live to current day.
  # last_year_obs <- sample((min(sim_df$year) + 2):max(sim_df$year))

  start_year_reconstruction <- min(first_year_obs)


  noise <- 1 / (1 + snr^2)
  fake_trees_full <- make_tree_p_t_linear(clim_df, w_p_true = w_p_true, w_t_true = w_t_true,
                                           n_trees = n_trees, sigma2_w = noise,
                                           beta0 = beta0, betat = betat, betap = betap, beta_age = beta_age,
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


stan_setup_p_t <- function(clim_df, trw_dfs,
                       start_year_calib, start_year_reconstruct = NULL,
                       w_param_precip = rep(1, 12), w_param_temp = rep(1, 12),
                       month_model = T, month_means_p = NULL, month_sds_p = NULL,
                       month_means_t = NULL, month_sds_t = NULL
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
    )


  precip_obs <- clim_obs$precip
  temp_obs <- clim_obs$temp


  # for evaluation, get monthly validation climate
  # e.g. vector of earliest tree year -> first instrumental observed year
  clim_unobs <- clim_df %>%
    filter(year < start_year_calib,
           year >= start_year_reconstruct)

  precip_unobs <- clim_unobs$precip
  temp_unobs <- clim_unobs$temp

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
      N_years_obs_precip = length(precip_obs) / 12,
      N_years_obs_temp = length(temp_obs) / 12,
      N_trees = n_trees,
      trw = trw,
      age = age,
      precip_obs = precip_obs,
      temp_obs = temp_obs,
      s = trw_lengths,
      start_years = first_years_obs,
      w_param_precip = w_param_precip,
      w_param_temp = w_param_temp,
      w_djf = c(1/3, 1/3, rep(0, 9), 1/3),
      w_jja = c(rep(0, 5), 1, 1, 1, rep(0, 4)) / 3,
      month_means_mod_t = month_means_t,
      month_sds_mod_t = month_sds_t,
      month_means_mod_p = month_means_p,
      month_sds_mod_p = month_sds_p
    )
  } else{
    df_fit <- list(
      N_years = max(last_years_obs) - min(first_years_obs) + 1,
      N_years_obs_precip = length(precip_obs),
      N_years_obs_temp = length(temp_obs),
      N_trees = n_trees,
      trw = trw,
      age = age,
      precip_obs = precip_obs,
      temp_obs = temp_obs,
      s = trw_lengths,
      start_years = first_years_obs
    )
  }



  df_fit
}





#' produce a barchart of the estimated seasonal window
#' one bar per month, height is posterior mean
#' @param fit is the output of stan(), with parameters w[1], ..., w[12]
fig_weight_pt <- function(fit){

  # extract weight posterior draws
  precip_weights <- fit %>%
    spread_draws(w_p[i]) %>%
    summarise_draws()
  temp_weights <- fit %>%
    spread_draws(w_t[i]) %>%
    summarise_draws()

  # plot weights
  fig_p <- precip_weights %>%
    mutate(i = factor(i, levels = c(10, 11, 12, 1:9))) %>% # water year ordering
    ggplot(aes(i, mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +   # 1/12 represents uninformative "uniform" weights
    labs(x = 'month', y = 'weights', title = 'precip')

  fig_t <- temp_weights %>%
    mutate(i = factor(i, levels = c(10, 11, 12, 1:9))) %>% # water year ordering
    ggplot(aes(i, mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +   # 1/12 represents uninformative "uniform" weights
    labs(x = 'month', y = 'weights', title = 'temp')


  fig_p / fig_t
}


#' function to extract the weighted yearly reconstruction from a stanfit object.
#' @param fit the output of a stan model with parameter 'clim_w_yearly'
#' @param clim_df a dataframe that contains instrumental monthly climate for a single location (one row per month-year)
#' @param start_year_reconstruct an integer, the year that we want to start producing a reconstruction
#' @param swap_log a logical flag for whether we take log(weighted climate average) or weighted log-climate average
get_yearly_pred_pt <- function(fit, clim_df, start_year_reconstruct, swap_log = F){

  clim_mis <- fit %>%
    spread_draws(clim_w_yearly[i]) %>%
    summarise_draws()


  mis_year_df <- clim_mis %>%
    ungroup() %>%
    mutate(year = start_year_reconstruct + i - 1,
           clim_mean_pred = mean
    )


  ## better method of the "truth": comparing weighted avg draw by draw --------

  month_lookup <- tibble(
    month = month.abb,
    month_no = 1:12
  )


  true_climate <- clim_df %>%
    # filter(between(year ,start_year_reconstruct, calib_start - 1)) %>%
    left_join(month_lookup, by = 'month')



  w_draws <- fit %>%
    spread_draws(w[i])

  if(swap_log){
    w_trueclim_draws <- w_draws %>%
      group_by(.chain, .iteration) %>%
      group_map(~ {
        left_join(true_climate, .x, by = c('month_no' = 'i')) %>%
          group_by(year) %>%
          summarise(water_year_w_clim = log(sum(clim * w)))
      }
      ) %>%
      bind_rows(.id = 'draw')

  } else{
    w_trueclim_draws <- w_draws %>%
      group_by(.chain, .iteration) %>%
      group_map(~ {
        left_join(true_climate, .x, by = c('month_no' = 'i')) %>%
          group_by(year) %>%
          summarise(water_year_w_clim = sum(clim * w))
      }
      ) %>%
      bind_rows(.id = 'draw')
  }



  w_trueclim_year <- w_trueclim_draws %>%
    group_by(year) %>%
    summarize(mean_w_trueclim = mean(water_year_w_clim),
              q5_w_clim = quantile(water_year_w_clim, .05),
              q95_w_clim = quantile(water_year_w_clim, .95)
    )


  wyear_df <- mis_year_df %>%
    left_join(w_trueclim_year, by = 'year')

  wyear_df
}


