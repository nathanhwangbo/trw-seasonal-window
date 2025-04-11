#####################################

# functions that help us make figures

#####################################



#' produce a barchart of the estimated seasonal window
#' one bar per month, height is posterior mean
#' @param fit is the output of stan(), with parameters w[1], ..., w[12]
fig_weight <- function(fit){

  # extract weight posterior draws
  clim_weights <- fit %>%
    spread_draws(w[i]) %>%
    summarise_draws()

  # plot weights
  clim_weights %>%
    mutate(i = factor(i, levels = c(10, 11, 12, 1:9))) %>% # water year ordering
    ggplot(aes(i, mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +   # 1/12 represents uninformative "uniform" weights
    labs(x = 'month', y = 'weights')
}




#' produce a plot of the time series reconstruction
#' @param wyear_df is a dataframe with one year per row, with at least five columns:
#'        one column called "year",
#'        one column with the mean reconstructed values (name specified in "pred_name"),
#'        one column with instrumental values (name specified in "truth_name"),
#'        one column with some quantile of reconstructed values (name specified in qlow_name),
#'        one column with some quantile of reconstructed values (name specified in qhigh_name)
#' @param start_year_obs an integer, the first year of instrumental climate data
#' @param start_year_calib an integer, the first year used for calibration
#' @param pred_name a string, the name of the column in wyear_df that contains mean reconstruction values
#' @param truth_name a string, the name of the column in wyear_df that contains instrumental climate values
#' @param qlow_name a string, the name of the column in wyear_df
#' @param display_years a numeric vector of length 2, the first and last year to be shown in the plot
#' @fig_only a logical flag for whether summary statistics should be returned alongside the ggplot object
fig_recon <- function(wyear_df, start_year_obs, start_year_calib,
                      pred_name = "clim_mean_pred", truth_name = "clim",
                      qlow_name = 'q5', qhigh_name = 'q95',
                      display_years = NULL, fig_only = F) {

  pred_var <- sym(pred_name)
  truth_var <- sym(truth_name)
  qlow_var <- sym(qlow_name)
  qhigh_var <- sym(qhigh_name)

  if(is.null(display_years)){
    display_years <- c(start_year_calib - 200, start_year_calib + 10)
  }

  ann_text <- wyear_df %>%
    filter(year < start_year_calib,
           year > start_year_obs) %>%
    summarize(
      cor_mean = cor(!!truth_var, !!pred_var),
      coverage9 = mean(between(!!truth_var, !!qlow_var, !!qhigh_var)),
      var_truth = var(!!truth_var),
      var_mean = var(!!pred_var)
    ) %>%
    round(2) %>%
    mutate(across(cor_mean:coverage9, ~format(round(.x, 2), nsmall = 2))) %>%
    mutate(
      label1 = str_glue('R: {cor_mean}'),
      label2 = str_glue('90% Coverage: {coverage9}')
    )


  fig <- wyear_df %>%
    filter(year > min(year)) %>%
    ggplot() +
    geom_ribbon(aes(x = year, ymin = !!qlow_var, ymax = !!qhigh_var), fill = 'plum', alpha = 0.7) +
    geom_line(aes(year, !!pred_var), color = 'purple3', size = 1) +
    geom_line(aes(year, !!truth_var), size = 1, alpha = 0.7) +
    lims(
      x = display_years
    ) +
    geom_text(data = ann_text,
              mapping = aes(x = -Inf,
                            y = Inf, label = label1),
              hjust = -0.05,
              vjust = 1.11,
              size = 24/.pt
    ) +
    geom_text(data = ann_text,
              mapping = aes(x = Inf,
                            y = Inf, label = label2),
              hjust = 1.01,
              vjust = 1.11,
              size = 24/.pt
    )

  if(fig_only){
    return(fig)
  }

  list(
    ann_text = ann_text,
    recon = fig
  )

}

#' function to extract a djf reconstruction from a stanfit object.
#' @param fit the output of a stan model with parameter 'clim_djf'
#' @param djf_df a dataframe that contains instrumental djf climate for a single location (one row per year)
#' @param start_year_reconstruct an integer, the year that we want to start producing a reconstruction
#' @param start_year_calib an integer, the first year we use for calibration.
get_djf_pred <- function(fit, djf_df, start_year_reconstruct, start_year_calib = 1950){

  djf_draws <- fit %>%
    spread_draws(clim_djf[i])

  clim_mis <- djf_draws %>%
    summarise_draws()


  mis_year_df <- clim_mis %>%
    ungroup() %>%
    mutate(year = start_year_reconstruct + i - 1,
           clim_mean_pred = mean
    )



  fig_df <- mis_year_df %>%
    left_join(djf_df, by = 'year')

  ## get CRPS


  set.seed(123)


  start_year_obs <- min(start_year_reconstruct, min(djf_df$year))

  # see which i correspond to reconstruction years
  validation_i <- mis_year_df %>%
    filter(year < start_year_calib,
           year > start_year_obs) %>%
    pull(i)

  # take 2 samples of size 1000 (can modify)
  n_samples <- 1000
  n_samples_tot <- 2 * n_samples
  rand_order <- sample(1:n_samples_tot, size = n_samples_tot, replace = F)
  s1 <- rand_order[1:n_samples]
  s2 <- setdiff(rand_order, s1)

  draws1 <- djf_draws %>%
    ungroup() %>%
    select(i, contains('clim')) %>%
    filter(i %in% validation_i) %>%
    group_by(i) %>%
    slice(s1) %>%
    group_split() %>%
    map(~deframe(.x)) %>%
    bind_cols(.name_repair = 'universal_quiet') %>%
    as.matrix.data.frame()

  draws2 <- djf_draws %>%
    ungroup() %>%
    select(i, contains('clim')) %>%
    filter(i %in% validation_i) %>%
    group_by(i) %>%
    slice(s2) %>%
    group_split() %>%
    map(~deframe(.x)) %>%
    bind_cols(.name_repair = 'universal_quiet') %>%
    as.matrix.data.frame()

  # assuming clim is the name of the true target
  truth_vec <- fig_df %>%
    filter(year < start_year_calib,
           year > start_year_obs) %>%
    pull(clim)

  yr_crps <- loo::crps(draws1, draws2, truth_vec)


  # combine it all together
  fig_df <- fig_df %>%
    mutate(crps = yr_crps$estimates[1])

  fig_df
}

#' function to extract the weighted yearly reconstruction from a stanfit object.
#' @param fit the output of a stan model with parameter 'clim_w_yearly'
#' @param clim_df a dataframe that contains instrumental monthly climate for a single location (one row per month-year)
#' @param start_year_reconstruct an integer, the year that we want to start producing a reconstruction
#' @param swap_log a logical flag for whether we take log(weighted climate average) or weighted log-climate average
get_yearly_pred <- function(fit, clim_df, start_year_reconstruct, swap_log = F){

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












