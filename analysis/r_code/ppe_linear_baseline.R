################################################################

### Pseudoproxy experiments using the "baseline model", illustrating the impact of misspecified seasonal windows.
### (Figure 1 in the manuscript)

###############################################################

source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))
options(mc.cores = parallelly::availableCores())

start_year_reconstruct <- 1600
start_year_calib <- 1950

# monthly climate model output (first run CESM1 last mil)
sim_df <- read_csv(here('data', 'sim_cesm1_mil_data.csv')) %>%
  filter(year >= start_year_reconstruct) %>%
  filter(year != 2004) # 2004 only has data up to sept, so remove for simplicity to have only years of full data


#' function to run a pseudoproxy experiment
#' @param w_true
#' @param w_guess
#' @param snr
#' @param calib_start_ppe
season_baseline <- function(sim_df_raw, w_true = c(1/3, 1/3, rep(0, 9), 1/3), w_guess = c(1/3, 1/3, rep(0, 9), 1/3),
                            snr=0.75, calib_start_ppe = 1950){
  # do transformations to climate first
  sim_precip_df <- sim_df_raw %>%
    mutate(
      log_precip = log(precip_mmmonth + 0.001),
      month = str_to_title(name)) %>%
    select(
      year,
      month,
      log_precip,
      precip = precip_mmmonth
    )

  # aggregate to get our reconstruction target
  clim_df_target_ppe <- process_climate(sim_precip_df,
                                     clim_name = 'precip',
                                     w_agg = w_guess, anom = F)


  # take log at yearly level
  clim_df_target_ppe <- clim_df_target_ppe %>%
    mutate(
      log_clim = log(clim + 0.001)
    )

  # pull out mean and standard deviation
  # so we can standardize only using the training period
  logclim_training <- clim_df_target_ppe %>%
    filter(year >= start_year_calib) %>%
    pull(log_clim)

  mean_logtarget = mean(logclim_training)
  sd_logtarget = sd(logclim_training)

  # standardize
  clim_df_target_ppe_std <- clim_df_target_ppe %>%
    mutate(
      clim = (log_clim - mean_logtarget) / sd_logtarget
    ) %>%
    select(-log_clim)

  # keep monthly to make trees
  clim_df_ppe <- process_climate(sim_precip_df,
                                 clim_name = 'precip', anom = F)

  # make the pseudoproxies
  trw_dfs_ppe <- format_trw_linear(clim_df_ppe,
                                   n_trees = 10,
                                   snr = snr,
                                   start_year_calib = calib_start_ppe,
                                   w_true = w_true,
                                   take_log = T,
                                   beta1 = 2, beta0 = 1, beta_age = -0.01
  )

  # format the pseudoproxies into a common format
  trw_dfs_mod_ppe <- process_trw(trw_dfs_ppe)

  # combine all data together in a nice list
  df_list_ppe_baseline <- stan_setup(clim_df_target_ppe_std, trw_dfs_mod_ppe,
                                     start_year_calib = calib_start_ppe, month_model = F)
  df_list_ppe_baseline$mean_logdjf = mean_logtarget
  df_list_ppe_baseline$sd_logdjf = sd_logtarget


  # fit baseline model (does not include model discrepancy term by default)
  fit_ppe_baseline <- stan(file = here("analysis", "stan_code", "baseline_model.stan"),
                           iter= 5000,
                           data = df_list_ppe_baseline,
                           seed = 123,
                           init_r = 1#,
                           # control = list(adapt_delta = 0.99)
                           )

  # reconstruction starts at the first year of pseudoproxy data
  recon_start_ppe_baseline <- min(df_list_ppe_baseline$start_years)

  # pull out reconstruction in standardized log space
  clim_draws_ppe_baseline <- fit_ppe_baseline %>%
    spread_draws(clim_yearly[i])
  clim_mis_ppe_baseline <- clim_draws_ppe_baseline%>%
    summarise_draws()

  # shared df with reconstruction and instrumental climate (one row per year)
  baseline_ppe_df <- clim_mis_ppe_baseline %>%
    ungroup() %>%
    mutate(year = recon_start_ppe_baseline + i - 1,
           clim_mean_pred = mean
    ) %>%
    left_join(clim_df_target_ppe_std, by = 'year')

  # pull out reconstruction in original units
  clim_mis_ppe_baseline_raw <- fit_ppe_baseline %>%
    spread_draws(clim_raw[i]) %>%
    summarise_draws()

  # shared df with reconstruction and instrumental climate
  baseline_raw_ppe_df <- clim_mis_ppe_baseline_raw %>%
    ungroup() %>%
    mutate(year = recon_start_ppe_baseline + i - 1,
           clim_mean_pred = mean
    ) %>%
    left_join(clim_df_target_ppe, by = 'year')


  # return
  list(
    fit = fit_ppe_baseline,
    ppe_df = baseline_ppe_df,
    ppe_raw_df = baseline_raw_ppe_df)
}

# truth is DJF

w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)
w_jja <- c(rep(0, 5), 1/3, 1/3, 1/3, rep(0, 4))
w_year <- rep(1/12, 12)
w_j<- c(1,rep(0, 11))


# we target DJF, and the pseudoproxies are sensitive to DJF
exact_target <- season_baseline(sim_df,
                                w_true = w_djf,
                                w_guess = w_djf,
                                calib_start_ppe = start_year_calib)

# we target the full year, and the pseudoproxies are sensitive to DJF
big_target <- season_baseline(sim_df,
                              w_true = w_djf,
                              w_guess = w_year,
                              calib_start_ppe = start_year_calib)

# we target July, and the pseudoproxies are sensitive to DJF
small_target <- season_baseline(sim_df,
                                w_true = w_djf,
                                w_guess = w_j,
                                calib_start_ppe = start_year_calib)

# we target JJA, and the pseudoproxies are sensitive to DJF
wrong_target <- season_baseline(sim_df,
                                w_true = w_djf,
                                w_guess = w_jja,
                                calib_start_ppe = start_year_calib)



baseline_list <- list(
  Exact =  exact_target,
  Big = big_target,
  Small = small_target,
  Wrong = wrong_target
)

# saveRDS(baseline_list, "baseline_list_eta.Rds")


# plot reconstructions for each in the standardized space
fig_baseline <- baseline_list %>%
  imap(~.x$ppe_df %>%
    fig_recon(start_year_calib = start_year_calib,
              start_year_obs = start_year_reconstruct, fig_only = T,
              display_years = c(1750, 1950)) +
      labs(title = .y,
           y = 'Standardized log precipitation')
  )



# combine all plots together
wrap_plots(fig_baseline) +
  plot_layout(axis_titles = 'collect') &
  scale_y_continuous(expand = c(0.2, .2))




