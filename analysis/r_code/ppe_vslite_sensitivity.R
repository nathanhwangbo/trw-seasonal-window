source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)


set.seed(1)

calib_start_ppe <- 1950
w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)

# monthly climate model output (first run CESM1 last mil)
sim_df <- read_csv(here('data', 'sim_cesm1_mil_data.csv')) %>%
  filter(year >= start_year_reconstruct) %>%
  filter(year != 2004) # 2004 only has data up to sept, so remove for simplicity to have only years of full data




#' function to test limits of weight parameter learning
#' in the case where the model is correctly specified
#' i.e. the trees are truly generated as a weighted linear combination of pdsi
#' this is basically a function-ized version of 1_analysis_ppe.R
#' @param snr signal-to-noise ratio of the tree. values between .25 and 1 are in a reasonable-ish range for proxies
#' @param w_true "true" seasonality of the trees, as weight vector
#' @return list of two items: (1) rstan model fit (2) "entropy" of the true weight vector (w_true) (3) kl divergence of the weight param prior vs posterior
sensitivity_vslite <- function(snr, sim_df_raw, suffix, w_true, weight_param = rep(1, 12), n_trees = 10,
                                  calib_start_ppe = 1800){

  message(str_glue("snr={snr}"))
  set.seed(1)
  sim_precip_df <- sim_df_raw %>%
    mutate(log_precip = log(precip_mmmonth + 0.001),
           month = str_to_title(name)) %>%
    select(
      year,
      month,
      log_precip,
      precip = precip_mmmonth
    )


  # ordered jan - december
  month_clim_ppe <- sim_precip_df %>%
    mutate(month = as_factor(month)) %>%
    group_by(month) %>%
    summarise(
      clim_mean = mean(log_precip),
      clim_sd = sd(log_precip)
    )

  month_means_ppe <- month_clim_ppe %>%
    pull(clim_mean)
  month_sds_ppe <- month_clim_ppe %>%
    pull(clim_sd)



  clim_df <- process_climate(sim_precip_df,
                             clim_name = 'log_precip',
                             anom = T)

  # process_climate is just used here to convert to water year
  vslite_clim_df <- sim_df_raw %>%
    mutate(year = year(time_new), month = month(time_new, label = T)) %>%
    select(year, month,
           precip_mmmonth, temp_c, lat) %>%
    process_climate(clim_name = '', anom = F, standardize_names = F) %>%
    select(-year) %>%
    rename(year = water_year)


  trw_dfs_ppe <- format_trw_vslite(vslite_clim_df,
                                   n_trees = n_trees,
                                   snr = snr,
                                   start_year_calib = calib_start_ppe,
                                   month_nums = which(w_true != 0),
                                   include_extras = T)

  trw_mod_dfs_ppe <- process_trw(trw_dfs_ppe$trw)

  df_list_ppe <- stan_setup(clim_df, trw_mod_dfs_ppe,
                            start_year_calib = calib_start_ppe,
                            w_param = weight_param, month_means = month_means_ppe,
                            month_sds = month_sds_ppe)

  fit <- stan(file = here("analysis", "stan_code", "proposed_model.stan"),
              iter= 1000,
              data = df_list_ppe)


  recon_start_ppe <- min(df_list_ppe$start_years)
  yearly_fig_ppe_df <- get_yearly_pred(fit, clim_df, recon_start_ppe)


  post_w <- fit %>%
    spread_draws(w[i]) %>%
    summarise_draws() %>%
    ungroup()

  post_w_params <- fit %>%
    spread_draws(rho_w, alpha_w, sigma_w) %>%
    summarise_draws()

  # this is assuming djf is w_true
  djf_df_ppe <- process_climate(sim_precip_df,
                                clim_name = 'precip',
                                w_agg = w_true, anom = F)

  djf_fig_ppe_df <- get_djf_pred(fit, djf_df_ppe, recon_start_ppe)


  list(
    yearly_fig_ppe_df = yearly_fig_ppe_df,
    djf_fig_ppe_df = djf_fig_ppe_df,
    post_w = post_w,
    post_w_params = post_w_params,
    weight_param = weight_param,
    trw_df = trw_mod_dfs_ppe
  )
}

snr_grid <- c(.25, .5, .75)

# true sensitivities
w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)


# priors
unif_prior <- rep(1, 12)
djf10_prior <- c(10,10,0,0,0,0,0,0,0,0,0,10)
# sondjf
halfyear10_prior <- c(10,10,rep(1, 6), rep(10, 4))



truedjf_unif_vslite <- map(snr_grid, ~sensitivity_vslite(snr = .x, sim_df = sim_df, suffix = "lobosppe_unif1prior_djftruth", weight_param = unif_prior, w_true = w_djf))
truedjf_djf10_vslite <- map(snr_grid, ~sensitivity_vslite(snr = .x, sim_df = sim_df, suffix = "lobosppe_djf10prior_djftruth", weight_param = djf10_prior, w_true = w_djf))
truedjf_halfyear10_vslite <- map(snr_grid, ~sensitivity_vslite(snr = .x, sim_df = sim_df, suffix = "lobosppe_djf10prior_djftruth", weight_param = halfyear10_prior, w_true = w_djf))

vslite_testlist <- list(
  truedjf_unif = truedjf_unif_vslite,
  truedjf_djf10 = truedjf_djf10_vslite,
  truedjf_halfyear10 = truedjf_halfyear10_vslite
)

saveRDS(vslite_testlist, here('temp_output', 'vslite_sensitivity', str_glue('vslite_testlist.Rds')))



(unif_fig_vslite <- truedjf_unif_vslite %>%
    map_dfr(~.x$post_w, .id = 'SNR') %>%
    mutate(snr = case_when(SNR == 1 ~ 0.25,
                           SNR == 2 ~ 0.5,
                           SNR == 3 ~ 0.75)) %>%
    ggplot(aes(factor(i), mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +
    labs(
      x = 'month',
      y = 'weight',
      title = 'unif(1) prior'
    ) +
    facet_wrap(~snr, nrow = 1)
)


(halfyear10_fig_vslite <- truedjf_halfyear10_vslite %>%
    map_dfr(~.x$post_w, .id = 'SNR') %>%
    mutate(snr = case_when(SNR == 1 ~ 0.25,
                           SNR == 2 ~ 0.5,
                           SNR == 3 ~ 0.75)) %>%
    ggplot(aes(factor(i), mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +
    labs(
      x = 'month',
      y = 'weight',
      title = 'halfyear(10) prior'
    ) +
    facet_wrap(~snr, nrow = 1)
)

(djf10_fig_vslite <- truedjf_djf10_vslite %>%
    map_dfr(~.x$post_w, .id = 'SNR') %>%
    mutate(snr = case_when(SNR == 1 ~ 0.25,
                           SNR == 2 ~ 0.5,
                           SNR == 3 ~ 0.75)) %>%
    ggplot(aes(factor(i), mean)) +
    geom_col() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1/12) +
    labs(
      x = 'month',
      y = 'weight',
      title = 'DJF(10) prior'
    ) +
    facet_wrap(~snr, nrow = 1)
)
