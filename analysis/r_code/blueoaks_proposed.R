source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)

set.seed(1)

calib_start_lobos <- 1950

monthly_precip_lobos <- read_csv(file = here('data', 'prism_casestudy.csv'))



# do transformations to climate first
lobos_precip_df <- monthly_precip_lobos %>%
  filter(year < 2024) %>%
  mutate(log_precip = log(prcp + 0.001)) %>%
  select(
    year,
    month,
    log_precip,
    precip = prcp
  )



clim_df_wateryear <- process_climate(lobos_precip_df,
                                    clim_name = 'log_precip',
                                    anom = T)

# ordered jan - december
month_clim_lobos <- lobos_precip_df %>%
  mutate(month = as_factor(month)) %>%
  group_by(month) %>%
  summarise(
    clim_mean = mean(log_precip),
    clim_sd = sd(log_precip)
  )

month_means_lobos <- month_clim_lobos$clim_mean
month_sds_lobos <- month_clim_lobos$clim_sd


trw_dfs_lobos <- format_trw_lobos()

trw_dfs_mod_lobos <- process_trw(trw_dfs_lobos)

unif_prior <- rep(1, 12)
halfyear_prior <- c(10,10,rep(1, 6), rep(10, 4))
df_list_lobos <- stan_setup(clim_df_wateryear, trw_dfs_mod_lobos,
                               start_year_calib = calib_start_lobos,
                               w_param = unif_prior, month_means = month_means_lobos,
                            month_sds = month_sds_lobos)




fit_lobos_new <- stan(file = here("analysis", "stan_code", "proposed_model.stan"),
                  iter= 5000,
                  data = df_list_lobos,
                  init_r = 1,
                  seed = 123#,
                  # control = list(adapt_delta = 0.99,
                  #                max_treedepth = 12)
)


#############################################################

# model evaluation in the space of the learned weighted average

##############################################################

recon_start_lobos <- min(df_list_lobos$start_years)
yearly_fig_df_lobos <- get_yearly_pred(fit_lobos, clim_df_wateryear, recon_start_lobos)


lobos_start_year_obs <- min(clim_df_wateryear$year)
fig_recon(yearly_fig_df_lobos, lobos_start_year_obs, calib_start_lobos,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim",
          display_years = c(1800, 1955))[[2]] +
  labs(title = 'blue oaks',
       y = 'log precip anom (sd)'
       )




lobos_w <- fig_weight(fit_lobos)


#############################################################

# model evaluation in the space of the DJF average

##############################################################

w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)
# just for evaluation
djf_df_lobos <- process_climate(lobos_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf, anom = F)

djf_fig_lobos_df <- get_djf_pred(fit_lobos, djf_df_lobos, recon_start_lobos)

fig_recon(djf_fig_lobos_df,
          lobos_start_year_obs, calib_start_lobos,
          pred_name = 'mean',
          display_years = c(1750, 1950)
)[[2]] +
  labs(title = 'blue oaks',
       subtitle = '',
       y = 'DJF mean precip (mm/month)')




