source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)

set.seed(1)

calib_start_ppe <- 1950
start_year_reconstruct <- 1600
w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)


# monthly climate model output (first run CESM1 last mil)
sim_df <- read_csv(here('data', 'sim_cesm1_mil_data.csv')) %>%
  filter(year >= start_year_reconstruct) %>%
  filter(year != 2004) # 2004 only has data up to sept, so remove for simplicity to have only years of full data


# do transformations to climate first
sim_precip_df <- sim_df %>%
  mutate(log_precip = log(precip_mmmonth + 0.001),
         month = str_to_title(name)) %>%
  select(
    year,
    month,
    log_precip,
    precip = precip_mmmonth
  )

clim_df_ppe <- process_climate(sim_precip_df,
                                    clim_name = 'log_precip')

# process_climate is just used here to convert to water year
vslite_clim_df <- sim_df %>%
  mutate(month = str_to_title(name)) %>%
  select(year, month,
         precip_mmmonth, temp_c, lat) %>%
  process_climate(clim_name = '', anom = F, standardize_names = F) #%>%
  # select(-year) #%>%
  # rename(year = water_year)

snr_ppe <- 0.75

trw_dfs_vslite <- format_trw_vslite(vslite_clim_df,
                                 n_trees = 10,
                                 snr = snr_ppe,
                                 start_year_calib = calib_start_ppe,
                                 month_nums = 1:12,
                                 include_extras = T,
                                 sm_limited = T)


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

trw_mod_dfs_ppe <- process_trw(trw_dfs_vslite$trw)


df_list_ppe <- stan_setup(clim_df_ppe, trw_mod_dfs_ppe,
                               start_year_calib = calib_start_ppe,
                          month_means = month_means_ppe,
                          month_sds = month_sds_ppe)




# this will take a while!
fit_ppe_vslite <- stan(file = here("analysis", "stan_code",
                                   "proposed_model_precip.stan"),
     iter= 1000,
     data = df_list_ppe,
     init_r = 1,
     seed = 123#,
     #control = list(adapt_delta = 0.9, max_treedepth = 12)
)




##################

# model evaluation

##################

## in weighted standardized anomaly space ----------------

recon_start_ppe <- min(df_list_ppe$start_years)
yearly_fig_ppe_df <- get_yearly_pred(fit_ppe_vslite, clim_df_ppe, recon_start_ppe)
yearly_fig_ppe_df <- yearly_fig_ppe_df %>%
  mutate(
    recon_start = recon_start_ppe,
    calib_start = calib_start_ppe
  )
# write_csv(yearly_fig_ppe_df, here('paper_output', 'yearly_fig_df_vslite_full.csv'))


fig_recon(yearly_fig_ppe_df,
          recon_start_ppe, calib_start_ppe,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim",
          fig_only = T
) +
  geom_ribbon(aes(x = year, ymin = q5_w_clim, ymax = q95_w_clim), alpha = 0.5) +
  labs(title = str_glue('vslite ppe, log, snr = {snr_ppe}'),
       y = 'w transformed anom precip')


## in djf mm/month precip space ----------------

# just for evaluation
djf_df_ppe <- process_climate(sim_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf)

djf_fig_ppe_df <- get_djf_pred(fit_ppe_vslite, djf_df_ppe, recon_start_ppe)

djf_fig_ppe_df <- djf_fig_ppe_df %>%
  mutate(
    recon_start = recon_start_ppe,
    calib_start = calib_start_ppe
  )


fig_recon(djf_fig_ppe_df,
          recon_start_ppe, calib_start_ppe,
          display_years = c(1750, 1950),
          fig_only = T
) +
  labs(
    title = '',
       y = 'DJF mean precip (mm/month)')

fig_weight(fit_ppe_vslite) +
  labs(title = 'Posterior Weights',
       subtitle = '')

# comparison to monthly correlations
trw_dfs_vslite$cor_by_month %>%
  mutate(month_no = factor(month_no, levels = c(10, 11, 12, 1:9))) %>%
  ggplot() +
  geom_col(aes(month_no, cor_mon)) +
  labs(
    title = "Monthly Correlation",
    x = 'month',
    y = 'cor(Precip, TRW)'
  )

