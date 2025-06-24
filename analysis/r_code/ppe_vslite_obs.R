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
  mutate(
         log_precip = log(precip_mmmonth + 0.001),
         month = str_to_title(name)) %>%
  select(
    year,
    month,
    log_precip,
    precip = precip_mmmonth
  )

# format instrumental climate data
clim_df_ppe <- process_climate(sim_precip_df,
                                    clim_name = 'log_precip')

# format instrumental climate data in original units
# process_climate is just used here to convert to water year
vslite_clim_df <- sim_df %>%
  mutate(month = str_to_title(name)) %>%
  select(year, month,
         precip_mmmonth, temp_c, lat) %>%
  process_climate(clim_name = '', anom = F, standardize_names = F)

# create truncated VS-lite pseudoproxies -- i.e. modified to only grow in DJF
snr_ppe <- 0.75
# breitenmoser et al table 2 means
b_t1 <- 4.64
b_t2 <- 16.34
# b_t2 <- 15.34 # pine
b_m1 <- 0.023
b_m2 <- 0.44
# b_m2 <- 0.5 # pine
trw_dfs_vslite2 <- format_trw_vslite(vslite_clim_df,
                                 n_trees = 20,
                                 snr = snr_ppe,
                                 start_year_calib = calib_start_ppe,
                                 month_nums = 1:12,
                                 # month_nums = c(12,1,2),
                                 include_extras = T,
                                 sm_limited = F,
                                 M1 = b_m1, M2 = b_m2,
                                 T1 = b_t1, T2 = b_t2)


# saveRDS(trw_dfs_vslite2, 'E:/Projects/old-blue-oaks/temp_output/trw_dfs_vslite2.Rds')
trw_dfs_vslite2 <- readRDS('E:/Projects/old-blue-oaks/temp_output/trw_dfs_vslite2.Rds')


# instrumental mean and sd for each month.
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

# format the trw data.
trw_mod_dfs_ppe <- process_trw(trw_dfs_vslite2$trw)

# format data to feed into stan
df_list_ppe <- stan_setup(clim_df_ppe, trw_mod_dfs_ppe,
                               start_year_calib = calib_start_ppe,
                          month_means = month_means_ppe,
                          month_sds = month_sds_ppe)


# fit model
fit_ppe_vslite <- stan(file = here("analysis",
                                   "stan_code",
                                   "proposed_model_precip.stan"),
     iter= 2000,
     data = df_list_ppe,
     init_r = 1,
     seed = 1234#,
     # control = list(adapt_delta = 0.9, max_treedepth = 12)
)

b= # saveRDS(fit_ppe_vslite, "E:/Projects/old-blue-oaks/temp_output/fit_ppe_vslite2.Rds")
fit_ppe_vslite <- readRDS("E:/Projects/old-blue-oaks/temp_output/fit_ppe_vslite2.Rds")

##################

# model evaluation

##################

# unless otherwise stated, the reconstruction starts with the oldest pseudoproxy.
recon_start_ppe <- min(df_list_ppe$start_years)

# plot weights (e.g. what was learned from the models)
fig_w <- fig_weight(fit_ppe_vslite)


# get weighted yearly average reconstruction
yearly_fig_ppe_df <- get_yearly_pred(fit_ppe_vslite, clim_df_ppe, recon_start_ppe)
yearly_fig_ppe_df <- yearly_fig_ppe_df %>%
  mutate(
    recon_start = recon_start_ppe,
    calib_start = calib_start_ppe
  )
# write_csv(yearly_fig_ppe_df, 'E:/Projects/old-blue-oaks/temp_output/yearly_fig_df_vslite2.csv')


# plot weighted yearly reconstruction
fig_recon(yearly_fig_ppe_df,
          recon_start_ppe, calib_start_ppe,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim",
          display_years = c(1750, 1950),
          fig_only = T
) +
  ## below: optionally show variability in weights
  # geom_ribbon(aes(x = year, ymin = q5_w_clim, ymax = q95_w_clim), alpha = 0.5) +
  labs(title = str_glue('vslite ppe, log, snr = {snr_ppe}'),
       y = 'w transformed anom precip')




# get instrumental djf climate
djf_df_ppe <- process_climate(sim_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf)

# get djf reconstruction from model
djf_fig_ppe_df <- get_djf_pred(fit_ppe_vslite, djf_df_ppe, recon_start_ppe)

djf_fig_ppe_df <- djf_fig_ppe_df %>%
  mutate(
    recon_start = recon_start_ppe,
    calib_start = calib_start_ppe
  )

# write_csv(djf_fig_ppe_df, 'E:/Projects/old-blue-oaks/temp_output/djf_fig_df_vslite2.csv')


# plot djf reconstruction
(fig_djf_recon <- fig_recon(djf_fig_ppe_df,
          recon_start_ppe, calib_start_ppe,
          display_years = c(1750, 1950), fig_only = T
) +
  labs(
    title = 'VS-Lite',
       y = 'DJF mean precip (mm/month)')
)

(fig_vslite <- fig_djf_recon + fig_w +
  plot_layout(widths = c(3,1)))


