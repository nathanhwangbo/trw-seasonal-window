source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

nc <- parallelly::availableCores()
options(mc.cores = nc)


set.seed(1)

w_djf <- c(1/3, 1/3, rep(0, 9), 1/3)
w_jja <- c(rep(0, 5), 1/3, 1/3,  1/3, rep(0, 4))

start_year_reconstruct <- 1700
calib_start_ppe <- 1950

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

# do transformations to climate first
sim_temp_df <- sim_df %>%
  mutate(month = str_to_title(name)) %>%
  select(
    year,
    month,
    temp = temp_c
  )


precip_df_ppe <- process_climate(sim_precip_df,
                               clim_name = 'log_precip',
                               anom = T) %>%
  rename('precip' = clim)

temp_df_ppe <- process_climate(sim_temp_df,
                               clim_name = 'temp',
                               anom = T) %>%
  rename('temp' = clim)



clim_df_ppe <-
  inner_join(precip_df_ppe, temp_df_ppe,
             by = c('year', 'month'))



####### for transforming data back into original units
# ordered jan - december
month_precip_ppe <- sim_precip_df %>%
  mutate(month = as_factor(month)) %>%
  group_by(month) %>%
  summarise(
    precip_mean = mean(log_precip),
    precip_sd = sd(log_precip)
  )

month_temp_ppe <- sim_temp_df %>%
  mutate(month = as_factor(month)) %>%
  group_by(month) %>%
  summarise(
    temp_mean = mean(temp),
    temp_sd = sd(temp)
  )


month_p_means_ppe <- month_precip_ppe %>%
  pull(precip_mean)
month_p_sds_ppe <- month_precip_ppe %>%
  pull(precip_sd)


month_t_means_ppe <- month_temp_ppe %>%
  pull(temp_mean)
month_t_sds_ppe <- month_temp_ppe %>%
  pull(temp_sd)

##### end


source(here('analysis', 'experimental', 'pt_helpers.R'))

snr_ppe <- 0.5

trw_dfs_ppe <- format_trw_p_t_linear(clim_df_ppe,
                                      n_trees = 20,
                                      snr = snr_ppe,
                                      start_year_calib = calib_start_ppe,
                                      w_p_true = w_djf, w_t_true = w_jja,
                                     beta0 = 2, betat = 2, betap = 2, beta_age = -1)

trw_mod_dfs_ppe <- process_trw(trw_dfs_ppe)

halfyear_prior <- c(10,10,rep(1, 6), rep(10, 4))
unif_prior <- rep(1,12)
df_list_ppe <- stan_setup_p_t(clim_df_ppe, trw_mod_dfs_ppe,
                               start_year_calib = calib_start_ppe,
                          w_param_precip = unif_prior,
                          w_param_temp = unif_prior,
                          month_means_t = month_t_means_ppe,
                          month_sds_t = month_t_sds_ppe,
                          month_means_p = month_p_means_ppe,
                          month_sds_p = month_p_sds_ppe)



# df_list_precip <- stan_setup(clim_df_ppe %>% rename(clim = precip), trw_mod_dfs_ppe,
#                           start_year_calib = calib_start_ppe,
#                           w_param = unif_prior, month_means = month_p_means_ppe,
#                           month_sds = month_p_sds_ppe)
#
# # this will take a while to run!
# fit_ppe_precip <- stan(file = here("analysis", "stan_code",
#                                    "proposed_model.stan"),
#                        iter= 1000,
#                        data = df_list_precip,
#                        init_r = 1,
#                        seed = 123#,
#                        #control = list(adapt_delta = 0.99)
# )


# init_r helps avoid non-posdef matrix initialization
# this will take a while to run!
fit_ppe_pt <- stan(file = here("analysis", "experimental",
                               "proposed_model_pt.stan"),
                 iter= 1000,
                 data = df_list_ppe,
                init_r = 1,
                seed = 123#,
                #control = list(adapt_delta = 0.99)
                )



## code if you want to save this output!
# saveRDS(fit_ppe_pt, 'E:/Projects/old-blue-oaks/temp_output/fit_ppe_pt.Rds')


fig_weight_pt(fit_ppe_pt) #+
  # labs(title = str_glue('linear ppe with log(precip + 0.001), snr = {snr_ppe}'))


# reconstruction starts with first year of trw data.
recon_start_ppe <- min(df_list_ppe$start_years)


## repeat, but get the reconstruction in terms of DJF precipitation (original units mm/month)

djf_p_df_ppe <- process_climate(sim_precip_df,
                              clim_name = 'precip',
                              w_agg = w_djf, anom = F)
start_year_obs <- min(djf_p_df_ppe$year)
djf_fig_ppe_df <- get_season_pred(fit_ppe_pt,
                               djf_p_df_ppe, varname = 'precip_djf',
                               recon_start_ppe)

djf_fig_ppe_df <- djf_fig_ppe_df %>%
  mutate(
    start_year_obs = start_year_obs,
    calib_start = calib_start_ppe
  )

# write_csv(djf_fig_ppe_df, 'E:/Projects/old-blue-oaks/temp_output/djf_p_fig_ppe_df.csv')


fig_recon(djf_fig_ppe_df,
          1600, calib_start_ppe,
          display_years = c(1750, 1950)
)[[2]] +
  labs(y = 'DJF mean precip (mm/month)')




#######

jja_t_df_ppe <- process_climate(sim_temp_df,
                                clim_name = 'temp',
                                w_agg = w_jja, anom = F)


jja_fig_ppe_df <- get_season_pred(fit_ppe_pt,
                                  jja_t_df_ppe, varname = 'temp_jja',
                                  recon_start_ppe)

jja_fig_ppe_df <- jja_fig_ppe_df %>%
  mutate(
    start_year_obs = start_year_obs,
    calib_start = calib_start_ppe
  )

# write_csv(jja_fig_ppe_df, 'E:/Projects/old-blue-oaks/temp_output/jja_t_fig_ppe_df.csv')


fig_recon(jja_fig_ppe_df,
          1600, calib_start_ppe,
          display_years = c(1750, 1950)
)[[2]] +
  labs(y = 'JJA mean temp (C)')



