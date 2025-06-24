# https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ak132x.rwl


source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

library(dplR)

nc <- parallelly::availableCores()
options(mc.cores = nc)

calib_start <- 1950
# recon_start <- 1950
unif_prior <- rep(1, 12)



firth_trw <- dplR::read.rwl(here('data', 'ak132x.rwl'))

firth_trw_df <- firth_trw %>%
  rownames_to_column('year') %>%
  mutate(year = as.numeric(year)) %>%
  pivot_longer(-year, names_to = 'tree_id', values_to = 'trw') %>%
  drop_na() %>%
  group_by(tree_id) %>%
  arrange(year, .by_group = T) %>%
  mutate(age = 1:n())

firth_trw_df_list <- group_split(firth_trw_df)
firth_trw_dfs_mod <- process_trw(firth_trw_df_list)



##################################

# get climate data
# note: the paper does some clever work to get these quantities
# i'm just going to simplify and use prism

# eg. white moutnian summit 37.634 ◦N, -118.256 ◦E

trw_lastyear <- max(firth_trw_df$year)

# 68.65 N, Longitude: -141.63, E
firth_lonlat <- c(-141.63, 68.65)

# there is some missing data in 1892 (and in earlier years)
firth_temp_df <- read_csv(here('data', 'berk_firth.csv')) %>%
  filter(year > 1892) %>%
  mutate(month = as_factor(month))

## eda -------------

ars_chron <- chron.ars(firth_trw) %>%
  rownames_to_column('year') %>%
  mutate(year = as.numeric(year)) %>%
  select(year, ars)

w_djf <- c(1, 1, rep(0, 9), 1) / 3
w_jja <- c(rep(0, 5), 1, 1, 1, rep(0, 4)) / 3
w_jas <- c(rep(0, 6), 1, 1, 1, rep(0, 3)) / 3

temp_chron <- process_climate(
  firth_temp_df,
  clim_name = 'tmean',
  anom = T, w_agg = w_jja
) %>%
  right_join(ars_chron, by = 'year')

cor(temp_chron$clim, temp_chron$ars, use = 'complete.obs')


temp_chron %>%
  mutate(clim = (clim - mean(clim, na.rm = T)) / sd (clim, na.rm = T),
         ars = (ars - mean(ars, na.rm = T)) / sd(ars, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(year, clim), color = 'blue') +
  geom_line(aes(year, ars), color = 'red') +
  lims(x = c(1800, 2000))


#####

firth_temp_df %>%
  left_join(ars_chron, by = 'year') %>%
  group_by(month) %>%
  summarise(r = cor(tmean, ars, use = 'complete.obs'))

###########################

# temp

##########################

temp_df <- process_climate(
  firth_temp_df,
  clim_name = 'tmean',
  anom = T,
  water_year = F
)


# ordered jan - december
month_temp_firth <- firth_temp_df %>%
  mutate(month = as_factor(month)) %>%
  group_by(month) %>%
  summarise(
    clim_mean = mean(tmean),
    clim_sd = sd(tmean)
  )

month_means_temp <- month_temp_firth$clim_mean
month_sds_temp <- month_temp_firth$clim_sd



## set up data ---------------


# halfyear_prior <- c(10,10,rep(1, 6), rep(10, 4))
input_list_firth <- stan_setup(temp_df, firth_trw_dfs_mod,
                                   start_year_calib = calib_start,
                                   w_param = unif_prior,
                                   month_means = month_means_temp,
                                   month_sds = month_sds_temp)

input_list_firth$w_jja <- w_jja

fit_firth <- stan(file = here("analysis",
                                  "stan_code",
                                  "proposed_model_firth.stan"),
                      iter= 1000,
                      data = input_list_firth,
                      init_r = 1,
                      seed = 123
)


# saveRDS(fit_firth, "E:/Projects/old-blue-oaks/temp_output/fit_firth.Rds")



## Evaluation --------------------

(w_learned_temp <- fig_weight(fit_firth))



recon_start <- min(input_list_firth$start_years)
yearly_fig_df <- get_yearly_pred(fit_firth, temp_df, recon_start)
yearly_fig_df <- yearly_fig_df %>%
  mutate(
    firth_start_year_obs = start_year_obs,
    calib_start_firth = calib_start
  )
# write_csv(yearly_fig_df, "E:/Projects/old-blue-oaks/temp_output/yearly_fig_df_firth.csv")



# start_year_obs <- min(temp_df$year)
start_year_obs <- 1900
fig_recon(yearly_fig_df, start_year_obs, calib_start,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim",
          display_years = c(1700, 2007))[[2]] +
  labs(title = 'eg',
       y = 'log precip anom (sd)'
  )



# just for evaluation
jja_df_firth <- process_climate(firth_temp_df,
                                clim_name = 'tmean',
                                w_agg = w_jja, anom = F)

jja_df_firth <- jja_df_firth %>%
  mutate(
    firth_start_year_obs = start_year_obs,
    calib_start_firth = calib_start
  )

jja_fig_firth_df <- get_season_pred(fit_firth, jja_df_firth,
                                    varname = 'clim_jja',
                                    recon_start)

# write_csv(jja_fig_firth_df, 'E:/Projects/old-blue-oaks/temp_output/jja_fig_firth_df.csv')


