# https://www.ncei.noaa.gov/pub/data/paleo/treering/measurements/northamerica/usa/ak132x.rwl


source(here::here("analysis", "r_code", "helpers", "analysis_helpers.R"))
source(here::here("analysis", "r_code", "helpers", "plot_helpers.R"))

library(dplR)

nc <- parallelly::availableCores()
options(mc.cores = nc)

calib_start <- 1900
# recon_start <- 1900
unif_prior <- rep(1, 12)



firth_trw <- dplR::read.rwl(here('data', 'ak132x.rwl'))

firth_trw_df <- firth_trw %>%
  rownames_to_column('year') %>%
  mutate(year = as.numeric(year)) %>%
  # filter(year >= calib_start) %>% # new!
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
  # filter(year > 1892) %>%
  mutate(month = as_factor(month))


## eda -------------

chron_firth <- firth_trw %>%
  detrend(method = 'ModHugershoff') %>%
  chron() %>%
  rownames_to_column('year') %>%
  mutate(year = as.numeric(year)) %>%
  select(year, chron = std)



w_djf <- c(1, 1, rep(0, 9), 1) / 12
w_jja <- c(rep(0, 5), 1, 1, 1, rep(0, 4)) / 12
w_jas <- c(rep(0, 6), 1, 1, 1, rep(0, 3)) / 12

temp_chron <- process_climate(
  firth_temp_df,
  clim_name = 'tmean',
  anom = T, w_agg = w_jja
) %>%
  right_join(chron_firth, by = 'year')


cor(temp_chron$clim, temp_chron$chron, use = 'complete.obs')


temp_chron %>%
  mutate(clim = (clim - mean(clim, na.rm = T)) / sd (clim, na.rm = T),
         chron = (chron - mean(chron, na.rm = T)) / sd(chron, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(year, clim), color = 'blue') +
  geom_line(aes(year, chron), color = 'red') +
  lims(x = c(1800, 2000))



###########################

# temp

##########################

temp_df_wateryear <- process_climate(
  firth_temp_df,
  clim_name = 'tmean',
  anom = T
) #%>%
  # filter(year >= calib_start) # new!


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
input_list_firth_chron <- stan_setup(temp_df_wateryear, firth_trw_dfs_mod,
                                   start_year_calib = calib_start,
                                   w_param = unif_prior,
                                   month_means = month_means_temp,
                                   month_sds = month_sds_temp)

input_list_firth_chron$N_trees <- NULL
input_list_firth_chron$age <- NULL
input_list_firth_chron$s <- NULL
input_list_firth_chron$start_years <- NULL
input_list_firth_chron$start_year <- min(chron_firth$year)
input_list_firth_chron$trw <- chron_firth$chron




fit_firth_chron <- stan(file = here("analysis",
                                        "experimental",
                                        "chron_model.stan"),
                            iter= 1000,
                            data = input_list_firth_chron,
                            init_r = 1,
                            seed = 123
)

# saveRDS(fit_firth_chron, "E:/Projects/old-blue-oaks/temp_output/fit_firth_chron.Rds")

## Evaluation --------------------

(w_learned_temp <- fig_weight(fit_firth_chron))

recon_start <- input_list_firth_chron$start_year
yearly_fig_df <- get_yearly_pred(fit_firth_chron, temp_df_wateryear, recon_start)

# start_year_obs <- min(temp_df_wateryear$year)
start_year_obs <- 1900
fig_recon(yearly_fig_df, start_year_obs, calib_start,
          pred_name = "clim_mean_pred",
          truth_name = "mean_w_trueclim",
          display_years = c(1700, 2007))[[2]] +
  labs(title = 'eg',
       y = 'log precip anom (sd)'
  )



