##########################

# Outputs:

# - Blue Oaks TRWs from Los Lobos and Figueroa Mountain  (blueoaks_long.csv)
# - CESM last millenium run time series of precip and temperature from nearest gridcell (sim_cesm1_mil_data.csv)
# - PRISM precipitation time series from nearest gridcell (prism_ts.csv)

# see data availability for data sources!

#########################


source(here::here('data_preprocessing', 'data_cleaning_helpers.R'))
library(SPEI) # get PET for the optional PDSI calculation
library(raster)
library(prism)# install.packages('prism')


####################################

# Read in observational tree ring data
# (los lobos)

####################################

#https://www.ncei.noaa.gov/access/paleo-search/study/search.json?dataPublisher=NOAA&dataTypeId=18&locations=Continent%3ENorth%20America%3EUnited%20States%20Of%20America%3ECalifornia&species=QUDG&headersOnly=true

trw_base_path <- here("data", "loslobos_rockspringsranch_figueroamountain_20220104",
                      "data","pub","data","paleo","treering",
                      "measurements","northamerica","usa")


# lonlat are taken from the file path, ca646-noaa.crn (for)
loslobos_lonlat <-  c(-119.24, 34.92)
figueroa_lonlat <- c(-120, 34.74)

raw_loslobos <- read_fwf(file.path(trw_base_path, "ca660.rwl"),
                         fwf_empty(file.path(trw_base_path,"ca656.rwl"),
                                   skip = 3, n = 1e3), skip = 3)
raw_figueroa <- read_fwf(file.path(trw_base_path, "ca656.rwl"),
                        fwf_empty(file.path(trw_base_path,"ca656.rwl"),
                                  skip = 3, n = 1e3), skip = 3)


#### Tree FIG84A looks problematic... it looks like there might be two different measurements for this tree? Why isn't the second one FIG84B???

# let's fix it up for now

figueroa_no_dup <- raw_figueroa %>%
  group_by(X1, X2) %>%
  mutate(is_dup = row_number(X1) == 2) %>%
  ungroup() %>%
  mutate(X1 = ifelse(is_dup, str_glue("{X1}_dup"), X1)) %>%
  select(-is_dup)

# manually fix FIG103, where the first and last row are missed for dup


# according to the documentation: "units are 0.01 mm if end-of-series marker is 999 and 0.001 mm if end-of-series marker is -9999"
# note: assuming that 9999 is same as -9999 for now
trw_df_long <- list(rockspring = raw_rockspring,
                    loslobos = raw_loslobos,
                    figueroa = figueroa_no_dup) %>%
  map_dfr(~trw_to_long(.x), .id = "site") %>%
  filter(abs(trw) != 9999)


# write_csv(trw_df_long, file = here("processed_data",
#                                    "blueoaks_long.csv"))


####################################

# Read in CESM1 last millenium data

####################################

#' get temperature, precip from a last millennium run, and estimate pdsi
#' @param mylat
#' @param mylon -- assumed to be in -180 to 180 format
get_sim_cesm1_mil <- function(mylat, mylon){

  ## Change these paths to wherever the data is downloaded locally!
  # The first ensemble member of control run looks something like
  # b.e11.BLMTRC5CN.f19_g16.001.clm2.h0.SOILWATER_10CM.085001-184912.nc

  temp_mil_paths <- list.files("C:/Users/natha/Box/climate/data/last_mil/TREFHT/", pattern = "b.e11.BLMTRC5CN.f19_g16.001.cam.h0.TREFHT.*", full.names = T)
  precip_mil_paths <- list.files("C:/Users/natha/Box/climate/data/last_mil/PRECT/", pattern = "b.e11.BLMTRC5CN.f19_g16.001.cam.h0.PRECT.*", full.names = T)

  # read in .nc files
  all_temps_mil <- map(temp_mil_paths, ~ tidync(.x))
  all_precips_mil <- map(precip_mil_paths, ~ tidync(.x))


  # for the first nc file (0850-1850) time is in days from 850 (getNcTime(all_temp_mils[1]))
  # for the second nc file (1850-), time is in days from 1850

  # but... the dates get messed up in some of the februaries, so instead.
  # instead, I use the rownumber to increment months
  sim_mil_dflist <- list(all_temps_mil, all_precips_mil) %>%
    map_depth(.depth = 2, ~hyper_filter(.x, lat = index == which.min(abs(lat - mylat)),
                                        lon = index == which.min(abs(lon - lon_to_360(mylon)))) %>%
                hyper_tibble() %>%
                mutate(lat = round(lat, 5))
    ) %>%
    map(function(y) list_modify(y,
                                y[[1]] %>%
                                  rowid_to_column("index") %>%
                                  mutate(time_new = my("01-0850") + months(index-1)),
                                y[[2]] %>%
                                  rowid_to_column("index") %>%
                                  mutate(time_new = my("01-1850") + months(index-1))
    )
    ) %>%
    map_depth(.depth = 1, ~bind_rows(.x) %>% select(-index))

  # 8.64e7 mm/day = 1 m/s
  # x m/s =  x * 8.64e7 mm/day * 30 day/ month
  sim_mil_df <- sim_mil_dflist %>%
    reduce(inner_join, by = c("lon", "lat", "time", "time_new")) %>%
    relocate(TREFHT, .after = PRECT) %>%
    mutate(temp_c = TREFHT - 273.15,
           precip_mmmonth = PRECT * 8.64e7 * 30,
           temp_scaled = (TREFHT - mean(TREFHT)) / sd(TREFHT),
           precip_scaled = (PRECT - mean(PRECT)) / sd(PRECT),
           avg_day = map_dbl(time_new, ~avg_monthly_daylength(68, .x)),
           solstice = get_solstice_formatted(year(time_new), hemi = 'north'),
           year = year(time_new),
           month = month(time_new)
    ) %>%
    mutate(
      water_year = ifelse(month >= 10, year + 1, year),
      avg_solstice = map_dbl(solstice, ~avg_monthly_daylength(68, .x)),
      PET = as.numeric(SPEI::thornthwaite(temp_c, 68))
    )

  # store month info as both abbreviation and number
  month_to_number <- enframe(toupper(month.abb), name = "month", value = "name")

  # optional pdsi calculation, if we want that for the analysis
  PDSI_mil <- pdsi::pdsi(awc = 100, lat = mylat, climate =
                           sim_mil_df %>%
                           transmute(year,
                                     month,
                                     temp = temp_c, prec = precip_mmmonth) %>%
                           as.data.frame(),
                         start = 1000, end = 2004, mode = "scpdsi") %>%
    pivot_longer(-YEAR, values_to = "scPDSI") %>%
    left_join(month_to_number, by = "name")


  # remove first year (no pdsi estimate is produced that year)
  sim_cesm1_mil_data <- PDSI_mil %>%
    left_join(sim_mil_df, by = c("YEAR" = "year", "month")) %>%
    rename(year = YEAR) %>%
    filter(water_year < 2005) # last year doesn't have full water year data

  sim_cesm1_mil_data
}



# the los lobos site
sim_cesm1_mil_data <- get_sim_cesm1_mil(mylon = -119.24, mylat = 34.92)

# write_csv(sim_cesm1_mil_data, file = here("processed_data",
#           "sim_cesm1_mil_data.csv"))

###################

# prism  data

####################


## set this directory to whether you would like the prism files to be downloaded!
# prism_set_dl_dir("C:/Users/natha/Box/climate/data/prism")
prism_set_dl_dir(here::here('data', 'prism'))
get_prism_monthlys(type = "ppt", year = 1895:2014, mon = 1:12, keepZip = FALSE)
get_prism_monthlys(type = "tmean", year = 1895:2014, mon = 1:12, keepZip = FALSE)

#' get prism data at a given lon-lat.
#' helper function taken and modified from the prism documentation
#' @param pd is the prism data, e.g. the output of prism_archive_subset()
#' @param location is a vector of (longitude, latitude) with longitude in -180 to 180 format
new_prism_slice <- function (pd, location)
{
  if (!is.null(dim(pd))) {
    stop("You must enter a vector of prism data, not a data frame.\n",
         "Try  prism_archive_subset().")
  }
  if (length(location) != 2 || !is.numeric(location)) {
    stop("`location` should be a numeric vector with length=2.")
  }
  ptype <- unique(pd_get_type(pd))
  if (length(ptype) != 1) {
    stop("`pd` includes multiple variables (", ptype, ").\n",
         "Please ensure that only one variable type is provided to `pd_slice()`.")
  }
  meta_d <- pd_get_date(pd)
  meta_names <- pd_get_name(pd)[1]
  param_name <- strsplit(meta_names, "-")[[1]][3]
  pstack <- pd_stack(pd)
  data <- unlist(raster::extract(pstack, matrix(location, nrow = 1),
                                 buffer = 10))
  data <- as_tibble(data)
  data$date <- as.Date(meta_d)
  data <- data[order(data$date), ]

  data
}

# get monthly precip data for each month
prism_precip <- prism_archive_subset(
  "ppt", "monthly", mon = 1:12
)

# extract data at particular gridcell
lobos_lonlat <- c(-119.24, 34.92)
prism_lobos <- new_prism_slice(prism_precip, lobos_lonlat) %>%
  transmute(
    year = year(date),
    month = month(date, label = T),
    prcp = value
  )

# write_csv(file = here('data', 'prism_casestudy.csv'))


### ---

# get monthly precip data for each month
prism_tmean <- prism_archive_subset(
  "tmean", "monthly", mon = 1:12
)


whitemt_lonlat <- c(-118.256, 37.634)
prism_precip_whitemt <- new_prism_slice(prism_precip, whitemt_lonlat) %>%
  transmute(
    year = year(date),
    month = month(date, label = T),
    prcp = value
  )
prism_temp_whitemt <- new_prism_slice(prism_tmean, whitemt_lonlat) %>%
  transmute(
    year = year(date),
    month = month(date, label = T),
    tmean = value
  )

prism_whitemt <- inner_join(prism_precip_whitemt,
           prism_temp_whitemt,
           by = c('year', 'month'))

# write_csv(prism_whitemt, file = here('data', 'prism_whitemtn.csv'))





#### firth data

# firth_lonlat <- c(-141.63, 68.65)

# # # using hadcrut temp
# get_HadCRUT_monthly_temp <- function(my_lon, my_lat, hadcrut_tidync){
#
#   # for some reason hyper tibble wasn't returning lon/lat/time
#   temp_raw <- hadcrut_tidync %>%
#     hyper_filter(latitude = index == which.min(abs(latitude - my_lat)),
#                  longitude = index == which.min(abs(longitude - my_lon))) %>%
#     hyper_tbl_cube()
#   #hyper_tibble()
#
#   # temp
#   temp_tibble <- enframe(temp_raw$mets$tas_mean, name = 'time1', value = "tas_mean") %>%
#     mutate(time = temp_raw$dims$time,
#            lon = temp_raw$dims$longitude,
#            lat = temp_raw$dims$latitude) %>%
#     mutate(time = myd("01-1850-01") + days(as.integer(time)))
#
#   message(str_glue({"lonlat is ({temp_tibble$lon[1]}, {temp_tibble$lat[1]})"}))
#
#
#   # go back and check to see if rounding the time to an integer is ok.
#   temp_tibble %>%
#     mutate(
#       year = year(time),
#       month = month(time, label = T),
#       trw_lon = my_lon,
#       trw_lat = my_lat
#     ) %>%
#     select(
#       year,
#       month,
#       tas_mean,
#       temp_lon = lon,
#       temp_lat = lat,
#       trw_lon,
#       trw_lat
#     ) %>%
#     # pivot_wider(names_from = 'month', values_from = 'tas_mean', id_cols = 'year') %>%
#     drop_na()
# }
#
#
# temp_path <- 'E:/Projects/old-blue-oaks/data/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc'
# hadcrut <- tidync(temp_path)
#
# hadcrut_firth <- get_HadCRUT_monthly_temp(firth_lonlat[1], firth_lonlat[2], hadcrut)
# do transformations to climate first
# firth_temp_df <- hadcrut_firth %>%
#   dplyr::select(
#     year,
#     month,
#     tmean = tas_mean
#   )
#

# gistemp -----------------
# gistemp_path <- 'E:/Projects/old-blue-oaks/data/gistemp250_GHCNv4.NC'
# gis <- tidync(gistemp_path)
#
# get_gis_monthly_temp_dec1994 <- function(my_lon, my_lat, gis_tidync){
#
#   # for some reason hyper tibble wasn't returning lon/lat/time
#   temp_raw <- gis_tidync %>%
#     hyper_filter(lat = index == which.min(abs(lat - my_lat)),
#                  lon = index == which.min(abs(lon - my_lon))) %>%
#     hyper_tbl_cube()
#   #hyper_tibble()
#
#   # temp
#   temp_tibble <- enframe(temp_raw$mets$tempanomaly, name = 'time1', value = "temp_anom") %>%
#     mutate(time = temp_raw$dims$time,
#            lon = temp_raw$dims$lon,
#            lat = round(temp_raw$dims$lat, 5)) %>%
#     mutate(time = myd("01-1800-01") + days(as.integer(time)))
#
#   message(str_glue({"lonlat is ({temp_tibble$lon[1]}, {temp_tibble$lat[1]})"}))
#
#
#   # go back and check to see if rounding the time to an integer is ok.
#   temp_tibble %>%
#     mutate(
#       year = year(time),
#       month = month(time, label = T),
#       trw_lon = my_lon,
#       trw_lat = my_lat
#     ) %>%
#     dplyr::select(
#       year,
#       month,
#       temp = temp_anom,
#       temp_lon = lon,
#       temp_lat = lat,
#       trw_lon,
#       trw_lat
#     ) %>%
#     mutate(temp = ifelse(year == 1994 & is.na(temp), 0, temp)) %>%
#     # pivot_wider(names_from = 'month', values_from = 'tas_mean', id_cols = 'year') %>%
#     drop_na(temp)
# }
# gistemp_firth <- get_gis_monthly_temp_dec1994(firth_lonlat[1], firth_lonlat[2], gis)
#
#
# # do transformations to climate first
# firth_temp_df <- gistemp_firth %>%
#   dplyr::select(
#     year,
#     month,
#     tmean = temp
#   )

# using berkely earth temp -----------------
temp_path <- "E:/Projects/old-blue-oaks/data/berkeleyearth_TAVG_LatLong1.nc"
berk <- tidync(temp_path)

ncmeta::nc_atts(temp_path, "time") %>%
  tidyr::unnest(cols = c(value))

#' from https://berkeleyearth.org/data/, global monthly land
get_berkeleyearth_monthly_temp <- function(my_lon, my_lat, berk_tidync){

  # for some reason hyper tibble wasn't returning lon/lat/time
  temp_raw <- berk_tidync %>%
    hyper_filter(latitude = index == which.min(abs(latitude - my_lat)),
                 longitude = index == which.min(abs(longitude - my_lon))) %>%
    hyper_tbl_cube()
  #hyper_tibble()

  # temp
  temp_tibble <- enframe(temp_raw$mets$temperature, name = 'time1', value = "temperature") %>%
    mutate(time = as.numeric(temp_raw$dims$time),
           lon = temp_raw$dims$longitude,
           lat = temp_raw$dims$latitude) %>%
    mutate(time = date_decimal(time))

  message(str_glue({"lonlat is ({temp_tibble$lon[1]}, {temp_tibble$lat[1]})"}))


  # go back and check to see if rounding the time to an integer is ok.
  temp_tibble %>%
    mutate(
      year = year(time),
      month = month(time, label = T),
      trw_lon = my_lon,
      trw_lat = my_lat
    ) %>%
    dplyr::select(
      year,
      month,
      temp = temperature,
      temp_lon = lon,
      temp_lat = lat,
      trw_lon,
      trw_lat
    ) %>%
    # pivot_wider(names_from = 'month', values_from = 'tas_mean', id_cols = 'year') %>%
    drop_na()
}
berk_firth <- get_berkeleyearth_monthly_temp(firth_lonlat[1],
                                    firth_lonlat[2],
                                    berk)


firth_temp_df <- berk_firth %>%
  dplyr::select(
    year,
    month,
    tmean = temp
  )
# write_csv(firth_temp_df, file = here('data', 'berk_firth.csv'))

