#####################################################################

# Common libraries and helper functions for cleaning data

#####################################################################

library(tidyverse)
library(here)
library(magrittr)
library(ncdf4)
library(lubridate)
library(tidync)
library(boxr)
library(patchwork)

theme_set(theme_bw(base_size = 20))



#' quick test function to check that the years match up correctly in the formatted df
#' if everything worked correctly, then all multiples of 10 should look correct (test_tens)
#' if everything worked correctly, then each tree should have sequential years from its min to max, no gaps
#' @param df is a long dataframe of tree ring widths, with column names "year", "year2", and "tree_id" (see long_df in the function trw_to_long)
#'            year is the original column, year2 is the corrected column
#' @return TRUE if both conditions are met, and the dataframe looks good.
check_year <- function(df){
  # test multiples of 10
  test_tens <- df %>%
    filter(year %% 10 == 0) %>%
    group_by(year) %>%
    slice(1) %$%
    all(year == year2)

  # test that no years are missing
  test_seq <- df %>%
    group_by(tree_id) %>%
    summarize(test = all.equal(year2, min(year2, na.rm = T):max(year2, na.rm = T))) %>%
    pull(test) %>%
    all()

  # return TRUE if both are correct
  all(test_tens, test_seq)
}

#' function to format the trw data into long format
#' @param raw_tree_df is the output of the files read in directly from the paleo data search
trw_to_long <- function(raw_tree_df){
  long_df <- raw_tree_df %>%
    pivot_longer(-c(X1, X2)) %>%
    select(tree_id = X1, year = X2, trw = value) %>%
    filter(!is.na(trw)) %>%
    group_by(tree_id) %>%
    mutate(years_since_start = row_number() - 1,
           start_year = year[1],
           year2 = start_year + years_since_start
    ) %>%
    ungroup()

  if(check_year(long_df)){
    formatted_df <- long_df %>%
      transmute(tree_id,
                year = year2,
                age = years_since_start + 1,
                trw)

    return(formatted_df)
  } else{
    stop("something went wrong with the formatting! check it out")
  }
}

#' check for duplicates in tree id-decades:
#' if there are duplicates, then return a list of the problimatic trees
check_for_dup <- function(raw_df){
  num_treedecade <- raw_df %>%
    group_by(X1,X2) %>%
    summarize(n = n(), .groups = "drop")

  is_dup <- all(num_treedecade$n == 1)

  if(is_dup){
    message("all good, no duplicates")
  } else{
    message("watch out! there are duplicates")
    num_treedecade %>%
      filter(n > 1)
  }
}



#' get time
#' @param nc is a netcdf object
getNcTime <- function(nc) {
  require(units)
  require(ncdf4)
  options(warn=1) #show warnings by default
  if (is.character(nc)) nc <- nc_open(nc)
  ncdims <- names(nc$dim) #get netcdf dimensions
  timevar <- ncdims[which(ncdims %in% c("time", "Time", "datetime", "Datetime", "date", "Date"))] #find (first) time variable
  if (length(timevar) > 1) {
    warning(paste("Found more than one time var. Using the first:", timevar[1]))
    timevar <- timevar[1]
  }
  if (length(timevar)!=1) stop("ERROR! Could not identify the correct time variable")
  times <- ncvar_get(nc, timevar) #get time data
  timeatt <- ncatt_get(nc, timevar) #get attributes
  timeunit <- timeatt$units

  timeunit

}

# ncatt_get()
# start is the first day of the month
avg_monthly_daylength <- function(lat, start){
  # list all days in the month
  days_in_month <- start + days(0:(days_in_month(start)-1))

  # get mean day length for the month
  mean(geosphere::daylength(lat = lat, doy = days_in_month))

}

#' quick little function to get lubridate formatted june/december of the year
#' useful for getting solstice month
#' note: need all of the year_vec to be the same hemisphere
get_solstice_formatted <- function(year_vec, hemi = "north"){
  if(hemi == 'north'){
    solstice_month <- "06"
  }else{
    solstice_month <- "12"
  }

  # if we're before the year 1000, need to add a 0 starting digit
  year_vec_char <- ifelse(nchar(year_vec) == 3,
                          as.character(year_vec) %>%
                            paste0("0", .),
                          as.character(year_vec)
                          )

  my(str_glue("{solstice_month}-{year_vec_char}"))
}


lon_to_180 <- function(lon360){
  ifelse(lon360 > 180, lon360 - 360, lon360)
}

lon_to_360 <- function(lon180){
  ifelse(lon180 < 0, 360 + lon180, lon180)
}



