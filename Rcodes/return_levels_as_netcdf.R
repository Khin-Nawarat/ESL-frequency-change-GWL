# Author: Khin Nawarat (2026, @IHE Delft)
# Saving return level as netcdf files


remove(list=ls())
library(ncdf4)

filedir <- "D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # Change to your local directory
source(paste0(filedir, "Rcode/allfunctions.r"))

grid_files <- c("common2grids.csv", "common3grids.csv")
outdir <- paste0(filedir,"Data/")
out_files  <- c("Return_levels_larger_set.nc", "Return_levels_smaller_set.nc")
whichdecade <- 2100 # change here to get return levels for other decades between 2020-2100
RLs <- c(10, 20, 50, 100, 200)
uncertainty_levels <- c("q0.05", "q0.5", "q0.95")
climate_states <- c("current", rasmussenScenarios)

for (k in 1:2) {
  gridcoo.commonlocs <- read.csv(paste0(filedir, "Rdatasets/", grid_files[k]))
  lon <- gridcoo.commonlocs$longitude
  lat <- gridcoo.commonlocs$latitude
  locid <- gridcoo.commonlocs$X
  nloc <- length(lon)

  dim_location <- ncdim_def(
    name = "location",
    units = "1",
    vals = seq_len(nloc),
    longname = "Index of coastal locations"
  )

  dim_RL <- ncdim_def(
    name = "return_period",
    units = "years",
    vals = RLs,
    longname = "Extreme sea level return period"
  )

  dim_uncertainty <- ncdim_def(
    name = "uncertainty",
    units = "1",
    vals = seq_along(uncertainty_levels),
    longname = "Uncertainty quantile index"
  )

  dim_climate <- ncdim_def(
    name = "GWL",
    units = "1",
    vals = seq_along(climate_states),
    longname = "Global warming level index"
  )

  var_lon <- ncvar_def(
    name = "longitude",
    units = "degrees_east",
    dim = list(dim_location),
    missval = NA_real_,
    longname = "Longitude of coastal location",
    prec = "double"
  )

  var_lat <- ncvar_def(
    name = "latitude",
    units = "degrees_north",
    dim = list(dim_location),
    missval = NA_real_,
    longname = "Latitude of coastal location",
    prec = "double"
  )

  var_locid <- ncvar_def(
    name = "location_id",
    units = "1",
    dim = list(dim_location),
    missval = NA_real_,
    longname = "Original grid location identifier",
    prec = "double"
  )

  var_rl <- ncvar_def(
    name = "return_level_value",
    units = "cm",
    dim = list(dim_location, dim_RL, dim_uncertainty, dim_climate),
    missval = -9999.0,
    longname = "Extreme sea level for a given return period  under present-day and future global warming levels",
    prec = "float",
    compression = 4
  )
  # create NetCDF

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  ncfile <- nc_create(
    paste0(outdir, out_files[k]),
    vars = list(var_lon, var_lat, var_locid, var_rl),
    force_v4 = TRUE
  )
  # write coordinates
  ncvar_put(ncfile, var_lon, lon)
  ncvar_put(ncfile, var_lat, lat)
  ncvar_put(ncfile, var_locid, locid)

  # labels
  ncatt_put(ncfile, "uncertainty", "labels", paste(uncertainty_levels, collapse = ","))
  ncatt_put(ncfile, "GWL", "labels", paste(climate_states, collapse = ","))
  ncatt_put(ncfile, "return_period", "labels", paste(as.character(RLs), collapse = ","))

  # global attributes
  ncatt_put(ncfile, 0, "title", "Extreme sea levels for different return-period and global warming levels")
  ncatt_put(ncfile, 0, "institution", "IHE Delft")
  ncatt_put(ncfile, 0, "source", "Full convolution (FC) method")
  ncatt_put(ncfile, 0, "Conventions", "CF-1.8")

  rl_data <- array(
    -9999.0,
    dim = c(nloc, length(RLs), length(uncertainty_levels), length(climate_states))
  )

  unc_index <- setNames(seq_along(uncertainty_levels), uncertainty_levels)
  climate_index <- setNames(seq_along(climate_states), climate_states)
  
  # write data
  for (RL in RLs) {
    # using the full-convolution projections generated in RData_fullconvolutions_fisherinfo files.
    load(paste0(filedir, "Rdatasets/RData_fullconvolutions_fisherinfo_", RL, "yr_", whichdecade))

    if (k == 1) {
      int_proj <- all2f.int.probproj
    } else {
      int_proj <- all3f.int.probproj
    }

    for (unc in uncertainty_levels) {
      row_name <- unc
      if (!(row_name %in% dimnames(int_proj)[[1]])) {
        stop("Missing quantile row: ", row_name)
      }
      for (cl in climate_states) {
        if (!(cl %in% dimnames(int_proj)[[2]])) {
          stop("Missing climate state: ", cl)
        }
        rl_data[, match(RL, RLs), unc_index[unc], climate_index[cl]] <-
          round(int_proj[row_name, cl, ], 1)
      }
    }
  }

  ncvar_put(ncfile, var_rl, rl_data)
  nc_close(ncfile)
}
