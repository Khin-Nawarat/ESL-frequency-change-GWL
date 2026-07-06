# Author: Khin Nawarat (2026, @IHE Delft)
# saving timing as netcdf files

remove(list=ls())
library(ncdf4)

filedir <- "D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # Change to your local directory


# datasets for two location sets
grid_files <- c("common2grids.csv", "common3grids.csv")
rdata_files <- c("RData_first_decade_FC_allRLs_q_wls2",
                 "RData_first_decade_FC_allRLs_q_wls3")
out_files <- c("Timing_FC_larger_set.nc",
               "Timing_FC_smaller_set.nc")

for (k in 1:2) {
  
gridcoo.commonlocs <- read.csv(paste0(filedir,"Rdatasets/",grid_files[k]))
load(paste0(filedir,"Rdatasets/",rdata_files[k]))  # loads first_decade

lon  <- gridcoo.commonlocs$longitude
lat  <- gridcoo.commonlocs$latitude
locid <- gridcoo.commonlocs$X
nloc <- length(lon)

# dimensions
dim_location <- ncdim_def(
  name  = "location",
  units = "1",
  vals  = seq_len(nloc),
  longname = "Index of coastal locations"
)
dim_RP <- ncdim_def(
  name = "return_period",
  units = "years",
  vals = c(10, 20, 50, 100, 200),
  longname = "Extreme sea level return period"
)

dim_unc <- ncdim_def(
  name = "uncertainty",
  units = "1",
  vals = 1:3,
  longname = "Uncertainty quantile index"
)

dim_GWL <- ncdim_def(
  name = "GWL",
  units = "1",
  vals = 1:6,   # six GWLs
  longname = "Global warming level index"
)


# coordinate variables
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
# main variable: first decade
var_firstdecade <- ncvar_def(
  name = "first_decade",
  units="year",
  dim = list(dim_location, dim_unc, dim_GWL, dim_RP),
  missval = -9999L,
  longname = "First decade when a given RP becomes annual (FC method)",
  prec = "short",
  compression = 4
)

# create NetCDF
ncfile <- nc_create(paste0(filedir,"Data_test/",out_files[k]),
                    vars = list(var_lon, var_lat, var_locid, var_firstdecade),
                    force_v4 = TRUE)

# write coordinates
ncvar_put(ncfile, var_lon, lon)
ncvar_put(ncfile, var_lat, lat)
ncvar_put(ncfile, var_locid, locid)

# labels
ncatt_put(ncfile, "uncertainty", "labels", paste(c("q0.05","q0.5","q0.95"), collapse = ","))
ncatt_put(ncfile, "GWL", "labels", paste(c("1.5C","2.0C","2.5C","3.0C","4.0C","5.0C"), collapse = ","))
ncatt_put(ncfile, "return_period", "labels", paste(c("10","20","50","100","200"), collapse = ","))

# global attributes
ncatt_put(ncfile, 0, "title", "First decade when ESL return levels become annual")
ncatt_put(ncfile, 0, "institution", "IHE Delft")
ncatt_put(ncfile, 0, "source", "with Full convolution (FC) method")

if (k == 1) {
  ncatt_put(
    ncfile, 0, "source",
    paste(
      "with full convolution (FC) methods.",
      "Present-day ESL estimates based on Vousdoukas et al. (2018) and Kirezci et al. (2020).",
      sep = "\n"
    )
  )
} else {
  ncatt_put(
    ncfile, 0, "source",
    paste(
      "with full convolution (FC) methods.",
      "Present-day ESL estimates based on Rasmussen et al. (2018), Vousdoukas et al. (2018), and Kirezci et al. (2020).",
      sep = "\n"
    )
  )
}

ncatt_put(ncfile, 0, "Conventions", "CF-1.8")

# write data
ncvar_put(ncfile, var_firstdecade, first_decade)

nc_close(ncfile)

}
