# Author: Khin Nawarat (2026, @IHE Delft)
# Compiling the results and saving as netcdf files

remove(list=ls())
filedir<-"D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # change to your local directory
source(paste0(filedir,"Rcode/allfunctions.r"))

grid_files <- c("common2grids.csv","common3grids.csv")
out_files  <- c("GWL_VS_FC_larger_set.nc","GWL_VS_FC_smaller_set.nc")

whichdecade <- 2100

library(ncdf4)
for (k in 1:2){
  
  gridcoo.commonlocs <- read.csv(paste0(filedir,"Rdatasets/",grid_files[k]))
  
  # location metadata
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
  name  = "return_period",
  units = "years",
  vals  = c(10, 20, 50, 100, 200),
  longname = "Extreme sea level return period"
)

dim_method <- ncdim_def(
  name  = "method",
  units = "1",
  vals  = c(1, 2),
  longname = "Method index"
)

dim_unc <- ncdim_def(
  name  = "uncertainty",
  units = "1",
  vals  = c(1, 2, 3),
  longname = "Uncertainty quantile index"
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

# main variable
var_GWL <- ncvar_def(
  name = "GWL",
  units = "1",
  dim = list(dim_location, dim_RP, dim_method, dim_unc),
  missval = -9999L,
  longname = "Global warming level causing annual extreme sea level event",
  prec = "short",   # values 1–7
  compression = 4
)

# create nc file
ncfile <- nc_create(
  paste0(filedir,"Data/",out_files[k]),
  vars=list(var_lon,var_lat,var_locid,var_GWL),
  force_v4=TRUE
)



# write coordinate data
ncvar_put(ncfile, var_lon, lon)
ncvar_put(ncfile, var_lat, lat)
ncvar_put(ncfile, var_locid, locid)

# method labels
ncatt_put(ncfile, "method", "labels", paste(c("VS","FC"), collapse=","))

# uncertainty labels
ncatt_put(ncfile, "uncertainty", "labels", paste(c("q0.05","q0.5","q0.95"), collapse=","))

# GWL category labels
ncatt_put(ncfile, "GWL", "category_labels"
          , paste(c("1.5","2","2.5","3","4","5",">5"), collapse=","))

# RL category labels
ncatt_put(ncfile, "return_period", "labels"
          , paste(c("10","20","50","100","200"), collapse=","))

# global attributes
ncatt_put(ncfile, 0, "title",
          "Global warming levels causing ESL RLs to annual levels")

ncatt_put(ncfile, 0, "institution",
          "IHE Delft")

if (k == 1) {
  ncatt_put(
    ncfile, 0, "source",
    paste(
      "Voting system (VS) and full convolution (FC) methods.",
      "Present-day ESL estimates based on Vousdoukas et al. (2018) and Kirezci et al. (2020).",
      sep = "\n"
    )
  )
} else {
  ncatt_put(
    ncfile, 0, "source",
    paste(
      "Voting system (VS) and full convolution (FC) methods.",
      "Present-day ESL estimates based on Rasmussen et al. (2018), Vousdoukas et al. (2018), and Kirezci et al. (2020).",
      sep = "\n"
    )
  )
}


ncatt_put(ncfile, 0, "Conventions", "CF-1.8")


# writing data inside the loop

# mapping indices
RL_vals <- c(10,20,50,100,200)
method_index <- c(VS = 1, FC = 2)
unc_index <- c(q0.05 = 1, q0.5 = 2, q0.95 = 3)


for ( RL in c(10,20,50,100,200)){
  load(filedir%&%"Rdatasets/RData_labels_"%&%RL%&%"to1_FC_VS_both_fisherinfo_"%&%whichdecade) 
   RL_index <- match (RL, RL_vals)
   for (method in c("VS", "FC")) {
     for (q in c("q0.05", "q0.5", "q0.95")) {
       prefix <- if (k == 1) {  # larger set
         ifelse(method=="FC","all2","V2")
       } else {                 # smaller set
         ifelse(method=="FC","all3","V3")
       }
       
       obj_name <- paste0(prefix, ".", RL, "to1.", q)
       if (!exists(obj_name)) {
         stop("Missing object: ", obj_name)
       }
       
       gwl_vec <- get(obj_name)
       print(obj_name)
       print(gwl_vec)
       ncvar_put(
         ncfile,
         var_GWL,
         vals  = gwl_vec,
         start = c(
           1,
           RL_index,
           method_index[method],
           unc_index[q]
         ),
         count = c(nloc, 1, 1, 1)
       )
     }
   }
   
}

nc_close(ncfile)
}



