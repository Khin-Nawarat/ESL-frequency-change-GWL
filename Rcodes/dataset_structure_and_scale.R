
# Author: Khin Nawarat (2026, @IHE Delft)
# Dataset structure and spatial scale (Figure 2)

remove(list = ls())
library(ncdf4)
library(maps)

filedir <- "D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/"   # change to your local directory

ncfile_name <- paste0(filedir, "Data/GWL_VS_FC_larger_set.nc")
nc <- nc_open(ncfile_name)

# -----------------------------
# Read variables from nc file
# -----------------------------
lon      <- ncvar_get(nc, "longitude")
lat      <- ncvar_get(nc, "latitude")
locid    <- ncvar_get(nc, "location_id")
gwl_data <- ncvar_get(nc, "GWL")

nc_close(nc)

gridcoo.commonlocs <- data.frame(
  longitude = lon,
  latitude  = lat,
  location_id = locid
)

# -----------------------------
# Confirm dimensions
# GWL dims = [location, return_level, method, uncertainty]
# -----------------------------
print(dim(gwl_data))

# return_level indices
# 1=10, 2=20, 3=50, 4=100, 5=200
idx_RL_10  <- 1
idx_RL_200 <- 5

# method indices
# 1=VS, 2=FC
idx_FC <- 2

# uncertainty indices
# 1=q0.05, 2=q0.5, 3=q0.95
idx_q50 <- 2

# -----------------------------
# Extract data for figure
# -----------------------------
gwl_10  <- gwl_data[, idx_RL_10,  idx_FC, idx_q50]
gwl_200 <- gwl_data[, idx_RL_200, idx_FC, idx_q50]

# Replace missval
gwl_10[gwl_10 < 0]   <- NA
gwl_200[gwl_200 < 0] <- NA

# -----------------------------
# Labels and colours
# 1:6 = GWL categories
# 7   = None / no transition within considered range
# -----------------------------
gwl_labels <- c("1.5°C", "2°C", "2.5°C", "3°C", "4°C", "5°C", "None")

colpal <- c("#A11643", "#D73C4A", "#FEAF61","#E7F593", "#AADFA2", "#62C4A3", "#5E4DA3")

# Convert values to plotting colour indices
# valid stored values are 1:7
to_col_index <- function(x) {
  out <- rep(7, length(x))   # default to "None" colour
  out[!is.na(x) & x %in% 1:7] <- x[!is.na(x) & x %in% 1:7]
  return(out)
}

col_10  <- colpal[to_col_index(gwl_10)]
col_200 <- colpal[to_col_index(gwl_200)]

# -----------------------------
# Plot function
# -----------------------------
plot_panel <- function(lon, lat, cols, main_title) {
  par(fg = "gray70")
  
  plot(data.frame(longitude = lon, latitude = lat),
       type = "n",
       las = 1,
       xlab = "Longitude",
       ylab = "Latitude",
       ylim = c(-90, 90),
       xlim = c(-180, 180),
       axes = FALSE,
       cex.main = 1.5,
       main = main_title)
  
  axis(1, at = seq(-200, 200, by = 50))
  axis(2, at = seq(-100, 100, by = 20), las = 1)
  box()
  
  map("world", add = TRUE, interior = FALSE, fill = TRUE, col = "gray70")
  points(lon, lat, pch = 19, col = cols, cex = 0.9)
  
  par(fg = 1)
}
# -----------------------------
# Output figure
# -----------------------------
jpeg(
  paste0(filedir, "pics/data_structure_and_spatial_scale.jpg"),
  quality = 100,
  width   = 1200,
  height  = 450,
  res     = 100
)

# top row: two panels
# bottom row: one shared legend spanning both columns
layout(
  matrix(c(1, 2,
           3, 3), nrow = 2, byrow = TRUE),
  heights = c(12, 2),
  widths  = c(1, 1)
)

# panel (a)
par(mar = c(4, 4, 2.5, 1.5))
plot_panel(
  gridcoo.commonlocs$longitude,
  gridcoo.commonlocs$latitude,
  col_10,
  "(a)"
)

# panel (b)
par(mar = c(4, 4, 2.5, 1.5))
plot_panel(
  gridcoo.commonlocs$longitude,
  gridcoo.commonlocs$latitude,
  col_200,
  "(b)"
)

# legend 
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

xpos <- seq(0.34, 0.66, length.out = 7)
y_dots <- 0.55
y_text <- 0.20

points(xpos, rep(y_dots, 7), pch = 19, col = colpal, cex = 3)
text(xpos, rep(y_text, 7), labels = gwl_labels, cex = 1.2)
dev.off()

