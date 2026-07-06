# Author: Khin Nawarat (2026, @IHE Delft)
# Check internal consistency and reproduce Figure 3

remove(list = ls())
library(ncdf4)

filedir <- "D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # Change to your local directory
ncfile_name <- paste0(filedir, "Data/GWL_VS_FC_larger_set.nc")

nc <- nc_open(ncfile_name)
gwl_data <- ncvar_get(nc, "GWL")
nc_close(nc)

# Dimensions
RL_vals <- c(10,20,50,100,200)
AEP <- 1 / RL_vals

idx_q50 <- 2
idx_VS  <- 1
idx_FC  <- 2

# GWL indices for 1.5–5°C
gwl_indices <- 1:6
gwl_labels  <- c("1.5°C","2°C","2.5°C","3°C","4°C","5°C")

colpal <- c(
  "#A11643", "#D73C4A", "#FEAF61",
  "#E7F593", "#AADFA2", "#62C4A3"
)

# Function: cumulative percentage (<= GWL)
compute_cum_percent <- function(method_idx) {
  
  mat <- matrix(NA, nrow = length(gwl_indices), ncol = length(RL_vals))
  
  for (g in seq_along(gwl_indices)) {
    for (r in seq_along(RL_vals)) {
      
      vec <- gwl_data[, r, method_idx, idx_q50]
      vec[vec == -9999] <- NA
      
      mat[g, r] <- mean(vec <= gwl_indices[g], na.rm = TRUE) * 100
    }
  }
  
  return(mat)
}

VS_mat <- compute_cum_percent(idx_VS)
FC_mat <- compute_cum_percent(idx_FC)

# ========================
# Plotting
# ========================
par(mfrow = c(1,2),
    mar = c(5.5,5.5,3,2))

cex_axis  <- 1.3
cex_lab   <- 1.4
cex_leg   <- 1.3
cex_title <- 1.5
lwd_curv  <- 2

# ---- VS panel ----
plot(RL_vals, VS_mat[1,], type = "b",
     col = colpal[1], pch = 16, lwd = lwd_curv,
     ylim = c(0,100), yaxs = "i",
     xaxt = "n", log = "x",            
     xlab = "Return period (years)",
     ylab = "Cumulative percentage of \nlocations becoming annual (%)",
     main = "(a)",
     cex.main = cex_title,
     cex.lab = cex_lab,
     cex.axis = cex_axis)

for (g in 2:6) {
  lines(RL_vals, VS_mat[g,],
        type = "b", col = colpal[g],
        pch = 16, lwd = lwd_curv)
}

axis(1, at = RL_vals, labels = RL_vals, cex.axis = cex_axis)


# ---- FC panel ----
plot(RL_vals, FC_mat[1,], type = "b",
     col = colpal[1], pch = 16, lwd = lwd_curv,
     ylim = c(0,100), yaxs = "i",
     xaxt = "n", log = "x",
     xlab = "Return period (years)",
     ylab = "",
     main = "(b)",
     cex.main = cex_title,
     cex.lab = cex_lab,
     cex.axis = cex_axis)

for (g in 2:6) {
  lines(RL_vals, FC_mat[g,],
        type = "b", col = colpal[g],
        pch = 16, lwd = lwd_curv)
}

axis(1, at = RL_vals, labels = RL_vals, cex.axis = cex_axis)

# ---- Legend ----
legend("bottomright",
       legend = gwl_labels,
       col = colpal,
       title = "Global warming level",
       lty = 1,lwd = lwd_curv,cex = cex_leg,
       pch = 16,
       bty = "n")

