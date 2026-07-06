# Author: Khin Nawarat  (2026, @IHE Delft)
# Different Monte Carlo realizations for 200-year return levels

remove(list = ls())

filedir <- "D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/"
output_dir <- paste0(filedir,"Rdatasets")
set.seed(42)

library(dplyr)
library(tidyr)
library(ggplot2)

# Source existing functions from your repository
source(paste0(filedir, "Rcode/allfunctions.r"))

# Read ESL datasets
rl1r <- read.csv(file.path(filedir, "CSV/Rasmussen/Rasmussen_ESLs.csv"))
rl1k <- read.csv(file.path(filedir, "CSV/Kirezci/Kirezci_ESLs.csv"))
rl1v <- read.csv(file.path(filedir, "CSV/Vousdoukas/Vousdoukas_ESLs.csv"))
rsites.gpd.parameters <- read.delim(file.path(filedir, "CSV/Rasmussen/gpd_parameters_uhawaii19.tsv"))


build_site_list <- function(data, dataset_name) {
  valid <- data %>% filter(shape > -0.5)
  if (nrow(valid) == 0) stop(sprintf("No valid sites for %s", dataset_name))
  valid %>% mutate(dataset = dataset_name)
}

sites_r <- build_site_list(rl1r, "Rasmussen")
sites_k <- build_site_list(rl1k, "Kirezci")
sites_v <- build_site_list(rl1v, "Vousdoukas")

# Set return level
RL <- 200

# Sample sizes to compare
grid_n <- c(1000, 5000, 10000)

# Compute site-level year length for Rasmussen
get_years <- function(site_id) {
  mm <- match(site_id, rsites.gpd.parameters$PSMSL_ID)
  if (is.na(mm)) return(NA_integer_)
  params <- rsites.gpd.parameters[mm, ]
  return(params$GPD_Record_End - params$GPD_Record_Start + 1)
}

sites_r <- sites_r %>% rowwise() %>% mutate(years = get_years(SLRID)) %>% ungroup()

compute_summary <- function(row, nsample) {
  if (row$dataset == "Rasmussen") {
    years <- row$years
  } else if (row$dataset == "Kirezci") {
    years <- 35
  } else if (row$dataset == "Vousdoukas") {
    years <- 170
  } else {
    stop("Unknown dataset")
  }
  if (is.na(years) || years <= 0) return(NULL)

  sim <- simulate.bivariate(nsample, row$scale, row$shape, row$lambda, years)
  if (any(is.na(sim))) return(NULL)

  TWL200 <- RL.N(RL, sim[, 2], sim[, 1], row$thresh, row$lambda)
  data.frame(
    dataset = row$dataset,
    SLRID = row$SLRID,
    nsample = nsample,
    mean = mean(TWL200, na.rm = TRUE),
    sd = sd(TWL200, na.rm = TRUE),
    p05 = quantile(TWL200, probs = 0.05, na.rm = TRUE),
    p50 = quantile(TWL200, probs = 0.50, na.rm = TRUE),
    p95 = quantile(TWL200, probs = 0.95, na.rm = TRUE),
    p99 = quantile(TWL200, probs = 0.99, na.rm = TRUE)
  )
}

# Build data frame for the selected site subset
all_selected <- bind_rows(sites_r, sites_k, sites_v) %>%
  dplyr::select(dataset, SLRID, longitude, latitude, shape, scale, lambda, thresh, years)

results <- do.call(rbind, lapply(seq_len(nrow(all_selected)), function(idx) {
  row <- all_selected[idx, ]
  do.call(rbind, lapply(grid_n, function(n) compute_summary(row, n)))
}))

results <- results %>% distinct(dataset, SLRID, nsample, .keep_all = TRUE)

if (is.null(results) || nrow(results) == 0) stop("No simulation results generated")

# Use the 10000 sample result as reference
reference <- results %>% filter(nsample == 10000) %>% dplyr::select(dataset, SLRID, p05_ref = p05, p50_ref = p50, p95_ref = p95, mean_ref = mean)

results <- results %>%
  left_join(reference, by = c("dataset", "SLRID")) %>%
  mutate(
    diff_mean = mean - mean_ref,
    abs_diff_mean = abs(diff_mean),
    diff_p05 = p05 - p05_ref,
    abs_diff_p05 = abs(diff_p05),
    diff_p50 = p50 - p50_ref,
    abs_diff_p50 = abs(diff_p50),
    diff_p95 = p95 - p95_ref,
    abs_diff_p95 = abs(diff_p95)
  )

# Write summary output
summary_csv <- file.path(output_dir, "mc_convergence_summary.csv")
write.csv(results, summary_csv, row.names = FALSE)
