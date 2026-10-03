# How far does the flow-turbidity relationship extend along the San Joaquin River, and is
# summer turbidity at RRI after the Golden Mussel detection unusual given upstream conditions?
#
# Stations, upstream to downstream: SJR (Vernalis), MSD (Mossdale), RRI (Stockton), TWI
# (Twitchell Island), ANH (Antioch). RRI comes from turbidity-depth-RRI.R (three depths);
# the others are single-sonde files in data/ and go through qc_single_sonde().
# Flow comes from CDEC, downloaded with
#   https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=<id>&SensorNums=<n>&dur_code=<D or E>&Start=2005-07-01&End=<today>
# and saved in data/:
#   VNS  Vernalis, daily mean (sensor 41)                    CDEC_VNS_flow_daily.csv
#   MSD  Mossdale, daily mean (sensor 41)                    CDEC_MSD_flow_daily.csv
#   SJJ  Jersey Point, 15-minute tidal flow (sensor 20),     CDEC_SJJ_flow_15min.csv
#        Godin-filtered to a daily net flow (just downstream of Twitchell Island)
# Antioch has stage but no flow gauge, so it is compared with SJJ, the closest flow record.
# All stations are compared with Vernalis flow (section 1, 3); section 1b and 3b use each
# station's own nearest flow.
#
# Vernalis is the only control (no mussels known there). Mossdale is NOT a control: CDFW reports
# Golden Mussels there, so a shift at MSD is expected to resemble one at RRI, not to be absent.
# To add a summer, extend the data files and rerun; years from `first_post_year` onward are
# treated as post-detection and compared with the earlier years.
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))
source(here("R", "flow_functions.R"))
if (!exists("wq_qc_final")) source(here("turbidity-depth-RRI.R"))

first_post_year <- 2025            # first summer after the detection (assumed Oct 2024; adjust if needed)
early_years <- 2008:2010           # dropped in the second baseline
river_order <- c("SJR", "MSD", "RRI", "TWI", "ANH")
min_valid_readings_per_day <- 48   # a day needs at least half its 15-minute readings
min_days_per_summer <- 60

# Single-sonde stations: conductivity ceiling is raised for the tidal, brackish stations
station_limits <- list(
  SJR = param_limits,
  MSD = param_limits,
  TWI = modifyList(param_limits, list(spcond_uscm = c(20, 30000))),
  ANH = modifyList(param_limits, list(spcond_uscm = c(20, 30000)))
)
stations <- imap(station_limits, ~ qc_single_sonde(here("data", paste0(.y, "_data.csv")), .x))

# Daily turbidity per station (median of the analysis column) ----------------------------------------
daily_median <- function(df, value) {
  df |>
    mutate(date = as_date(date_time)) |>
    group_by(date) |>
    summarise(turb = if (sum(!is.na(.data[[value]])) >= min_valid_readings_per_day)
                       median(.data[[value]], na.rm = TRUE) else NA_real_, .groups = "drop")
}
# RRI: median across its three depths, when at least two have a valid day
rri_daily <- wq_qc_final |>
  group_by(depth_m) |>
  group_modify(~ daily_median(.x, "turb_fnu_analysis")) |>
  group_by(date) |>
  summarise(turb = if (sum(!is.na(turb)) >= 2) median(turb, na.rm = TRUE) else NA_real_, .groups = "drop") |>
  mutate(station = "RRI")

daily <- bind_rows(map2_dfr(stations, names(stations), ~ mutate(daily_median(.x, "turb_fnu_analysis"), station = .y)),
                   rri_daily) |>
  filter(!is.na(turb), turb > 0) |>
  mutate(station = factor(station, levels = river_order))

# Vernalis flow (used for every station in sections 1-3). Log scale needs positive flow.
flow <- read_cdec_daily(here("data", "CDEC_VNS_flow_daily.csv")) |> filter(flow > 0)

# Each station's own nearest flow (sections 1b, 3b). `x` is the covariate used in the models:
# log10(flow) where flow is always positive (river stations), flow in 1000 cfs where it can
# be negative (tidally filtered flow at SJJ in summer is near zero, often below).
msd_flow <- read_cdec_daily(here("data", "CDEC_MSD_flow_daily.csv"))
sjj_flow <- godin_daily_flow(here("data", "CDEC_SJJ_flow_15min.csv"))
local_flow <- list(
  VNS = transmute(flow, date, x = log10(flow)),
  MSD = msd_flow |> filter(flow > 0) |> transmute(date, x = log10(flow)),
  SJJ = sjj_flow |> transmute(date, x = flow / 1000)
)
local_flow_of <- c(SJR = "VNS", MSD = "MSD", RRI = "VNS", TWI = "SJJ", ANH = "SJJ")
cat("Flow records (days with data) and summer days at or below zero:\n")
print(tibble(gauge = c("VNS", "MSD", "SJJ"),
             days = c(nrow(flow), nrow(msd_flow), nrow(sjj_flow)),
             summer_days_le_0 = c(0, sum(month(msd_flow$date) %in% 7:9 & msd_flow$flow <= 0),
                                  sum(month(sjj_flow$date) %in% 7:9 & sjj_flow$flow <= 0))))

cat("Valid QC share of readings at the single-sonde stations:\n")
imap_dfr(stations, ~ .x |> summarise(station = .y, pct_valid = round(100 * mean(qaqc_flag == "Valid"), 1))) |> print()

# 1. Flow-turbidity relationship by station --------------------------------------------------------------
# Daily turbidity against Vernalis flow from `lag` days earlier (Spearman).
lag_correlations <- expand_grid(station = river_order, lag = c(0, 2, 4, 6, 8, 10, 14)) |>
  pmap_dfr(function(station, lag) {
    x <- daily |> filter(station == !!station) |> inner_join(mutate(flow, date = date + lag), by = "date")
    summer <- month(x$date) %in% 7:9
    tibble(station, lag,
           all_days = cor(x$turb, x$flow, method = "spearman"),
           summer_days = cor(x$turb[summer], x$flow[summer], method = "spearman"))
  })
cat("\nSpearman correlation of daily turbidity with Vernalis flow, all days:\n")
lag_correlations |>
  select(station, lag, all_days) |>
  pivot_wider(names_from = lag, values_from = all_days, names_prefix = "lag") |>
  mutate(across(-station, ~ round(.x, 2))) |>
  print()
cat("\nJuly-September days only:\n")
lag_correlations |>
  select(station, lag, summer_days) |>
  pivot_wider(names_from = lag, values_from = summer_days, names_prefix = "lag") |>
  mutate(across(-station, ~ round(.x, 2))) |>
  print()

cat("\nSummer log-log regression at lag 0, log10(turbidity) ~ log10(flow):\n")
map_dfr(river_order, function(s) {
  x <- daily |> filter(station == s, month(date) %in% 7:9) |> inner_join(flow, by = "date")
  m <- lm(log10(turb) ~ log10(flow), x)
  tibble(station = s, days = nrow(x), slope = round(coef(m)[2], 2), r2 = round(summary(m)$r.squared, 2))
}) |> print()

# 1b. Each station against its own nearest flow gauge ---------------------------------------------------------------
# Spearman is unaffected by the log/linear choice of x. RRI and SJR use Vernalis, as above.
own_flow_lags <- expand_grid(station = river_order, lag = c(0, 1, 2, 4, 8)) |>
  pmap_dfr(function(station, lag) {
    g <- local_flow_of[[station]]
    x <- daily |> filter(station == !!station) |> inner_join(mutate(local_flow[[g]], date = date + lag), by = "date")
    summer <- month(x$date) %in% 7:9
    tibble(station, gauge = g, lag,
           all_days = cor(x$turb, x$x, method = "spearman"),
           summer_days = cor(x$turb[summer], x$x[summer], method = "spearman"))
  })
cat("\nSpearman correlation of daily turbidity with each station's own flow gauge (lag in days):\n")
own_flow_lags |>
  pivot_longer(c(all_days, summer_days), names_to = "days", values_to = "rho") |>
  mutate(rho = round(rho, 2)) |>
  pivot_wider(names_from = lag, values_from = rho, names_prefix = "lag") |>
  arrange(factor(station, levels = river_order), days) |>
  print(n = Inf)

cat("\nSummer log10(turbidity) ~ own flow at lag 0 (x = log10 cfs, or 1000 cfs for SJJ):\n")
map_dfr(river_order, function(s) {
  g <- local_flow_of[[s]]
  x <- daily |> filter(station == s, month(date) %in% 7:9) |> inner_join(local_flow[[g]], by = "date")
  m <- lm(log10(turb) ~ x, x)
  tibble(station = s, gauge = g, days = nrow(x), slope = round(coef(m)[2], 3), r2 = round(summary(m)$r.squared, 2))
}) |> print()

# Daily tidal range at SJJ (tidal energy; spring-neap) against turbidity at the two tidal stations
cat("\nSummer Spearman correlation of turbidity with SJJ daily tidal range (lag 0):\n")
map_dfr(c("TWI", "ANH"), function(s) {
  x <- daily |> filter(station == s, month(date) %in% 7:9) |> inner_join(sjj_flow, by = "date")
  tibble(station = s, days = nrow(x), rho = round(cor(x$turb, x$tidal_range, method = "spearman"), 2))
}) |> print()

# 2. Summer medians along the river and downstream/upstream ratios ------------------------------------------------
summer_medians <- daily |>
  filter(month(date) %in% 7:9) |>
  group_by(station, year = year(date)) |>
  summarise(days = n(), median = median(turb), .groups = "drop") |>
  filter(days >= min_days_per_summer)

cat("\nSummer median turbidity (FNU) by station and year:\n")
summer_medians |>
  select(-days) |>
  mutate(median = round(median, 1)) |>
  pivot_wider(names_from = station, values_from = median) |>
  select(year, any_of(river_order)) |>
  arrange(year) |>
  print(n = Inf)

summer_wide <- summer_medians |> select(-days) |> pivot_wider(names_from = station, values_from = median)
station_pairs <- list(c("RRI", "MSD"), c("RRI", "SJR"), c("TWI", "MSD"), c("ANH", "MSD"), c("MSD", "SJR"))
ratio_check <- expand_grid(post_year = unique(summer_medians$year[summer_medians$year >= first_post_year]),
                           pair = station_pairs) |>
  pmap_dfr(function(post_year, pair) {
    r <- summer_wide |> transmute(year, ratio = .data[[pair[1]]] / .data[[pair[2]]]) |> filter(is.finite(ratio))
    base <- filter(r, year < first_post_year)
    post <- filter(r, year == post_year)
    tibble(post_year, pair = paste(pair[1], "/", pair[2]), baseline_years = nrow(base),
           baseline_median = round(median(base$ratio), 2),
           baseline_range = paste(round(min(base$ratio), 2), "-", round(max(base$ratio), 2)),
           post_ratio = round(post$ratio, 2),
           pct_baseline_lower = round(100 * mean(base$ratio < post$ratio)))
  })
cat("\nSummer downstream/upstream turbidity ratio, post-detection summers vs earlier years:\n")
print(ratio_check, n = Inf, width = Inf)

# 3. Prediction from Vernalis turbidity and flow ----------------------------------------------------------------------
# log10(turbidity at a station) ~ log10(Vernalis turbidity, lagged) + flow covariate, summer days.
# The flow covariate is log10 Vernalis flow here; section 3b swaps in each station's own gauge.
# Each baseline year is predicted from a model fit to the other baseline years (leave-one-year-out),
# so baseline residuals are as honest as the post-detection ones.
vernalis_daily <- daily |> filter(station == "SJR") |> select(date, vernalis_turb = turb)

model_data <- function(target, lag, gauge = "VNS") {
  daily |>
    filter(station == target, month(date) %in% 7:9) |>
    inner_join(mutate(vernalis_daily, date = date + lag), by = "date") |>
    inner_join(local_flow[[gauge]], by = "date") |>
    transmute(year = year(date), ly = log10(turb), lx = log10(vernalis_turb), lf = x) |>
    filter(is.finite(ly), is.finite(lx), is.finite(lf))
}

# Lag that fits RRI best over the baseline years (the choice barely matters: R-squared is flat)
lag_fit <- map_dfr(0:8, function(lag) {
  d <- model_data("RRI", lag) |> filter(year < first_post_year)
  tibble(lag, r2 = summary(lm(ly ~ lx + lf, d))$r.squared)
})
best_lag <- lag_fit$lag[which.max(lag_fit$r2)]
cat("\nVernalis lag used (days):", best_lag, "\n")

residual_test <- function(target, baseline_years, post_year, gauge = "VNS") {
  d <- model_data(target, best_lag, gauge)
  base_resid <- map_dfr(baseline_years, function(y) {
    test <- filter(d, year == y)
    if (nrow(test) < 40) return(NULL)
    fit <- lm(ly ~ lx + lf, filter(d, year %in% setdiff(baseline_years, y)))
    tibble(year = y, mean_resid = mean(test$ly - predict(fit, test)))
  })
  fit <- lm(ly ~ lx + lf, filter(d, year %in% baseline_years))
  test <- filter(d, year == post_year)
  post_resid <- mean(test$ly - predict(fit, test))
  tibble(station = target, flow_gauge = gauge, post_year, baseline_years = nrow(base_resid),
         sd_resid_log10 = round(sd(base_resid$mean_resid), 3),
         resid_log10 = round(post_resid, 3),
         pct_vs_expected = round(100 * (10^post_resid - 1)),
         pct_baseline_lower = round(100 * mean(base_resid$mean_resid < post_resid)),
         z = round((post_resid - mean(base_resid$mean_resid)) / sd(base_resid$mean_resid), 2),
         model_r2 = round(summary(fit)$r.squared, 2))
}

post_years <- sort(unique(summer_medians$year[summer_medians$year >= first_post_year]))
baselines <- list("all earlier years" = function(s) setdiff(sort(unique(model_data(s, best_lag)$year)), post_years),
                  "without 2008-2010" = function(s) setdiff(sort(unique(model_data(s, best_lag)$year)), c(post_years, early_years)))
targets <- c("RRI", "MSD", "ANH", "TWI")   # MSD has mussels (CDFW), so it is a second affected station, not a placebo
prediction_test <- expand_grid(target = targets, baseline = names(baselines), post_year = post_years) |>
  pmap_dfr(function(target, baseline, post_year) {
    mutate(residual_test(target, baselines[[baseline]](target), post_year), baseline = baseline)
  })
cat("\nPost-detection summers vs prediction from Vernalis turbidity and flow:\n")
print(prediction_test, n = Inf, width = Inf)
# 3b. The same test using each station's own nearest flow gauge instead of Vernalis flow
own_targets <- tibble(target = c("MSD", "TWI", "ANH"), gauge = c("MSD", "SJJ", "SJJ"))
own_flow_test <- expand_grid(own_targets, baseline = names(baselines), post_year = post_years) |>
  pmap_dfr(function(target, gauge, baseline, post_year) {
    mutate(residual_test(target, baselines[[baseline]](target), post_year, gauge), baseline = baseline)
  })
cat("\nSame test with each station's own flow gauge (MSD: Mossdale; TWI, ANH: Jersey Point net flow):\n")
print(own_flow_test, n = Inf, width = Inf)
cat("\nMinimum departure from expectation that would stand out at 2 standard deviations (RRI, all earlier years):",
    round(100 * (1 - 10^(-2 * filter(prediction_test, station == "RRI", baseline == "all earlier years")$sd_resid_log10[1]))),
    "% lower for one summer\n")
