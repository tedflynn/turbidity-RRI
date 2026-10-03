# How far does the flow-turbidity relationship extend along the San Joaquin River, and is
# summer turbidity at RRI after the Golden Mussel detection unusual given upstream conditions?
#
# Stations, upstream to downstream: SJR (Vernalis), MSD (Mossdale), RRI (Stockton), TWI
# (Twitchell Island), ANH (Antioch). RRI comes from turbidity-depth-RRI.R (three depths);
# the others are single-sonde files in data/ and go through qc_single_sonde().
# Flow is the daily mean at Vernalis (CDEC VNS, sensor 41), downloaded from
#   https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=VNS&SensorNums=41&dur_code=D&Start=2005-07-01&End=<today>
# and saved as data/CDEC_VNS_flow_daily.csv.
#
# Vernalis is the control (no mussels known there); whether Mossdale is mussel-free is unconfirmed.
# To add a summer, extend the data files and rerun; years from `first_post_year` onward are
# treated as post-detection and compared with the earlier years.
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))
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

flow <- read_csv(here("data", "CDEC_VNS_flow_daily.csv"), show_col_types = FALSE,
                 col_types = cols(.default = col_character())) |>
  transmute(date = as_date(ymd_hm(`DATE TIME`)), flow = suppressWarnings(as.numeric(VALUE))) |>
  filter(!is.na(flow), flow > 0)

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
# log10(turbidity at a station) ~ log10(Vernalis turbidity, lagged) + log10(Vernalis flow), summer days.
# Each baseline year is predicted from a model fit to the other baseline years (leave-one-year-out),
# so baseline residuals are as honest as the post-detection ones.
vernalis_daily <- daily |> filter(station == "SJR") |> select(date, vernalis_turb = turb)

model_data <- function(target, lag) {
  daily |>
    filter(station == target, month(date) %in% 7:9) |>
    inner_join(mutate(vernalis_daily, date = date + lag), by = "date") |>
    inner_join(flow, by = "date") |>
    transmute(year = year(date), ly = log10(turb), lx = log10(vernalis_turb), lf = log10(flow)) |>
    filter(is.finite(ly), is.finite(lx), is.finite(lf))
}

# Lag that fits RRI best over the baseline years (the choice barely matters: R-squared is flat)
lag_fit <- map_dfr(0:8, function(lag) {
  d <- model_data("RRI", lag) |> filter(year < first_post_year)
  tibble(lag, r2 = summary(lm(ly ~ lx + lf, d))$r.squared)
})
best_lag <- lag_fit$lag[which.max(lag_fit$r2)]
cat("\nVernalis lag used (days):", best_lag, "\n")

residual_test <- function(target, baseline_years, post_year) {
  d <- model_data(target, best_lag)
  base_resid <- map_dfr(baseline_years, function(y) {
    test <- filter(d, year == y)
    if (nrow(test) < 40) return(NULL)
    fit <- lm(ly ~ lx + lf, filter(d, year %in% setdiff(baseline_years, y)))
    tibble(year = y, mean_resid = mean(test$ly - predict(fit, test)))
  })
  fit <- lm(ly ~ lx + lf, filter(d, year %in% baseline_years))
  test <- filter(d, year == post_year)
  post_resid <- mean(test$ly - predict(fit, test))
  tibble(station = target, post_year, baseline_years = nrow(base_resid),
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
targets <- c("RRI", "MSD", "ANH", "TWI")   # MSD is a placebo: it should not shift if it is mussel-free
prediction_test <- expand_grid(target = targets, baseline = names(baselines), post_year = post_years) |>
  pmap_dfr(function(target, baseline, post_year) {
    mutate(residual_test(target, baselines[[baseline]](target), post_year), baseline = baseline)
  })
cat("\nPost-detection summers vs prediction from Vernalis turbidity and flow:\n")
print(prediction_test, n = Inf, width = Inf)
cat("\nMinimum departure from expectation that would stand out at 2 standard deviations (RRI, all earlier years):",
    round(100 * (1 - 10^(-2 * filter(prediction_test, station == "RRI", baseline == "all earlier years")$sd_resid_log10[1]))),
    "% lower for one summer\n")
