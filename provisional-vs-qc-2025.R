# How different are CDEC's real-time (not QC'd) turbidity readings from the QC'd station files?
# 2025 is the one year where both exist for every station, so it shows how much the provisional
# 2026 numbers could differ from what the finished record will say.
#
# 1. Raw readings at matching 15-minute timestamps: how often identical, ratio when not, gaps.
# 2. After the same QC: daily and July-September medians from the CDEC download vs the station file.
#
# Needs data/CDEC_<station>_turb_2025.csv (CDEC CSVDataServlet, sensor 27 NTU or 221 FNU, 2025-01-01
# to 2026-01-01). The QC'd side is the station files, or for RRI the 1 m sonde record. The station-file
# daily medians come from flow-turbidity-network.R (sourced if `daily` is not in memory).
# Leaves `provisional_ratio` (CDEC / QC'd summer 2025 median, per station) for later use.
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))
if (!exists("daily")) source(here("flow-turbidity-network.R"))
if (!exists("wq_qc_final")) source(here("turbidity-depth-RRI.R"))

stations <- c("SJR", "MSD", "RRI", "PPT", "TWI", "BLP", "ANH", "TRN", "OSJ", "MOK")
tidal_cond <- modifyList(param_limits, list(spcond_uscm = c(20, 30000)))
limits_for <- function(s) if (s %in% c("TWI", "ANH", "BLP")) tidal_cond else param_limits
in_2025 <- function(x) filter(x, year(date_time) == 2025)

qcd_raw <- function(s) {
  x <- if (s == "RRI") {
    wq_qc_final |> filter(depth_m == 1) |> transmute(date_time, turb = turb_fnu)
  } else {
    read_station_file(here("data", paste0(s, "_data.csv"))) |> transmute(date_time, turb = turb_fnu)
  }
  in_2025(x)
}
cdec_raw <- function(s) {
  read_station_file(here("data", paste0("CDEC_", s, "_turb_2025.csv"))) |> transmute(date_time, turb = turb_fnu) |> in_2025()
}

# 1. Raw agreement ---------------------------------------------------------------------------------------------
agreement <- map_dfr(stations, function(s) {
  j <- full_join(rename(qcd_raw(s), qcd = turb), rename(cdec_raw(s), cdec = turb), by = "date_time")
  both <- filter(j, !is.na(qcd), !is.na(cdec))
  tibble(station = s, n_both = nrow(both), only_qcd = sum(!is.na(j$qcd) & is.na(j$cdec)),
         only_cdec = sum(is.na(j$qcd) & !is.na(j$cdec)),
         pct_identical = round(100 * mean(abs(both$cdec - both$qcd) <= 0.005), 1),
         pct_differ_10pct = round(100 * mean(abs(both$cdec - both$qcd) > 0.1 * pmax(both$qcd, 0.5)), 1),
         median_ratio = round(median(both$cdec[both$qcd > 0.5] / both$qcd[both$qcd > 0.5]), 3),
         spearman = round(cor(both$cdec, both$qcd, method = "spearman"), 3),
         qcd_p99 = round(quantile(both$qcd, 0.99), 1), cdec_p99 = round(quantile(both$cdec, 0.99), 1),
         qcd_max = round(max(both$qcd), 1), cdec_max = round(max(both$cdec), 1))
})
cat("Raw readings, 2025, matching timestamps (CDEC real-time vs QC'd file):\n")
print(agreement, width = Inf)

# 2. Same QC on both: daily and summer medians -------------------------------------------------------------------
daily_median <- function(df, value) {
  df |> mutate(date = as_date(date_time)) |> group_by(date) |>
    summarise(turb = if (sum(!is.na(.data[[value]])) >= 48) median(.data[[value]], na.rm = TRUE) else NA_real_, .groups = "drop")
}
cdec_daily <- map_dfr(stations, function(s) {
  q <- qc_station_data(read_station_file(here("data", paste0("CDEC_", s, "_turb_2025.csv"))) |> in_2025(), limits_for(s))
  if (s == "RRI") q <- prepare_analysis_turbidity(q)
  mutate(daily_median(q, "turb_fnu_analysis"), station = s)
})
qcd_daily <- daily |> filter(year(date) == 2025) |> mutate(station = as.character(station))
both_daily <- inner_join(rename(qcd_daily, qcd = turb), rename(cdec_daily, cdec = turb), by = c("station", "date")) |>
  filter(!is.na(qcd), !is.na(cdec), qcd > 0, cdec > 0)

summer_cmp <- both_daily |>
  mutate(summer = month(date) %in% 7:9) |>
  group_by(station) |>
  summarise(days = n(),
            all_year_median_ratio = round(median(cdec / qcd), 3),
            summer_days = sum(summer),
            summer_median_qcd = round(median(qcd[summer]), 2), summer_median_cdec = round(median(cdec[summer]), 2),
            summer_ratio = round(median(cdec[summer]) / median(qcd[summer]), 3),
            pct_days_differ_10pct = round(100 * mean(abs(cdec / qcd - 1) > 0.1)),
            .groups = "drop") |>
  mutate(station = factor(station, levels = stations)) |>
  arrange(station)
cat("\nAfter the same QC, daily medians: CDEC vs QC'd file, 2025\n")
print(summer_cmp, width = Inf)

cat("\nMonthly median of the daily CDEC / QC'd ratio (is the difference seasonal?):\n")
both_daily |>
  group_by(station, month = month(date)) |>
  summarise(ratio = round(median(cdec / qcd), 2), .groups = "drop") |>
  mutate(station = factor(station, levels = stations)) |>
  pivot_wider(names_from = month, values_from = ratio) |>
  arrange(station) |>
  print(width = Inf)

provisional_ratio <- summer_cmp |> transmute(station = as.character(station), ratio = summer_ratio)
