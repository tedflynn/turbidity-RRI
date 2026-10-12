# Inventory of the raw EDI continuous water quality downloads (raw-EDI-data/<program>/<station>_*.csv).
#
# For each file: rows, period, what each parameter covers, where turbidity starts and ends, how many
# summers (Jul-Sep) have enough turbidity to use, and a clock check. The clock check matters because the
# analysis assumes Pacific Standard Time: in summer, water temperature peaks in mid to late afternoon,
# so the median hour of the daily maximum should fall around 15-17 h if the clock is PST. A value near
# 16 + 8 = 24 would mean UTC, and a value an hour later than expected would mean daylight time.
#
# Output: data/edi_inventory.csv (one row per file) and data/edi_summer_days.csv (valid summer days per
# station and year). These count raw readings before any QC; qc_jobs / run-qc.R does the QC.
library(tidyverse)
library(lubridate)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))

min_readings_per_day <- 48      # same rule as the daily medians in the analysis
files <- list.files(here("raw-EDI-data"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)

profile_one <- function(path) {
  x <- read_station_file(path)
  station <- x$station[1]
  has <- function(v) !is.na(x[[v]])
  tz_hour <- x |>
    filter(month %in% 6:9, has("watertemp_c")) |>
    mutate(date = as_date(date_time), hour = hour(date_time) + minute(date_time) / 60) |>
    group_by(date) |>
    filter(sum(!is.na(watertemp_c)) >= 80) |>             # near-complete days only
    slice_max(watertemp_c, n = 1, with_ties = FALSE) |>
    ungroup() |>
    summarise(n = n(), hour = median(hour))
  turb <- x[has("turb_fnu"), ]
  daily <- x |>
    filter(month %in% 7:9) |>
    group_by(year, date = as_date(date_time)) |>
    summarise(n = sum(!is.na(turb_fnu)), .groups = "drop") |>
    filter(n >= min_readings_per_day) |>
    count(year, name = "valid_summer_days") |>
    mutate(station = station, .before = 1)
  inv <- tibble(
    station = station, program = basename(dirname(path)), file = basename(path),
    rows = nrow(x), first = as_date(min(x$date_time)), last = as_date(max(x$date_time)),
    turb_first = if (nrow(turb)) as_date(min(turb$date_time)) else NA_Date_,
    turb_last = if (nrow(turb)) as_date(max(turb$date_time)) else NA_Date_,
    n_turbidity = nrow(turb),
    pct_turbidity = round(100 * mean(has("turb_fnu")), 1),
    pct_temp = round(100 * mean(has("watertemp_c")), 1), pct_cond = round(100 * mean(has("spcond_uscm")), 1),
    pct_do = round(100 * mean(has("do_mgl")), 1), pct_ph = round(100 * mean(has("ph")), 1),
    pct_chl = round(100 * mean(has("chlor_ugl")), 1),
    temp_peak_hour = round(tz_hour$hour, 1), temp_peak_days = tz_hour$n,
    median_turbidity = round(median(turb$turb_fnu), 1), p99_turbidity = round(quantile(turb$turb_fnu, 0.99, names = FALSE), 1),
    duplicate_times = sum(duplicated(x$date_time)))
  list(inventory = inv, summer = daily)
}

res <- map(files, profile_one)
inventory <- map_dfr(res, "inventory") |> arrange(program, station)
summer_days <- map_dfr(res, "summer")

write_csv(inventory, here("data", "edi_inventory.csv"))
write_csv(summer_days, here("data", "edi_summer_days.csv"))

cat("Files:", nrow(inventory), "\n\nPeriod and turbidity coverage:\n")
inventory |> select(program, station, first, last, turb_first, turb_last, pct_turbidity, median_turbidity, p99_turbidity) |>
  print(n = Inf, width = Inf)
cat("\nShare of rows with each parameter (%):\n")
inventory |> select(station, pct_temp, pct_cond, pct_do, pct_ph, pct_chl, pct_turbidity) |> print(n = Inf)
cat("\nClock check: median hour of the daily maximum water temperature in Jun-Sep (expect about 15-17 if PST):\n")
inventory |> select(station, temp_peak_hour, temp_peak_days) |> arrange(temp_peak_hour) |> print(n = Inf)
cat("\nDuplicate timestamps:\n"); print(filter(inventory, duplicate_times > 0) |> select(station, duplicate_times))
