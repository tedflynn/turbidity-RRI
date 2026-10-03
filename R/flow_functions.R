# Readers for CDEC flow downloads (data/CDEC_*.csv from the CSVDataServlet).
#
# read_cdec_daily(): daily mean flow (sensor 41). CDEC writes "---" for missing days, and a
#   few sentinel values (-99999) occur; anything below -10,000 cfs is treated as missing.
#   Negative daily means are real at tidal stations in low flow and are kept.
# godin_daily_flow(): net (tidally filtered) daily flow from 15-minute tidal flow (sensor 20).
#   Applies a Godin filter (boxcar of 24 h, 24 h, 25 h) to remove the tides, after treating
#   readings beyond `max_abs` cfs as bad. Gaps up to 6 h are interpolated before filtering,
#   and a day needs at least 80 filtered values and 70% real readings to count. Also returns
#   the day's tidal range (max - min of the raw flow), a proxy for tidal energy.
read_cdec_daily <- function(path) {
  read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
    transmute(date = as_date(ymd_hm(`DATE TIME`)), flow = suppressWarnings(as.numeric(VALUE))) |>
    filter(!is.na(date), !is.na(flow), flow > -10000)
}

godin_daily_flow <- function(path, max_abs = 250000) {
  boxcar <- function(v, n) as.numeric(stats::filter(v, rep(1 / n, n), sides = 2))
  godin <- function(x) boxcar(boxcar(boxcar(x, 96), 96), 100)

  d <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character())) |>
    transmute(t = ymd_hm(`DATE TIME`, tz = "Etc/GMT+8"), value = suppressWarnings(as.numeric(VALUE))) |>
    filter(!is.na(t)) |>
    distinct(t, .keep_all = TRUE) |>
    arrange(t)
  grid <- tibble(t = seq(min(d$t), max(d$t), by = "15 min")) |> left_join(d, by = "t")
  grid$value[!is.finite(grid$value) | abs(grid$value) > max_abs] <- NA

  grid |>
    mutate(filtered = godin(imputeTS::na_interpolation(value, maxgap = 24)), real = !is.na(value)) |>
    group_by(date = as_date(t)) |>
    summarise(flow = mean(filtered, na.rm = TRUE), real_pct = mean(real), n_filtered = sum(!is.na(filtered)),
              tidal_range = if (any(real)) max(value, na.rm = TRUE) - min(value, na.rm = TRUE) else NA_real_,
              .groups = "drop") |>
    mutate(flow = if_else(n_filtered >= 80 & real_pct >= 0.7, flow, NA_real_)) |>
    filter(!is.na(flow)) |>
    select(date, flow, tidal_range)
}
