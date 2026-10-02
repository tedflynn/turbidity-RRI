# Relate the QC flags to the field visits at which sondes were exchanged.
# 1. Do flagged episodes start/end at exchanges, and how many visits do they span?
#    (A fault that survives several exchanges points to the sonde position or equipment,
#     not to ordinary fouling, which an exchange resets.)
# 2. How large is the turbidity step across each exchange, and have we missed any large ones?
#    Compared with the field sheet's own dirty vs clean readings.
#
# Needs data/RRI_fieldsheet_merged.csv (field visit sheet; gitignored) and the QC output
# from turbidity-depth-RRI.R (sourced below if not already in memory).
library(tidyverse)
library(lubridate)
library(here)

if (!exists("wq_qc_final")) source(here("turbidity-depth-RRI.R"))

tz <- "Etc/GMT+8"
set.seed(7)

# Field visits ------------------------------------------------------------------------
# Field sheet times carry a "Z" suffix but are local clock time (PST): visits fall between 7 and 13 h.
# Where departure is missing, assume a 75-minute visit (the median).
visits <- read_csv(here("data", "RRI_fieldsheet_merged.csv"), show_col_types = FALSE,
                   col_types = cols(.default = col_character())) |>
  mutate(arr = force_tz(ymd_hms(ArrivalTime, quiet = TRUE), tz),
         dep = force_tz(ymd_hms(DepartureTime, quiet = TRUE), tz)) |>
  filter(arr >= min(wq_qc_final$date_time), arr <= max(wq_qc_final$date_time)) |>
  arrange(arr) |>
  mutate(dep = if_else(is.na(dep) | dep < arr, arr + 75 * 60, dep))
visit_times <- as.numeric(visits$arr)
cat("Field visits within the sonde record:", nrow(visits), "\n")

# 1. Flagged episodes vs visits ------------------------------------------------------------
# One flag column per source; an episode is a run of flagged readings with gaps under 24 h.
flag_sources <- c(turb_drift = "drift_flag", turb_erratic = "erratic_flag", turb_flat = "flat_flag",
                  ph = "ph_drift_flag", temp = "temp_drift_flag", cond = "cond_drift_flag",
                  do = "do_drift_flag", chl = "chl_drift_flag")
flagged_values <- c("Fail: Depth Drift", "Fail: Erratic Period", "Fail: Flat Line")

episodes <- imap_dfr(flag_sources, function(col, source_name) {
  wq_qc_final |>
    filter(.data[[col]] %in% flagged_values) |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(episode = cumsum(c(TRUE, diff(as.numeric(date_time)) > 24 * 3600))) |>
    group_by(depth_m, episode) |>
    summarise(start = min(date_time), end = max(date_time), .groups = "drop") |>
    mutate(source = source_name)
}) |>
  mutate(days = as.numeric(difftime(end, start, units = "days"))) |>
  filter(days >= 1)

# Whole visits inside an episode, ignoring 36 h at each end (flag buffers and boundary effects)
trim <- 36 * 3600
episodes <- episodes |>
  mutate(visits_inside = map2_int(start, end, ~ sum(visit_times > as.numeric(.x) + trim &
                                                      visit_times < as.numeric(.y) - trim)))

cat("\nEpisodes of 1+ day, by source:\n")
episodes |>
  group_by(source) |>
  summarise(episodes = n(), median_days = round(median(days), 1),
            share_spanning_1_visit = round(mean(visits_inside >= 1), 2),
            share_spanning_3_visits = round(mean(visits_inside >= 3), 2)) |>
  print(n = Inf)

cat("\nLongest episodes and the number of whole visits inside them:\n")
episodes |>
  arrange(desc(days)) |>
  transmute(source, depth_m, start = as_date(start), end = as_date(end),
            days = round(days, 1), visits_inside) |>
  head(14) |>
  print(n = Inf)

# Are episode starts/ends near a visit more often than chance? Null: shift the whole visit
# schedule by a random 8-60 days (either direction), which keeps its spacing.
near_a_visit <- function(times, visit_times_) {
  mean(map_lgl(times, ~ any(abs(as.numeric(.x) - visit_times_) < 36 * 3600)))
}
null <- map_dfr(1:200, function(i) {
  shift <- sample(c(-1, 1), 1) * runif(1, 8, 60) * 86400
  tibble(end = near_a_visit(episodes$end, visit_times + shift),
         start = near_a_visit(episodes$start, visit_times + shift))
})
cat(sprintf("\nEpisode ENDS within 36 h of a visit: observed %.2f; chance mean %.2f, 95th percentile %.2f\n",
            near_a_visit(episodes$end, visit_times), mean(null$end), quantile(null$end, 0.95)))
cat(sprintf("Episode STARTS within 36 h of a visit: observed %.2f; chance mean %.2f, 95th percentile %.2f\n",
            near_a_visit(episodes$start, visit_times), mean(null$start), quantile(null$start, 0.95)))

# 2. Turbidity step across each exchange -------------------------------------------------------
# step = median of the 6 h before arrival minus median of the 6 h after departure (+15 min).
# net_step removes what the other two depths did at the same visit (tide, storms).
by_depth <- wq_qc_final |> arrange(depth_m, date_time) |> group_split(depth_m, .keep = TRUE)
names(by_depth) <- map_chr(by_depth, ~ as.character(.x$depth_m[1]))

window_rows <- function(x, from, to) {
  times <- as.numeric(x$date_time)
  first <- findInterval(as.numeric(from) - 1, times) + 1
  last <- findInterval(as.numeric(to), times)
  if (last < first) integer(0) else first:last
}

steps <- map_dfr(seq_len(nrow(visits)), function(k) {
  map_dfr(c(1, 3, 6), function(d) {
    x <- by_depth[[as.character(d)]]
    before <- x$turb_fnu_qc[window_rows(x, visits$arr[k] - 6 * 3600, visits$arr[k])]
    after <- x$turb_fnu_qc[window_rows(x, visits$dep[k] + 900, visits$dep[k] + 6 * 3600 + 900)]
    day_before <- window_rows(x, visits$arr[k] - 24 * 3600, visits$arr[k])
    tibble(visit = k, arr = visits$arr[k], depth_m = d,
           before = if (sum(!is.na(before)) >= 12) median(before, na.rm = TRUE) else NA_real_,
           after = if (sum(!is.na(after)) >= 12) median(after, na.rm = TRUE) else NA_real_,
           flagged_before = any(!is.na(x$analysis_exclusion[day_before])))
  })
}) |>
  mutate(step = before - after, ratio = before / pmax(after, 0.1)) |>
  group_by(visit) |>
  mutate(other_step = map_dbl(seq_along(step), ~ median(step[-.x], na.rm = TRUE)),
         net_step = step - other_step) |>
  ungroup()

cat("\nTurbidity change across an exchange (before - after, FNU), by depth:\n")
steps |>
  filter(!is.na(step)) |>
  group_by(depth_m) |>
  summarise(n = n(), p05 = quantile(step, 0.05), median = median(step), p95 = quantile(step, 0.95),
            pct_over_5 = 100 * mean(step > 5), pct_over_20 = 100 * mean(step > 20)) |>
  mutate(across(-c(depth_m, n), ~ round(.x, 2))) |>
  print()

# Fouling-type steps: more than 5 FNU and 1.5x, and more than 3 FNU beyond what the other depths did
large_steps <- steps |>
  filter(!is.na(net_step), step > 5, ratio > 1.5, net_step > 3) |>
  mutate(summer = month(arr) %in% 7:9)

cat("\nLarge fouling-type steps, and whether our QC flagged the 24 h before:\n")
large_steps |>
  group_by(depth_m) |>
  summarise(n = n(), already_flagged = sum(flagged_before), not_flagged = sum(!flagged_before),
            summer_not_flagged = sum(!flagged_before & summer)) |>
  print()

cat("\nLarge steps not flagged by our QC, July-September (the analysis window):\n")
large_steps |>
  filter(!flagged_before, summer) |>
  arrange(desc(net_step)) |>
  transmute(date = as_date(arr), depth_m, before = round(before, 1), after = round(after, 1),
            net_step = round(net_step, 1)) |>
  print(n = Inf)

# Field sheet's own dirty vs clean readings for the primary sonde (surface = 1m, middle = 3m, bottom = 6m)
as_num <- function(x) suppressWarnings(as.numeric(x))
field_steps <- bind_rows(
  tibble(arr = visits$arr, depth_m = 1, field_step = as_num(visits$DirtyPrimaryTurbiditySurface) - as_num(visits$CleanPrimaryTurbiditySurface)),
  tibble(arr = visits$arr, depth_m = 3, field_step = as_num(visits$DirtyPrimaryTurbidityMiddle) - as_num(visits$CleanPrimaryTurbidityMiddle)),
  tibble(arr = visits$arr, depth_m = 6, field_step = as_num(visits$DirtyPrimaryTurbidityBottom) - as_num(visits$CleanPrimaryTurbidityBottom))
) |>
  filter(!is.na(field_step))

compare <- inner_join(filter(steps, !is.na(step)), field_steps, by = c("arr", "depth_m"))

cat("\nContinuous-record step vs field-sheet dirty minus clean:\n")
compare |>
  group_by(depth_m) |>
  summarise(n = n(), spearman = round(cor(step, field_step, method = "spearman"), 2),
            median_field = round(median(field_step), 2), median_record = round(median(step), 2)) |>
  print()

cat("\nVisits where the field sheet shows heavy fouling (dirty - clean > 10 FNU) and what the record shows:\n")
compare |>
  filter(field_step > 10) |>
  transmute(date = as_date(arr), depth_m, field_dirty_minus_clean = round(field_step, 1),
            record_step = round(step, 1), flagged_by_us = flagged_before) |>
  print(n = Inf)
