# Golden Mussel clarity analysis: turning QC'd turbidity into the series that the flow-turbidity
# comparison uses. Nothing here is general QC. The general, any-station QC is in R/qc_functions.R and
# R/station_qc.R (run it with run-qc.R); the three-depth RRI special case is in R/rri_multidepth_qc.R.
# What belongs here is a decision about the clarity question: which flagged readings to leave out of
# the comparison, the summer ceiling, how readings become a daily value, and how the 2026 provisional
# data are stitched on.
#
# Requires R/qc_functions.R, R/station_qc.R and R/rri_multidepth_qc.R to be sourced first.

# Analysis-ready turbidity ------------------------------------------------------
# Builds turb_fnu_analysis for the summer clarity comparison, leaving
# turb_fnu_final untouched. Excludes (1) points flagged for review above and
# (2) any summer reading above `ceiling` FNU together with a buffer of
# `pad_readings` (24 = 6 hrs) on each side, so the ramp-up and decay of the
# same event go with it. The ceiling is a judgment call, not a sensor limit:
# summer turbidity here is normally <30 FNU (99th percentile 17-26), so
# readings above 100 are treated as transient disturbances (vessel wake,
# dredging, instrument problems) rather than ambient clarity. analysis_exclusion
# records why each point was dropped.
prepare_analysis_turbidity <- function(df, ceiling = 100, pad_readings = 24, months = 7:9) {
  df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      high_hit = !is.na(turb_fnu_final) & turb_fnu_final > ceiling & month %in% months,
      in_high_episode = month %in% months & near_flag(high_hit, pad_readings),
      analysis_exclusion = case_when(
        qaqc_flag %in% c("Flagged: Uncorroborated Excursion (review)",
                         "Flagged: Erratic Period (review)",
                         "Flagged: Depth Drift (review)",
                         "Flagged: Flat Line (review)") ~ qaqc_flag,
        in_high_episode ~ paste0("Summer episode > ", ceiling, " FNU"),
        TRUE ~ NA_character_
      ),
      turb_fnu_analysis = if_else(is.na(analysis_exclusion), turb_fnu_final, NA_real_)
    ) |>
    ungroup() |>
    select(-high_hit, -in_high_episode)
}

# Daily median of column `value`, for days with at least `min_valid` readings (48 = half of a 15-minute day)
daily_median <- function(df, value, min_valid = 48) {
  df |>
    mutate(date = as_date(date_time)) |>
    group_by(date) |>
    summarise(turb = if (sum(!is.na(.data[[value]])) >= min_valid)
                       median(.data[[value]], na.rm = TRUE) else NA_real_, .groups = "drop")
}

# Daily turbidity for every station, one row per station and day (columns date, turb, station).
#
# Single-sonde stations: the stored QC (run-qc.R), minus the flagged "(review)" readings. No summer
# ceiling is applied to them (ceiling = Inf); only the RRI record has it.
# RRI: the three-depth QC (rri_qc_final()). With use_cdec_2026 the series is the 1 m sonde, extended
#   with the provisional CDEC 1 m download after the last sonde reading, with the same 100 FNU summer
#   ceiling; otherwise it is the median across the three depths, when at least two have a valid day.
# use_cdec_2026 = FALSE also drops the provisional CDEC readings appended to the other stations' files.
# A table of the QC outcome per station is attached as attribute "qc_summary".
build_daily_turbidity <- function(stations, use_cdec_2026 = TRUE, min_valid = 48, station_levels = stations) {
  qc_summary <- list()
  single <- map(setdiff(stations, "RRI"), function(s) {
    qc <- load_station_qc(s)
    if (!use_cdec_2026) qc <- filter(qc, source == "primary")
    qc_summary[[s]] <<- tibble(station = s, pct_valid = round(100 * mean(qc$qaqc_flag == "Valid", na.rm = TRUE), 1),
                               first = as_date(min(qc$date_time[!is.na(qc$turb_fnu)])),
                               readings_with_turbidity = sum(!is.na(qc$turb_fnu)))
    qc |>
      prepare_analysis_turbidity(ceiling = Inf) |>
      daily_median("turb_fnu_analysis", min_valid) |>
      mutate(station = s)
  })

  rri_qc <- rri_qc_final()
  if (use_cdec_2026) {
    last_sonde <- max(rri_qc$date_time)
    rri_2026 <- load_station_qc("RRI_cdec_2026") |>
      filter(date_time > last_sonde) |>
      prepare_analysis_turbidity()
    rri <- bind_rows(filter(rri_qc, depth_m == 1), rri_2026) |>
      daily_median("turb_fnu_analysis", min_valid) |>
      mutate(station = "RRI")
  } else {
    rri <- rri_qc |>
      group_by(depth_m) |>
      group_modify(~ daily_median(.x, "turb_fnu_analysis", min_valid)) |>
      group_by(date) |>
      summarise(turb = if (sum(!is.na(turb)) >= 2) median(turb, na.rm = TRUE) else NA_real_, .groups = "drop") |>
      mutate(station = "RRI")
  }

  daily <- bind_rows(bind_rows(single), rri) |>
    filter(!is.na(turb), turb > 0) |>
    mutate(station = factor(station, levels = station_levels))
  attr(daily, "qc_summary") <- bind_rows(qc_summary)
  daily
}
