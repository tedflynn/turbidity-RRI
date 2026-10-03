# Golden Mussel analysis, special case: the RRI sonde string (1, 3 and 6 m in one housing).
#
# Almost no other station has sondes at several depths, so the cross-depth steps do not belong in
# the general QC (R/qc_functions.R, R/station_qc.R). They are kept here, with the function that
# runs the full three-depth sequence for RRI and the loader for the RRI files.
#
# Cross-depth steps (need all of depth_m = 1, 3, 6 in one data frame):
#   flag_depth_drift(), flag_param_depth_drift(), corroborate_cross_depth()
# Entry points:
#   load_rri_raw()       three sonde files -> one data frame, clock fixed to PST, water-year type joined
#   qc_rri_three_depth() the whole QC sequence (general steps plus the cross-depth ones)
#   rri_qc_final()       the QC'd data frame, from the cache in data/qc/ when present
#
# Requires R/qc_functions.R (near_flag() and the general steps) to be sourced first.

# Step 2c: Cross-depth drift detection --------------------------------------
# Biofouling shows up as a slow one-sensor climb (e.g. 1m rising from ~9 to
# 73 FNU over Sept 27-Oct 2 2018, then resetting at a service visit) that
# steps 2 and 2b can't see: no spikes, and the local baseline just follows
# the ramp. All three sondes share one housing, so each depth's 24 hr median
# level is compared with the mean of the other depths' levels. A depth is
# flagged when it stays more than `diff_thresh` FNU AND `ratio_thresh` times
# above them for at least `min_readings` in a row (96 = 24 hrs), with a buffer
# of `pad_readings` (96 = 24 hrs) either side to catch the ramp-up. Upward
# only: fouling reads high. level_24h/other_level_24h are kept for inspection.
# Defaults (8 FNU, 3x) were chosen so steady season-long offsets between depths
# (e.g. 1m ~10 vs ~5 FNU in summer 2009-2010, 6m ~18 FNU in 2011, with
# chlorophyll unchanged) count as normal depth differences, not drift; at the
# looser 5 FNU / 2x they would exclude 30-47% of those summers.
flag_depth_drift <- function(df, window = 48, diff_thresh = 8, ratio_thresh = 3,
                             min_readings = 96, pad_readings = 96) {
  lev <- df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(level_24h = slide_dbl(turb_fnu_qc, median, na.rm = TRUE, .before = window, .after = window)) |>
    ungroup()

  wide <- lev |>
    select(date_time, depth_m, level_24h) |>
    pivot_wider(names_from = depth_m, values_from = level_24h, names_prefix = "lvl_")
  lvl_cols <- grep("^lvl_", names(wide), value = TRUE)

  other <- map_dfr(sort(unique(lev$depth_m)), function(d) {
    m <- rowMeans(wide[setdiff(lvl_cols, paste0("lvl_", d))], na.rm = TRUE)
    tibble(date_time = wide$date_time, depth_m = d, other_level_24h = if_else(is.nan(m), NA_real_, m))
  })

  lev |>
    left_join(other, by = c("date_time", "depth_m")) |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      candidate = !is.na(level_24h) & !is.na(other_level_24h) &
        (level_24h - other_level_24h) > diff_thresh &
        level_24h / other_level_24h > ratio_thresh,
      run_id = consecutive_id(candidate)
    ) |>
    group_by(depth_m, run_id) |>
    mutate(sustained = candidate & n() >= min_readings) |>
    group_by(depth_m) |>
    mutate(
      in_drift = near_flag(sustained, pad_readings),
      drift_flag = case_when(
        is.na(turb_fnu_qc) ~ NA_character_,
        in_drift ~ "Fail: Depth Drift",
        TRUE ~ "Pass"
      )
    ) |>
    ungroup() |>
    select(-candidate, -run_id, -sustained, -in_drift)
}

# Step 2e: Cross-depth drift check for another parameter (used for pH) ---------
# Same idea as step 2c, but built for a bounded, roughly symmetric variable
# like pH. Each depth's 24 hr median level is compared with the MEDIAN of all
# three depths' levels (not the mean of the other two): when one sensor fails,
# the median stays with the two healthy sensors, so the blame lands on the odd
# one out instead of spreading to its neighbors. Needs all three depths present.
# Two-sided. A depth is flagged when its level differs from the reference by
# more than `diff_thresh` (the parameter's own units) AND, if `ratio_thresh` > 1,
# by more than that multiplicative factor, for `min_readings` in a row (192 =
# 48 hrs), plus `pad_readings` (48 = 12 hrs) either side. Thresholds are set
# beyond the 99th percentile of normal depth-to-depth differences (some real
# stratification exists: 1m runs ~0.5 C and ~0.7 mg/L DO above 6m at the 95th
# percentile). Settings used for each parameter are in turbidity-depth-RRI.R.
# Adds <prefix>_level_24h, <prefix>_ref_24h, <prefix>_drift_flag.
flag_param_depth_drift <- function(df, param = "ph_qc", prefix = "ph", diff_thresh = 0.5,
                                   ratio_thresh = 1, window = 48, min_readings = 192,
                                   pad_readings = 48) {
  lev <- df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(.lvl = slide_dbl(.data[[param]], median, na.rm = TRUE, .before = window, .after = window)) |>
    ungroup()

  wide <- lev |>
    select(date_time, depth_m, .lvl) |>
    pivot_wider(names_from = depth_m, values_from = .lvl, names_prefix = "d")
  a <- wide[[2]]; b <- wide[[3]]; c <- wide[[4]]
  ref <- pmax(pmin(a, b), pmin(pmax(a, b), c))   # median of three; NA if any is NA
  ref_tbl <- tibble(date_time = wide$date_time, .ref = ref)

  out <- lev |>
    left_join(ref_tbl, by = "date_time") |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      candidate = !is.na(.lvl) & !is.na(.ref) & abs(.lvl - .ref) > diff_thresh &
        pmax(.lvl, .ref) / pmax(pmin(.lvl, .ref), 1e-6) > ratio_thresh,
      run_id = consecutive_id(candidate)
    ) |>
    group_by(depth_m, run_id) |>
    mutate(sustained = candidate & n() >= min_readings) |>
    group_by(depth_m) |>
    mutate(
      in_drift = near_flag(sustained, pad_readings),
      .flag = case_when(is.na(.data[[param]]) ~ NA_character_, in_drift ~ "Fail: Depth Drift", TRUE ~ "Pass")
    ) |>
    ungroup() |>
    select(-candidate, -run_id, -sustained, -in_drift)

  names(out)[names(out) == ".lvl"] <- paste0(prefix, "_level_24h")
  names(out)[names(out) == ".ref"] <- paste0(prefix, "_ref_24h")
  names(out)[names(out) == ".flag"] <- paste0(prefix, "_drift_flag")
  out
}

# Step 4b: Cross-depth corroboration ---------------------------------------
# All three sondes share one housing at this site, so a real water-column
# event (vessel-wake resuspension, a storm) should register in turbidity at
# the other depths too, even when it doesn't move conductivity or DO. A
# disturbance confined to one depth's sensor is a stronger single-sensor-
# artifact signal. Reuses the per-depth modified_z from step 2 -- no new
# rolling-window computation needed.
# pad_readings extends the comparison window a bit on each side (in units of
# 15-min readings) to allow for a short lag between depths.
corroborate_cross_depth <- function(df, z_thresh = 3.5, pad_readings = 1) {
  z_wide <- df |>
    select(date_time, depth_m, modified_z) |>
    arrange(date_time) |>
    pivot_wider(names_from = depth_m, values_from = modified_z, names_prefix = "z_depth_")

  # Timestamps are a confirmed regular 15-min grid, so a run's row range in
  # z_wide can be found by index arithmetic instead of a per-run table scan.
  t0 <- min(z_wide$date_time)
  idx_of <- function(t) as.integer(round(as.numeric(difftime(t, t0, units = "mins")) / 15)) + 1L

  runs <- df |>
    filter(duration_class == "Sustained Excursion") |>
    group_by(depth_m, run_id) |>
    summarise(start = min(date_time), end = max(date_time), .groups = "drop")

  depth_cols <- paste0("z_depth_", sort(unique(df$depth_m)))

  runs$depth_corroborated <- purrr::pmap_lgl(runs, function(depth_m, run_id, start, end) {
    i0 <- max(1L, idx_of(start) - pad_readings)
    i1 <- min(nrow(z_wide), idx_of(end) + pad_readings)
    other_cols <- setdiff(depth_cols, paste0("z_depth_", depth_m))
    vals <- unlist(z_wide[i0:i1, other_cols])
    any(abs(vals) > z_thresh, na.rm = TRUE)
  })

  df |> left_join(runs |> select(depth_m, run_id, depth_corroborated), by = c("depth_m", "run_id"))
}

# Loading and QC ------------------------------------------------------------------------------------------

# Water-year types from CDEC (data/WY-type.csv), San Joaquin version. Returns water_year and wy_type.
read_sj_wy_types <- function(path = here("data", "WY-type.csv")) {
  read_csv(path, show_col_types = FALSE) |>
    select("Year", "Yr-type-SJR") |>
    rename("wy_type" = "Yr-type-SJR") |>
    rename("water_year" = "Year")
}

# The three RRI sonde files (data provided QC'd by Andrew Tran in C-EMP), one data frame with depth_m,
# calendar columns and the San Joaquin water-year type.
load_rri_raw <- function() {
  df_1m <- read_csv(here("data", "RRI_data_1m.csv"), show_col_types = FALSE) |> mutate(depth_m = 1)
  df_3m <- read_csv(here("data", "RRI_data_3m.csv"), show_col_types = FALSE) |> mutate(depth_m = 3)
  df_6m <- read_csv(here("data", "RRI_data_6m.csv"), show_col_types = FALSE) |> mutate(depth_m = 6)

  bind_rows(df_1m, df_3m, df_6m) |>
    mutate(
      # The raw strings carry a "-08:00" suffix, but the clock values are really UTC:
      # daily temperature/DO/pH cycles, field-sheet exchange times and EMP grab samples
      # all line up only after subtracting 8 hours. read_csv() parses the string clock
      # as if it were PST, so shift back 8 h, then express in fixed PST (no DST).
      date_time = with_tz(date_time - hours(8), tzone = "Etc/GMT+8"),
      year = year(date_time),
      month = month(date_time),
      day = day(date_time),
      water_year = if_else(month >= 10, year + 1, year)   # water year begins October 1st
    ) |>
    left_join(read_sj_wy_types(), by = "water_year") |>
    # Factor levels in wetness order so graphs do not plot alphabetically
    mutate(wy_type = factor(wy_type, levels = c("W", "AN", "BN", "D", "C")))
}

# Secondary QC for the three-depth RRI record. Field-based primary QC (Wagner et al. 2006 method) was
# applied upstream; these steps catch residual statistical outliers. Every step adds columns without
# dropping rows. Returns the data frame with qaqc_flag, turb_fnu_final and the analysis columns
# (turb_fnu_analysis, analysis_exclusion; summer ceiling per prepare_analysis_turbidity()).
qc_rri_three_depth <- function(raw, ceiling = 100) {
  # Basic QC: remove erroneous negative sonde turbidity readings
  wq <- raw |> filter(turb_fnu >= 0 | is.na(turb_fnu))

  wq <- wq |>
    check_gross_range(min_val = 0, max_val = 4000) |>                  # 1: physical range
    check_param_ranges() |>                                              # 1b: other parameters
    flag_robust_outliers(window = 12, z_thresh = 3.5, noise_floor = 5) |> # 2: rolling median / MAD
    # 2b: sustained disturbances longer than step 2's window (flagged for review, not removed)
    flag_erratic_periods(window = 48, jump_thresh = 5) |>
    # 2c: one sensor climbing away from the other two, e.g. biofouling that resets at a service visit
    flag_depth_drift(diff_thresh = 8, ratio_thresh = 3) |>
    # 2d: 12+ identical readings in a row = stuck sensor
    flag_flat_lines(min_run = 12, pad_readings = 24)

  # 2e: each of the other parameters drifting away from the other two depths (24 hr median vs the median
  # of all three depths, sustained 48+ hrs). Adds <prefix>_drift_flag columns for any analysis of those
  # parameters; does not change turbidity decisions. Thresholds sit beyond the 99th percentile of normal
  # depth-to-depth differences (see flag_param_depth_drift()).
  drift_settings <- list(
    ph   = list(param = "ph_qc",          diff_thresh = 0.5, ratio_thresh = 1),
    temp = list(param = "watertemp_c_qc", diff_thresh = 1.0, ratio_thresh = 1),
    cond = list(param = "spcond_uscm_qc", diff_thresh = 100, ratio_thresh = 1.25),
    do   = list(param = "do_mgl_qc",      diff_thresh = 1.5, ratio_thresh = 1),
    chl  = list(param = "chlor_ugl_qc",   diff_thresh = 5,   ratio_thresh = 2)
  )
  for (p in names(drift_settings)) {
    s <- drift_settings[[p]]
    wq <- flag_param_depth_drift(wq, param = s$param, prefix = p,
                                 diff_thresh = s$diff_thresh, ratio_thresh = s$ratio_thresh)
  }

  wq |>
    classify_run_duration(max_spike_len = 2) |>                                    # 3: spike vs sustained
    corroborate_sustained_excursions(window = 12, corrob_z_thresh = 2) |>          # 4: other channels
    # 4b: a real event in one housing should show up at more than one depth
    corroborate_cross_depth(z_thresh = 3.5, pad_readings = 1) |>
    finalize_qc_flags() |>                                                         # 5: final decision
    prepare_analysis_turbidity(ceiling = ceiling, pad_readings = 24)               # 6: analysis column
}

# The QC'd RRI record. The QC takes a few minutes, so the result is cached in data/qc/ (gitignored);
# refresh = TRUE reruns it, which is needed after the sonde files or any QC setting change.
rri_qc_final <- function(refresh = FALSE) {
  cache <- here("data", "qc", "RRI_three_depth_qc.rds")
  if (!refresh && file.exists(cache)) return(readRDS(cache))
  out <- qc_rri_three_depth(load_rri_raw())
  dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, cache)
  out
}
