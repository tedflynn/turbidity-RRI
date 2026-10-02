# Secondary (automated/statistical) QC for sonde turbidity data.
# This runs AFTER primary field-based QC (Wagner et al. 2006 method, sensor
# logs, calibration checks) -- the goal here is to catch residual outliers,
# not to redo field QC.
#
# Each function below is one step of the pipeline and returns the input data
# frame with new columns added (nothing is dropped), so you can inspect the
# result of each step before moving to the next.

# Helper: TRUE for every row within +/- pad rows of a TRUE in `flag`. Uses a
# cumulative-sum window difference so it stays fast on ~600k-row series.
# `flag` must be logical with no NAs and ordered by time within one depth.
near_flag <- function(flag, pad) {
  n <- length(flag)
  cs <- cumsum(flag)
  idx <- seq_len(n)
  lo <- idx - pad - 1L
  (cs[pmin(idx + pad, n)] - if_else(lo < 1L, 0L, cs[pmax(lo, 1L)])) > 0
}

# Step 1: Gross/physical range check -------------------------------------
# Flags values that are physically impossible for the sensor, not values
# that are merely unusual. max_val defaults to the YSI EXO2 turbidity smart
# sensor's rated range (0-4000 FNU) -- confirm this against your sensor's
# calibration sheet, since a non-default sensor range would change it.
check_gross_range <- function(df, min_val = 0, max_val = 4000) {
  df |>
    mutate(
      gross_range_flag = case_when(
        is.na(turb_fnu) ~ NA_character_,
        turb_fnu < min_val | turb_fnu > max_val ~ "Fail: Gross Range",
        TRUE ~ "Pass"
      ),
      turb_fnu_qc = if_else(gross_range_flag == "Fail: Gross Range", NA_real_, turb_fnu)
    )
}

# Step 2: Robust local outlier detection ---------------------------------
# Turbidity is heavy-tailed (mostly low, occasional large real spikes), so
# this uses a rolling median/MAD (robust to outliers) rather than
# mean/SD, scored as a modified z-score (Iglewicz & Hoaglin, 1993):
#   modified_z = 0.6745 * (x - rolling_median) / rolling_MAD
# noise_floor guards against flagging trivial deviations during quiet
# periods when the local MAD is near zero.
flag_robust_outliers <- function(df, window = 12, z_thresh = 3.5, noise_floor = 5) {
  df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      roll_median = slide_dbl(turb_fnu_qc, median, na.rm = TRUE, .before = window, .after = window),
      roll_mad = slide_dbl(turb_fnu_qc, mad, na.rm = TRUE, .before = window, .after = window),
      abs_dev = abs(turb_fnu_qc - roll_median),
      modified_z = 0.6745 * (turb_fnu_qc - roll_median) / roll_mad,
      spike_flag = case_when(
        is.na(turb_fnu_qc) ~ NA_character_,
        abs_dev > noise_floor & abs(modified_z) > z_thresh ~ "Fail: Outlier",
        TRUE ~ "Pass"
      )
    ) |>
    ungroup()
}

# Step 2b: Erratic-period detection -----------------------------------------
# Step 2 compares each point to a +/-3 hr rolling baseline, so it goes blind
# when a disturbance lasts longer than the window: the baseline itself fills
# with disturbed data. This step instead measures how erratic the series is,
# using the median absolute 15-min change over a centered 24 hr window.
# Quiet data moves ~0.4 FNU per step; the multi-week erratic periods found in
# 3m (Sept-Oct 2008) and 6m (Sept 2014) move ~90 FNU per step. The median
# (not mean) keeps a few isolated spikes from tripping the flag.
flag_erratic_periods <- function(df, window = 48, jump_thresh = 5) {
  df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      abs_jump = abs(turb_fnu_qc - lag(turb_fnu_qc)),
      med_jump = slide_dbl(abs_jump, median, na.rm = TRUE, .before = window, .after = window),
      erratic_flag = case_when(
        is.na(turb_fnu_qc) ~ NA_character_,
        !is.na(med_jump) & med_jump > jump_thresh ~ "Fail: Erratic Period",
        TRUE ~ "Pass"
      )
    ) |>
    ungroup()
}

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

# Step 3: Duration-aware classification ----------------------------------
# A single anomalous reading that immediately returns to baseline is a much
# stronger glitch candidate than a run of several consecutive flagged
# readings (which looks more like a real, sustained event). max_spike_len
# is in units of 15-min readings (2 = 30 minutes).
classify_run_duration <- function(df, max_spike_len = 2) {
  df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      is_outlier = !is.na(spike_flag) & spike_flag == "Fail: Outlier",
      run_id = consecutive_id(is_outlier)
    ) |>
    group_by(depth_m, run_id) |>
    mutate(run_length = if_else(is_outlier, n(), NA_integer_)) |>
    ungroup() |>
    mutate(
      duration_class = case_when(
        !is_outlier ~ NA_character_,
        run_length <= max_spike_len ~ "Isolated Spike",
        TRUE ~ "Sustained Excursion"
      )
    )
}

# Step 4: Cross-parameter corroboration ----------------------------------
# For "Sustained Excursion" runs only: a real event (storm, resuspension)
# should show up in other sonde channels too (e.g. conductivity drops with
# freshwater dilution); a turbidity-only excursion with no corroborating
# signal elsewhere is a stronger fouling/sensor-artifact candidate.
# Reuses the same robust-z approach from step 2 on spcond_uscm and do_mgl.
corroborate_sustained_excursions <- function(df, window = 12, corrob_z_thresh = 2) {
  df <- df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      cond_roll_median = slide_dbl(spcond_uscm, median, na.rm = TRUE, .before = window, .after = window),
      cond_roll_mad = slide_dbl(spcond_uscm, mad, na.rm = TRUE, .before = window, .after = window),
      cond_z = 0.6745 * (spcond_uscm - cond_roll_median) / cond_roll_mad,
      do_roll_median = slide_dbl(do_mgl, median, na.rm = TRUE, .before = window, .after = window),
      do_roll_mad = slide_dbl(do_mgl, mad, na.rm = TRUE, .before = window, .after = window),
      do_z = 0.6745 * (do_mgl - do_roll_median) / do_roll_mad,
      row_corroborated = (!is.na(cond_z) & abs(cond_z) > corrob_z_thresh) |
        (!is.na(do_z) & abs(do_z) > corrob_z_thresh)
    ) |>
    ungroup()

  run_corrob <- df |>
    filter(duration_class == "Sustained Excursion") |>
    group_by(depth_m, run_id) |>
    summarise(run_corroborated = any(row_corroborated), .groups = "drop")

  df |> left_join(run_corrob, by = c("depth_m", "run_id"))
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

# Step 6: Analysis-ready turbidity ------------------------------------------
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
                         "Flagged: Depth Drift (review)") ~ qaqc_flag,
        in_high_episode ~ paste0("Summer episode > ", ceiling, " FNU"),
        TRUE ~ NA_character_
      ),
      turb_fnu_analysis = if_else(is.na(analysis_exclusion), turb_fnu_final, NA_real_)
    ) |>
    ungroup() |>
    select(-high_hit, -in_high_episode)
}

# Step 5: Final QC decision -----------------------------------------------
# Gross-range failures and isolated spikes are near-certain glitches, so
# those are removed (set to NA). A sustained excursion counts as
# corroborated if EITHER another channel (step 4) OR another depth
# (step 4b) shows a coincident signal; only excursions confined to a single
# depth's turbidity reading alone are flagged for manual review rather than
# removed, since deleting a real event (storm, vessel-wake resuspension)
# would bias downstream analysis.
finalize_qc_flags <- function(df) {
  df |>
    mutate(
      any_corroborated = run_corroborated | depth_corroborated,
      qaqc_flag = case_when(
        gross_range_flag == "Fail: Gross Range" ~ "Removed: Gross Range",
        erratic_flag == "Fail: Erratic Period" ~ "Flagged: Erratic Period (review)",
        drift_flag == "Fail: Depth Drift" ~ "Flagged: Depth Drift (review)",
        duration_class == "Isolated Spike" ~ "Removed: Isolated Spike",
        duration_class == "Sustained Excursion" & !any_corroborated ~ "Flagged: Uncorroborated Excursion (review)",
        duration_class == "Sustained Excursion" & any_corroborated ~ "Kept: Corroborated Event",
        TRUE ~ "Valid"
      ),
      turb_fnu_final = if_else(
        qaqc_flag %in% c("Removed: Gross Range", "Removed: Isolated Spike"),
        NA_real_,
        turb_fnu_qc
      )
    )
}
