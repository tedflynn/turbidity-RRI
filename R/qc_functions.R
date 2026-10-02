# Secondary (automated/statistical) QC for sonde turbidity data.
# This runs AFTER primary field-based QC (Wagner et al. 2006 method, sensor
# logs, calibration checks) -- the goal here is to catch residual outliers,
# not to redo field QC.
#
# Each function below is one step of the pipeline and returns the input data
# frame with new columns added (nothing is dropped), so you can inspect the
# result of each step before moving to the next.

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
