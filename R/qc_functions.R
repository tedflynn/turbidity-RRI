# Secondary (automated/statistical) QC for sonde turbidity data.
# This runs AFTER primary field-based QC (Wagner et al. 2006 method, sensor
# logs, calibration checks) -- the goal here is to catch residual outliers,
# not to redo field QC.
#
# Each function below is one step of the pipeline and returns the input data
# frame with new columns added (nothing is dropped), so you can inspect the
# result of each step before moving to the next.
#
# These steps work on one sonde at a time (a data frame with station, depth_m,
# date_time and the parameter columns; depth_m is just a grouping key and is 1
# for a station with one sonde). Project-specific decisions are kept out of
# here: the three-depth RRI steps are in R/rri_multidepth_qc.R and the choice
# of which flagged readings to leave out of the clarity comparison is in
# R/mussel_analysis.R. R/station_qc.R chains these steps for any station file.

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

# Step 1b: Range screen for the other sonde parameters ---------------------
# Same idea as step 1, applied to temperature, pH, conductivity, DO and
# chlorophyll. Limits are "impossible for this river or this sensor", not
# "unusual": they are deliberately wide. Out-of-range values are set to NA in
# a new <param>_qc column (originals untouched) with a range_flag_<param>
# column recording why. Defaults, with what they catch in this record:
#   watertemp_c 5-35 C      9999 codes; 42-100 C readings; 0.1-0.4 C readings (3m Sept 2009, 1m Jun 2022)
#   ph          5.5-10      0.9-5.5 and 10-277 (2009 3m failure, Mar 2008, Feb 2021)
#   spcond_uscm 20-5000     1-14 uS/cm and 5,597-134,716 uS/cm
#   do_mgl      0-25 mg/L   41-1,012 mg/L (Feb 2021)
#   chlor_ugl   0-400 ug/L  EXO total-algae sensor range; ~7,400 readings above it,
#                           mostly long blocks (3m Sep-Oct 2008, 1m May-Jun 2009, 6m Mar-Apr 2014)
# Confirm the sensor-range ones (chlorophyll, temperature) against your sensor
# spec sheets and the site-based ones (pH, conductivity) against your site knowledge.
param_limits <- list(
  watertemp_c = c(5, 35),
  ph          = c(5.5, 10),
  spcond_uscm = c(20, 5000),
  do_mgl      = c(0, 25),
  chlor_ugl   = c(0, 400)
)

check_param_ranges <- function(df, limits = param_limits) {
  for (p in names(limits)) {
    x <- df[[p]]
    bad <- !is.na(x) & (x < limits[[p]][1] | x > limits[[p]][2])
    df[[paste0("range_flag_", p)]] <- case_when(is.na(x) ~ NA_character_, bad ~ "Fail: Range", TRUE ~ "Pass")
    df[[paste0(p, "_qc")]] <- if_else(bad, NA_real_, x)
  }
  df
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

# Step 2d: Flat-line (stuck sensor) detection ----------------------------------
# A sensor that stops updating repeats one value. At 0.1 FNU resolution, short
# repeats are normal in clear water (runs of 6-9 identical readings at 1-6 FNU
# are common), so only runs of `min_run` or more (12 = 3 hrs) are flagged,
# plus `pad_readings` (24 = 6 hrs) either side so neighboring stuck stretches
# (e.g. 1m, Jan 3-11 2016, 27 hrs at exactly 2.8 FNU; 1427 FNU for 5+ hrs on
# Dec 15 2015) are caught together. NA breaks a run.
flag_flat_lines <- function(df, min_run = 12, pad_readings = 24) {
  df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      run_id = consecutive_id(turb_fnu_qc),
      is_flat = !is.na(turb_fnu_qc) & ave(seq_along(run_id), run_id, FUN = length) >= min_run,
      flat_flag = case_when(
        is.na(turb_fnu_qc) ~ NA_character_,
        near_flag(is_flat, pad_readings) ~ "Fail: Flat Line",
        TRUE ~ "Pass"
      )
    ) |>
    ungroup() |>
    select(-run_id, -is_flat)
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
# Reuses the same robust-z approach from step 2 on spcond_uscm_qc and do_mgl_qc (range-screened in step 1b).
corroborate_sustained_excursions <- function(df, window = 12, corrob_z_thresh = 2) {
  df <- df |>
    group_by(depth_m) |>
    arrange(date_time, .by_group = TRUE) |>
    mutate(
      cond_roll_median = slide_dbl(spcond_uscm_qc, median, na.rm = TRUE, .before = window, .after = window),
      cond_roll_mad = slide_dbl(spcond_uscm_qc, mad, na.rm = TRUE, .before = window, .after = window),
      cond_z = 0.6745 * (spcond_uscm_qc - cond_roll_median) / cond_roll_mad,
      do_roll_median = slide_dbl(do_mgl_qc, median, na.rm = TRUE, .before = window, .after = window),
      do_roll_mad = slide_dbl(do_mgl_qc, mad, na.rm = TRUE, .before = window, .after = window),
      do_z = 0.6745 * (do_mgl_qc - do_roll_median) / do_roll_mad,
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

# Step 5: Final QC decision -----------------------------------------------
# Gross-range failures and isolated spikes are near-certain glitches, so
# those are removed (set to NA). A sustained excursion counts as
# corroborated if another channel (step 4) shows a coincident signal; only
# excursions confined to the turbidity reading alone are flagged for manual
# review rather than removed, since deleting a real event (storm, vessel-wake
# resuspension) would bias downstream analysis.
# Stations with sondes at several depths (RRI) add two inputs from
# R/rri_multidepth_qc.R: drift_flag (one depth drifting away from the others)
# and depth_corroborated (the excursion also shows at the other depths, which
# counts as corroboration). Without those columns they are ignored.
finalize_qc_flags <- function(df) {
  depth_drift <- if ("drift_flag" %in% names(df)) df$drift_flag else NA_character_
  depth_corr <- if ("depth_corroborated" %in% names(df)) df$depth_corroborated else FALSE
  df |>
    mutate(
      any_corroborated = run_corroborated | .env$depth_corr,
      qaqc_flag = case_when(
        gross_range_flag == "Fail: Gross Range" ~ "Removed: Gross Range",
        erratic_flag == "Fail: Erratic Period" ~ "Flagged: Erratic Period (review)",
        .env$depth_drift == "Fail: Depth Drift" ~ "Flagged: Depth Drift (review)",
        flat_flag == "Fail: Flat Line" ~ "Flagged: Flat Line (review)",
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
