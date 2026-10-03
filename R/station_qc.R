# Secondary QC for a single-sonde continuous water quality file (columns: station, date, time,
# dissolvedoxygen, fluorescence, ph, spc, turbidity, watertemperature), using the functions in
# R/qc_functions.R. Source that file first (it also defines param_limits).
#
# Runs: turbidity range, range screen for the other parameters, robust outliers, erratic
# periods, flat lines, spike vs sustained classification and cross-parameter corroboration
# (conductivity and DO). Cannot run the cross-depth steps (one sonde), and does not apply the
# RRI-specific 100 FNU summer ceiling. File clocks must already be PST (checked per file).
#
# `limits` are the parameter range limits (see param_limits); override per station, e.g. a
# brackish station needs a much higher conductivity ceiling.
qc_single_sonde <- function(path, limits = param_limits) {
  raw <- read_csv(path, show_col_types = FALSE) |>
    transmute(station, depth_m = 1,
              date_time = ymd_hms(paste(date, time), tz = "Etc/GMT+8"),
              year = year(date_time), month = month(date_time),
              do_mgl = dissolvedoxygen, chlor_ugl = fluorescence, ph,
              spcond_uscm = spc, turb_fnu = turbidity, watertemp_c = watertemperature)

  qc <- raw |>
    check_gross_range() |>
    check_param_ranges(limits) |>
    flag_robust_outliers() |>
    flag_erratic_periods() |>
    flag_flat_lines() |>
    classify_run_duration() |>
    corroborate_sustained_excursions()

  # finalize_qc_flags() expects the depth-based columns; with one sonde they carry no signal
  qc$drift_flag <- NA_character_
  qc$depth_corroborated <- FALSE
  qc <- finalize_qc_flags(qc)

  qc |>
    mutate(analysis_exclusion = if_else(
      qaqc_flag %in% c("Flagged: Erratic Period (review)", "Flagged: Flat Line (review)",
                       "Flagged: Uncorroborated Excursion (review)"), qaqc_flag, NA_character_),
      turb_fnu_analysis = if_else(is.na(analysis_exclusion), turb_fnu_final, NA_real_))
}
