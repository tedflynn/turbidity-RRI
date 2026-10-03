# Secondary QC for single-sonde continuous water quality files: a reader for the file layouts
# in use, and qc_single_sonde(), which runs the shared QC steps.
#
# Reads a station file into the QC schema. Two layouts are recognised:
#   - station, date, time, dissolvedoxygen, fluorescence, ph, spc, turbidity, watertemperature
#   - DWR layout: Station_Code, Date_Time (month/day/year hh:mm:ss), then any of Water_Temperature,
#     Specific_Conductance, DO, Turbidity, pH, Chlorophyll
#   - CDEC CSVDataServlet download of one sensor (STATION_ID, SENSOR_NUMBER, DATE TIME, VALUE, ...):
#     turbidity only (sensors 27 and 221), real-time and not QC'd. CDEC writes "---" for missing.
# Parameters a station does not measure become NA columns (the QC steps skip them). Both layouts
# are Pacific Standard Time; this was checked per file from the daily temperature and DO cycles.
read_station_file <- function(path) {
  d <- read_csv(path, show_col_types = FALSE, col_types = cols(.default = col_character()))
  num <- function(nm) if (nm %in% names(d)) suppressWarnings(as.numeric(d[[nm]])) else rep(NA_real_, nrow(d))
  if ("SENSOR_NUMBER" %in% names(d)) {
    out <- tibble(station = d$STATION_ID, depth_m = 1,
                  date_time = ymd_hm(d$`DATE TIME`, tz = "Etc/GMT+8", quiet = TRUE),
                  do_mgl = NA_real_, chlor_ugl = NA_real_, ph = NA_real_, spcond_uscm = NA_real_,
                  turb_fnu = num("VALUE"), watertemp_c = NA_real_)
  } else if ("Date_Time" %in% names(d)) {
    out <- tibble(station = d$Station_Code, depth_m = 1,
                  date_time = mdy_hms(d$Date_Time, tz = "Etc/GMT+8", quiet = TRUE),
                  do_mgl = num("DO"), chlor_ugl = num("Chlorophyll"), ph = num("pH"),
                  spcond_uscm = num("Specific_Conductance"), turb_fnu = num("Turbidity"),
                  watertemp_c = num("Water_Temperature"))
  } else {
    out <- tibble(station = d$station, depth_m = 1,
                  date_time = ymd_hms(paste(d$date, d$time), tz = "Etc/GMT+8", quiet = TRUE),
                  do_mgl = num("dissolvedoxygen"), chlor_ugl = num("fluorescence"), ph = num("ph"),
                  spcond_uscm = num("spc"), turb_fnu = num("turbidity"),
                  watertemp_c = num("watertemperature"))
  }
  out |>
    filter(!is.na(date_time)) |>
    mutate(year = year(date_time), month = month(date_time), .after = date_time)
}

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
#
# `extra_paths` are further files (e.g. a CDEC download of the latest months) whose readings after
# the last time in `path` are appended before the QC. qc_station_data() takes the data directly.
qc_single_sonde <- function(path, limits = param_limits, extra_paths = character()) {
  raw <- read_station_file(path)
  for (p in extra_paths) {
    raw <- bind_rows(raw, filter(read_station_file(p), date_time > max(raw$date_time)))
  }
  qc_station_data(raw, limits)
}

qc_station_data <- function(raw, limits = param_limits) {
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
