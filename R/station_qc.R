# General secondary QC for continuous water quality stations with one sonde: a reader for the
# file layouts in use, qc_station_data() / qc_single_sonde(), which chain the steps in
# R/qc_functions.R, and the helpers that run QC for a list of stations and store the results
# (run_station_qc(), load_station_qc(); driver script run-qc.R). Nothing here depends on the
# Golden Mussel question: that analysis reads the stored QC output (see R/mussel_analysis.R).
#
# Reads a station file into the QC schema. Three layouts are recognised:
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

# Secondary QC for a single-sonde continuous water quality file (any of the layouts above), using the
# functions in R/qc_functions.R. Source that file first (it also defines param_limits).
#
# Runs: turbidity range, range screen for the other parameters, robust outliers, erratic
# periods, flat lines, spike vs sustained classification and cross-parameter corroboration
# (conductivity and DO). Cannot run the cross-depth steps (one sonde). File clocks must
# already be PST (checked per file).
#
# `limits` are the parameter range limits (see param_limits); override per station, e.g. a
# brackish station needs a much higher conductivity ceiling.
#
# `extra_paths` are further files (e.g. a CDEC download of the latest months) whose readings after
# the last time in `path` are appended before the QC. The `source` column records where each
# reading came from ("primary" for `path`, "provisional" for the extra files). qc_station_data()
# takes the data directly.
qc_single_sonde <- function(path, limits = param_limits, extra_paths = character()) {
  raw <- read_station_file(path) |> mutate(source = "primary")
  for (p in extra_paths) {
    raw <- bind_rows(raw, filter(read_station_file(p), date_time > max(raw$date_time)) |> mutate(source = "provisional"))
  }
  qc_station_data(raw, limits)
}

# Adds the QC columns (qaqc_flag, turb_fnu_final, and the flag columns of each step) to a data frame
# in the schema of read_station_file(). Readings are flagged, not deleted, except for gross-range
# failures and isolated spikes, which are set to NA in turb_fnu_final.
qc_station_data <- function(raw, limits = param_limits) {
  raw |>
    check_gross_range() |>
    check_param_ranges(limits) |>
    flag_robust_outliers() |>
    flag_erratic_periods() |>
    flag_flat_lines() |>
    classify_run_duration() |>
    corroborate_sustained_excursions() |>
    finalize_qc_flags()
}

# Station QC jobs and stored results ------------------------------------------------------------------

# Conductivity ceiling for the brackish, tidal stations
tidal_cond <- modifyList(param_limits, list(spcond_uscm = c(20, 30000)))

# One entry per QC'd series: `file` is the station file in data/, `limits` its parameter range limits,
# `extra` further files whose later readings are appended (the provisional CDEC downloads of 2026).
# RRI has its own three-depth QC (R/rri_multidepth_qc.R); only its provisional CDEC 1 m download is here.
cdec_2026 <- function(s) paste0("CDEC_", s, "_turb_2026.csv")
qc_jobs <- c(
  setNames(lapply(c("SJR", "MSD", "PPT", "TRN", "MOK", "OSJ"),
                  function(s) list(file = paste0(s, "_data.csv"), limits = param_limits, extra = cdec_2026(s))),
           c("SJR", "MSD", "PPT", "TRN", "MOK", "OSJ")),
  setNames(lapply(c("TWI", "ANH", "BLP"),
                  function(s) list(file = paste0(s, "_data.csv"), limits = tidal_cond, extra = cdec_2026(s))),
           c("TWI", "ANH", "BLP")),
  list(RRI_cdec_2026 = list(file = cdec_2026("RRI"), limits = param_limits, extra = character()))
)

qc_store <- function(...) here("data", "qc", ...)

# Run the QC for each job and store the result in data/qc/<name>_qc.rds (gitignored). Returns a one-row-
# per-series summary, also written to data/qc/qc_summary.csv. `only` limits the run to some jobs by name.
# Extra files that are not in data/ are skipped.
run_station_qc <- function(jobs = qc_jobs, only = names(jobs)) {
  dir.create(qc_store(), showWarnings = FALSE, recursive = TRUE)
  out <- imap(jobs[only], function(job, nm) {
    extra <- file.path(here("data"), job$extra)
    qc <- qc_single_sonde(here("data", job$file), job$limits, extra[file.exists(extra)])
    saveRDS(qc, qc_store(paste0(nm, "_qc.rds")))
    has_turb <- !is.na(qc$turb_fnu)
    tibble(series = nm, pct_valid = round(100 * mean(qc$qaqc_flag == "Valid", na.rm = TRUE), 1),
           first = as_date(min(qc$date_time[has_turb])), last = as_date(max(qc$date_time[has_turb])),
           readings_with_turbidity = sum(has_turb))
  })
  summary <- bind_rows(out)
  write_csv(summary, qc_store("qc_summary.csv"))
  summary
}

# Load a stored QC result (run run-qc.R, or run_station_qc(), first)
load_station_qc <- function(name) {
  f <- qc_store(paste0(name, "_qc.rds"))
  if (!file.exists(f)) stop("No stored QC for ", name, " (", f, "). Run run-qc.R first.", call. = FALSE)
  readRDS(f)
}
