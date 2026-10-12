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
qc_single_sonde <- function(path, limits = param_limits, extra_paths = character(), erratic_jump = 5) {
  raw <- read_station_file(path) |> mutate(source = "primary")
  for (p in extra_paths) {
    raw <- bind_rows(raw, filter(read_station_file(p), date_time > max(raw$date_time)) |> mutate(source = "provisional"))
  }
  qc_station_data(raw, limits, erratic_jump)
}

# Adds the QC columns (qaqc_flag, turb_fnu_final, and the flag columns of each step) to a data frame
# in the schema of read_station_file(). Readings are flagged, not deleted, except for gross-range
# failures and isolated spikes, which are set to NA in turb_fnu_final.
#
# `erratic_jump` is the erratic-period threshold in FNU (flag_erratic_periods(): median 15-minute change over
# a day). The default of 5 was tuned on RRI, where the typical step is about 0.4 FNU, so it means "about
# ten times the station's normal noise". In the turbid, tidal stations (Chipps Island, Suisun Bay) the normal
# step is 2 to 5 FNU, and 5 flags most of the record. "auto" sets the threshold to ten times the station's
# own typical step, never below 5 (auto_erratic_jump()). The value used is kept as attr(, "erratic_jump").
qc_station_data <- function(raw, limits = param_limits, erratic_jump = 5) {
  if (identical(erratic_jump, "auto")) erratic_jump <- auto_erratic_jump(raw$turb_fnu)
  out <- raw |>
    check_gross_range() |>
    check_param_ranges(limits) |>
    flag_robust_outliers() |>
    flag_erratic_periods(jump_thresh = erratic_jump) |>
    flag_flat_lines() |>
    classify_run_duration() |>
    corroborate_sustained_excursions() |>
    finalize_qc_flags()
  attr(out, "erratic_jump") <- erratic_jump
  out
}

# Is any parameter a steady ramp through time? A river or tidal station's turbidity, conductivity, temperature
# and so on rise and fall but do not climb steadily for years, so a rank correlation with time near +1 or -1
# means a column was sorted or otherwise damaged (Chipps Island's specific conductance: 0.86). Long real
# declines show up moderately (turbidity at Old River at Tracy: -0.6). Returns one row per parameter with
# `suspect` TRUE at or above `threshold`; the data are only checked, not changed.
qc_integrity <- function(x, params = c("turb_fnu", "watertemp_c", "ph", "spcond_uscm", "do_mgl", "chlor_ugl"), threshold = 0.8) {
  t <- as.numeric(x$date_time)
  map_dfr(params, function(p) {
    ok <- !is.na(x[[p]])
    if (sum(ok) < 1000) return(NULL)
    tibble(param = p, n = sum(ok), rho_time = cor(x[[p]][ok], t[ok], method = "spearman"))
  }) |> mutate(suspect = abs(rho_time) >= threshold)
}

# The same check for every stored QC series (or those named in `series`)
check_stored_integrity <- function(series = names(qc_jobs)) {
  map_dfr(series, \(s) qc_integrity(load_station_qc(s)) |> mutate(series = s, .before = 1))
}

# Ten times the station's median absolute 15-minute change in turbidity, never below `floor`
auto_erratic_jump <- function(turb, multiple = 10, floor = 5) {
  max(floor, multiple * median(abs(diff(turb)), na.rm = TRUE))
}

# Station QC jobs and stored results ------------------------------------------------------------------

# Conductivity ceilings: freshwater stations keep param_limits (20-5000 uS/cm); the tidal, brackish west Delta
# stations get 30000; the Suisun Bay and Chipps Island stations reach seawater range (about 41000 at Martinez)
# and get 60000.
tidal_cond <- modifyList(param_limits, list(spcond_uscm = c(20, 30000)))
bay_cond   <- modifyList(param_limits, list(spcond_uscm = c(20, 60000)))

# One entry per QC'd series. `file` is the station file, relative to the project folder; `limits` its parameter
# range limits; `extra` further files whose later readings are appended (the provisional CDEC downloads of
# 2026). RRI has its own three-depth QC (R/rri_multidepth_qc.R); only its provisional CDEC 1 m download is here.
cdec_2026 <- function(s) file.path("data", paste0("CDEC_", s, "_turb_2026.csv"))
edi_file <- function(program, s) {
  file.path("raw-EDI-data", program, paste0(s, if (program == "C-EMP-CWQ") "_data.csv" else "_POR.csv"))
}
make_jobs <- function(stations, limits, file_fn, extra_fn = function(s) character(), erratic_jump = 5) {
  setNames(lapply(stations, function(s)
    list(file = file_fn(s), limits = limits, extra = extra_fn(s), erratic_jump = erratic_jump)), stations)
}
data_file <- function(s) file.path("data", paste0(s, "_data.csv"))

qc_jobs <- c(
  # The ten stations of the flow-turbidity analysis (files in data/, with the provisional CDEC 2026 readings)
  make_jobs(c("SJR", "MSD", "PPT", "TRN", "MOK", "OSJ"), param_limits, data_file, cdec_2026),
  make_jobs(c("TWI", "ANH", "BLP"), tidal_cond, data_file, cdec_2026),
  list(RRI_cdec_2026 = list(file = cdec_2026("RRI"), limits = param_limits, extra = character())),

  # Further stations, from the EDI downloads in raw-EDI-data/ (record to 2025-12-31, no 2026 readings yet)
  make_jobs(c("SRH"), param_limits, \(s) edi_file("C-EMP-CWQ", s)),                       # Sacramento River at Hood
  make_jobs(c("RVB", "SSI"), tidal_cond, \(s) edi_file("C-EMP-CWQ", s)),                  # Rio Vista, Sherman Island
  make_jobs(c("FRK"), param_limits, \(s) edi_file("C-EMP-CWQ", s)),                       # Franks Tract
  # Suisun Bay: turbid and tidal, so the erratic-period threshold scales with each station's noise.
  # MAL (Chipps Island) is left out: in the EDI file its specific conductance for 2008-2024 is a steady ramp from
  # about 130 to 17,500 uS/cm (the column looks sorted), and 2025 is different, so the record cannot be QC'd.
  make_jobs(c("MRZ", "GZL", "HON", "RYC"), bay_cond, \(s) edi_file("C-EMP-CWQ", s), erratic_jump = "auto"),
  make_jobs(c("BET", "HLT", "HOL", "OBI", "ORQ"), param_limits, \(s) edi_file("NCRO-Central-Delta", s)),
  make_jobs(c("FAL", "FCT", "TSL"), tidal_cond, \(s) edi_file("NCRO-Central-Delta", s)),
  make_jobs(c("GLC", "GLE", "MHO", "MRU", "MRX", "ODM", "OH1", "ORI", "ORM", "ORX", "TWA", "VCU", "WCI"),
            param_limits, \(s) edi_file("NCRO-South-Delta", s))
)

qc_store <- function(...) here("data", "qc", ...)

# Run the QC for each job and store the result in data/qc/<name>_qc.rds (gitignored). Returns a one-row-
# per-series summary, also written to data/qc/qc_summary.csv. `only` limits the run to some jobs by name.
# Extra files that do not exist are skipped.
run_station_qc <- function(jobs = qc_jobs, only = names(jobs)) {
  dir.create(qc_store(), showWarnings = FALSE, recursive = TRUE)
  out <- imap(jobs[only], function(job, nm) {
    extra <- if (length(job$extra)) here(job$extra) else character()
    qc <- qc_single_sonde(here(job$file), job$limits, extra[file.exists(extra)], job$erratic_jump %||% 5)
    saveRDS(qc, qc_store(paste0(nm, "_qc.rds")))
    has_turb <- !is.na(qc$turb_fnu)
    # Share of the readings that have a turbidity value and passed every test. (Rows with no turbidity value
    # carry qaqc_flag "Valid" because there was nothing to flag, so they are left out of the denominator.)
    tibble(series = nm, erratic_jump = round(attr(qc, "erratic_jump"), 1),
           pct_valid = round(100 * mean(qc$qaqc_flag[has_turb] == "Valid"), 1),
           first = as_date(min(qc$date_time[has_turb])), last = as_date(max(qc$date_time[has_turb])),
           readings_with_turbidity = sum(has_turb))
  })
  summary <- bind_rows(out)
  f <- qc_store("qc_summary.csv")
  all_rows <- if (file.exists(f)) bind_rows(summary, filter(read_csv(f, show_col_types = FALSE), !series %in% summary$series)) else summary
  write_csv(all_rows[order(match(all_rows$series, names(jobs))), ], f)
  summary
}

# Load a stored QC result (run run-qc.R, or run_station_qc(), first)
load_station_qc <- function(name) {
  f <- qc_store(paste0(name, "_qc.rds"))
  if (!file.exists(f)) stop("No stored QC for ", name, " (", f, "). Run run-qc.R first.", call. = FALSE)
  readRDS(f)
}
