# Summaries of the stored QC results (data/qc/<series>_qc.rds, from run-qc.R) for the QC report
# (reports/qc-report.qmd). Each series is read once and reduced to small tables, so ten stations of
# 15-minute data (about 7 million readings) never need to sit in memory together.
#
# Requires R/station_qc.R (load_station_qc(), qc_jobs) and the tidyverse.

# Outcome categories, in the order they are stacked and coloured in the report. "No turbidity value" is
# not a QC decision: the sonde logged nothing (or the reading was blank), and finalize_qc_flags() calls
# those rows "Valid" because there is nothing to flag. The report separates them so that "valid" always
# means a turbidity reading that passed every test.
qc_outcomes <- tribble(
  ~outcome,                                       ~label,                      ~colour,
  "Valid",                                        "Valid",                     "#B8C4CE",
  "Kept: Corroborated Event",                     "Corroborated event (kept)", "#4C9F70",
  "Flagged: Uncorroborated Excursion (review)",   "Uncorroborated excursion",  "#E69F00",
  "Flagged: Flat Line (review)",                  "Flat line",                 "#0072B2",
  "Flagged: Erratic Period (review)",             "Erratic period",            "#CC79A7",
  "Removed: Isolated Spike",                      "Isolated spike (removed)",  "#D55E00",
  "Removed: Gross Range",                         "Gross range (removed)",     "#7A0403",
  "No turbidity value",                           "No turbidity value",        "#F2F2F2"
)
qc_outcome_label <- function(x) factor(x, levels = qc_outcomes$outcome, labels = qc_outcomes$label)
qc_outcome_colours <- setNames(qc_outcomes$colour, qc_outcomes$label)

# Parameters screened: column in the station file, the range-flag column, and a display label
qc_params <- tribble(
  ~param,        ~flag_col,                  ~label,
  "turb_fnu",    "gross_range_flag",         "Turbidity",
  "watertemp_c", "range_flag_watertemp_c",   "Water temperature",
  "ph",          "range_flag_ph",            "pH",
  "spcond_uscm", "range_flag_spcond_uscm",   "Specific conductance",
  "do_mgl",      "range_flag_do_mgl",        "Dissolved oxygen",
  "chlor_ugl",   "range_flag_chlor_ugl",     "Chlorophyll"
)

# Runs of consecutive readings with the same flagged outcome, one row per run (an "event"). Isolated spikes
# are included (they are runs of one or two readings). `baseline_fnu` is the median of the valid turbidity
# readings in the 24 hours before the run, to show how far the event stood above normal.
qc_events_one <- function(x, name) {
  x <- x |> arrange(date_time)
  runs <- x |>
    mutate(run = consecutive_id(outcome), i = row_number()) |>
    group_by(run) |>
    summarise(outcome = first(outcome), start = first(date_time), end = last(date_time),
              first_i = first(i), n_readings = n(), peak_fnu = suppressWarnings(max(turb_fnu, na.rm = TRUE)),
              median_fnu = median(turb_fnu, na.rm = TRUE), .groups = "drop") |>
    filter(!outcome %in% c("Valid", "No turbidity value"))
  if (nrow(runs) == 0) return(tibble())
  ok <- x$outcome == "Valid"
  turb <- x$turb_fnu
  runs$baseline_fnu <- map_dbl(runs$first_i, function(i) {
    idx <- max(1, i - 96):max(1, i - 1)
    median(turb[idx][ok[idx]], na.rm = TRUE)
  })
  runs |> mutate(series = name, hours = n_readings * 0.25, .before = 1) |> select(-run, -first_i)
}

# Reduce one stored QC result to the tables the report uses
qc_summarise_series <- function(name) {
  x <- load_station_qc(name) |>
    mutate(outcome = if_else(is.na(turb_fnu), "No turbidity value", qaqc_flag),
           date = as_date(date_time))
  if (!"source" %in% names(x)) x$source <- "primary"

  # A file can start long before its turbidity record (BLP's begins in 1997, its turbidity in 2010). The
  # turbidity outcomes below cover the period with turbidity; the range screen of the other sensors uses
  # the whole file.
  has_turb <- !is.na(x$turb_fnu)
  x_all <- x
  x <- x |> filter(date_time >= min(date_time[has_turb]), date_time <= max(date_time[has_turb]))

  expected <- as.integer(floor(as.numeric(difftime(max(x$date_time), min(x$date_time), units = "mins")) / 15)) + 1L
  scope <- tibble(series = name, first = as_date(min(x$date_time)), last = as_date(max(x$date_time)),
                  rows = nrow(x), expected_rows = expected, with_turbidity = sum(!is.na(x$turb_fnu)),
                  provisional_rows = sum(x$source == "provisional"))

  outcome <- x |> count(year, month, source, outcome) |> mutate(series = name, .before = 1)

  # Range screens: how many readings each sensor reported, and how many fell outside the physical limits
  params <- map_dfr(seq_len(nrow(qc_params)), function(k) {
    p <- qc_params$param[k]; f <- qc_params$flag_col[k]
    tibble(year = x_all$year, measured = !is.na(x_all[[p]]),
           fail = !is.na(x_all[[f]]) & x_all[[f]] %in% c("Fail: Range", "Fail: Gross Range")) |>
      group_by(year) |>
      summarise(measured = sum(measured), range_fail = sum(fail), .groups = "drop") |>
      mutate(series = name, param = p, .before = 1)
  })

  # Could a sustained excursion have been corroborated at all? Needs conductivity or DO at that time.
  exc <- x |> filter(duration_class == "Sustained Excursion")
  corroboration <- tibble(series = name, excursion_readings = nrow(exc),
                          with_channels = sum(!is.na(exc$spcond_uscm_qc) | !is.na(exc$do_mgl_qc)),
                          corroborated = sum(exc$run_corroborated %in% TRUE))

  daily <- x |>
    group_by(date) |>
    summarise(readings_with_turbidity = sum(!is.na(turb_fnu)),
              median_final = if (sum(!is.na(turb_fnu_final)) >= 48) median(turb_fnu_final, na.rm = TRUE) else NA_real_,
              .groups = "drop") |>
    mutate(series = name, .before = 1)

  flag_days <- x |>
    filter(!outcome %in% c("Valid", "No turbidity value", "Kept: Corroborated Event")) |>
    group_by(date, outcome) |>
    summarise(n = n(), peak_fnu = max(turb_fnu, na.rm = TRUE), .groups = "drop") |>
    mutate(series = name, .before = 1)

  list(scope = scope, outcome = outcome, params = params, corroboration = corroboration,
       daily = daily, flag_days = flag_days, events = qc_events_one(x, name))
}

# All series: a list with one bound table per element. `series` defaults to every stored QC result.
qc_report_tables <- function(series = names(qc_jobs)) {
  parts <- map(series, qc_summarise_series)
  names(parts) <- series
  map(set_names(names(parts[[1]])), \(tbl) bind_rows(map(parts, tbl)))
}
