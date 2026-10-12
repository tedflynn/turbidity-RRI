# Quality control for the continuous water quality stations, separate from any analysis.
#
# Runs the general single-sonde QC (R/qc_functions.R, chained in R/station_qc.R) for every station in
# `qc_jobs` (R/station_qc.R) and stores each result in data/qc/<series>_qc.rds, plus a one-line-per-series
# summary in data/qc/qc_summary.csv. Nothing is deleted: every reading is kept with a qaqc_flag and the
# flag columns of each step; analysis scripts decide what to use.
#
# Run this after adding or changing a data file or any QC setting, then rerun the analysis scripts. To add
# a station, put its file in data/, add an entry to `qc_jobs` (with parameter range limits suited to the
# site) and run this script. The RRI three-depth record has its own QC (R/rri_multidepth_qc.R) and is
# stored the same way by rri_qc_final(refresh = TRUE).
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))

qc_summary <- run_station_qc()
cat("QC outcome by series (pct_valid = share of readings with no flag or removal):\n")
print(qc_summary)

# Flag counts by station
cat("\nReadings by QC flag:\n")
qc_jobs |>
  names() |>
  map_dfr(~ load_station_qc(.x) |> count(qaqc_flag) |> mutate(series = .x)) |>
  pivot_wider(names_from = qaqc_flag, values_from = n, values_fill = 0) |>
  print(width = Inf)

# A report of these results (counts, coverage, events, review queue, per-station plots) is in reports/qc-report.qmd:
#   quarto render reports/qc-report.qmd

# File integrity: any parameter that rises or falls steadily with time looks like a sorted or damaged column
cat("\nParameters that are a steady ramp through time (|Spearman rho with time| >= 0.8); none expected:\n")
integrity <- check_stored_integrity()
print(filter(integrity, suspect))
