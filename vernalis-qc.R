# Secondary QC for the continuous water quality record at Vernalis (station "SJR"), using the
# same functions as the RRI pipeline (R/qc_functions.R).
#
# Single sonde, so the cross-depth steps (2c depth drift, 2e parameter drift, 4b cross-depth
# corroboration) cannot run, and the RRI-specific summer ceiling of 100 FNU is NOT applied:
# river turbidity here legitimately reaches hundreds of FNU in storms.
# What does run: step 1 turbidity range, step 1b range screen for the other parameters,
# step 2 robust outliers, step 2b erratic periods, step 2d flat lines, step 3 spike vs sustained,
# step 4 cross-parameter corroboration (conductivity and DO), step 5 final flags.
#
# Needs data/SJR_data.csv. Its clock was checked against daily temperature, DO and pH cycles
# and is true Pacific Standard Time (unlike the RRI file, which needed an 8 h correction).
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))

# Rename to the column names the QC functions expect. depth_m is a constant 1: one sonde.
# `fluorescence` is treated as chlorophyll in ug/L (the 0-400 range screen assumes this;
# values run 0.4-174, similar to the RRI chlorophyll).
vernalis <- read_csv(here("data", "SJR_data.csv"), show_col_types = FALSE) |>
  transmute(depth_m = 1,
            date_time = ymd_hms(paste(date, time), tz = "Etc/GMT+8"),
            year = year(date_time), month = month(date_time),
            do_mgl = dissolvedoxygen, chlor_ugl = fluorescence, ph,
            spcond_uscm = spc, turb_fnu = turbidity, watertemp_c = watertemperature)

vern_qc <- vernalis |>
  check_gross_range() |>
  check_param_ranges() |>
  flag_robust_outliers() |>
  flag_erratic_periods() |>
  flag_flat_lines() |>
  classify_run_duration() |>
  corroborate_sustained_excursions()

# finalize_qc_flags() expects the depth-based columns; with one sonde they are all "no signal"
vern_qc$drift_flag <- NA_character_
vern_qc$depth_corroborated <- FALSE
vern_qc <- finalize_qc_flags(vern_qc)

# Analysis column: removes the same things as the RRI analysis column except the summer ceiling
vern_qc <- vern_qc |>
  mutate(analysis_exclusion = if_else(
    qaqc_flag %in% c("Flagged: Erratic Period (review)", "Flagged: Flat Line (review)",
                     "Flagged: Uncorroborated Excursion (review)"), qaqc_flag, NA_character_),
    turb_fnu_analysis = if_else(is.na(analysis_exclusion), turb_fnu_final, NA_real_))

vern_qc |>
  count(qaqc_flag) |>
  mutate(pct = round(100 * n / sum(n), 2)) |>
  print()

vern_qc |>
  summarise(across(starts_with("range_flag_"), ~ sum(.x == "Fail: Range", na.rm = TRUE))) |>
  print()

vern_qc |>
  group_by(year) |>
  summarise(erratic = sum(erratic_flag == "Fail: Erratic Period", na.rm = TRUE),
            flat = sum(flat_flag == "Fail: Flat Line", na.rm = TRUE),
            isolated_spikes = sum(duration_class == "Isolated Spike", na.rm = TRUE),
            uncorroborated = sum(qaqc_flag == "Flagged: Uncorroborated Excursion (review)", na.rm = TRUE)) |>
  print(n = Inf)
