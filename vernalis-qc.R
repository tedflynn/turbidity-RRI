# Secondary QC for the continuous water quality record at Vernalis (station "SJR"), using the
# same functions as the RRI pipeline (R/qc_functions.R) through qc_single_sonde() in R/station_qc.R.
#
# Single sonde, so the cross-depth steps cannot run, and the RRI-specific summer ceiling of
# 100 FNU is not applied: river turbidity here legitimately reaches hundreds of FNU in storms.
#
# Needs data/SJR_data.csv. Its clock was checked against daily temperature, DO and pH cycles
# and is true Pacific Standard Time (unlike the RRI file, which needed an 8 h correction).
# `fluorescence` is treated as chlorophyll in ug/L (confirmed).
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))

vern_qc <- qc_single_sonde(here("data", "SJR_data.csv"))

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
