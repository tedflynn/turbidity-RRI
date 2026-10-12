# Which EDI stations have enough turbidity to anchor a pre-detection baseline?
#
# Reads data/edi_summer_days.csv (from edi-inventory.R: summer days, Jul-Sep, with at least 48 turbidity
# readings, per station and year; raw readings, before QC) and counts, per station, the usable pre-detection
# summers (a summer is usable with at least `min_days_per_summer` such days, the same rule as the analysis),
# whether they include both wet and dry years, and whether 2025 is usable. Also compares the EDI copy of
# the stations that were already in data/ with those files, to confirm they are the same record.
#
# Output: data/edi_summer_screen.csv.
library(tidyverse)
library(lubridate)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))

min_days_per_summer <- 60
last_pre_summer <- 2024           # detection was October 2024, so the 2025 summer is the first post-detection one

summer <- read_csv(here("data", "edi_summer_days.csv"), show_col_types = FALSE)
inventory <- read_csv(here("data", "edi_inventory.csv"), show_col_types = FALSE)
wy <- read_csv(here("data", "WY-type.csv"), show_col_types = FALSE) |>
  select(year = 1, wy_type = last_col()) |> mutate(wet = wy_type %in% c("W", "AN"))

usable <- summer |>
  filter(valid_summer_days >= min_days_per_summer) |>
  left_join(wy, by = "year")

longest_run <- function(y) { if (length(y) == 0) return(0L); max(rle(cumsum(c(1, diff(sort(y)) != 1)))$lengths) }

screen <- inventory |>
  select(station, program, turb_first, turb_last, pct_turbidity) |>
  left_join(
    usable |> filter(year <= last_pre_summer) |>
      group_by(station) |>
      summarise(usable_pre_summers = n(), first_usable = min(year), wet_summers = sum(wet, na.rm = TRUE),
                dry_summers = sum(!wet, na.rm = TRUE), longest_run = longest_run(year), .groups = "drop"),
    by = "station") |>
  left_join(usable |> filter(year == 2025) |> transmute(station, usable_2025 = TRUE), by = "station") |>
  mutate(across(c(usable_pre_summers, wet_summers, dry_summers, longest_run), ~ replace_na(.x, 0L)),
         usable_2025 = replace_na(usable_2025, FALSE),
         both_wet_and_dry = wet_summers >= 2 & dry_summers >= 2) |>
  arrange(desc(usable_pre_summers))
write_csv(screen, here("data", "edi_summer_screen.csv"))

cat("Usable pre-detection summers (", min_days_per_summer, "+ days with turbidity), 2024 and earlier:\n", sep = "")
print(screen |> select(-program, -turb_last), n = Inf, width = Inf)

cat("\nSummers that fall short, for stations with at least 8 usable ones (valid days in the summer):\n")
summer |>
  filter(station %in% screen$station[screen$usable_pre_summers >= 8], year >= 2008, year <= 2025, valid_summer_days < min_days_per_summer) |>
  arrange(station, year) |> print(n = Inf)

# Are the EDI copies the same record as the files already in data/? -----------------------------------------
cat("\nEDI copy against the station file already in data/ (turbidity at matching timestamps):\n")
same <- map_dfr(c("SJR", "MSD", "PPT", "TWI", "ANH", "BLP", "TRN", "MOK", "OSJ"), function(s) {
  old <- read_station_file(here("data", paste0(s, "_data.csv"))) |> select(date_time, old = turb_fnu)
  p <- inventory$program[inventory$station == s]
  f <- list.files(here("raw-EDI-data", p), pattern = paste0("^", s, "_"), full.names = TRUE)
  new <- read_station_file(f) |> select(date_time, new = turb_fnu)
  j <- inner_join(old, new, by = "date_time")
  both <- filter(j, !is.na(old), !is.na(new))
  tibble(station = s, rows_old = nrow(old), rows_edi = nrow(new), matched = nrow(j), both_present = nrow(both),
         identical = sum(both$old == both$new), max_abs_diff = max(abs(both$old - both$new)),
         only_in_edi = sum(!new$date_time %in% old$date_time & !is.na(new$new)),
         only_in_old = sum(!old$date_time %in% new$date_time & !is.na(old$old)))
})
print(same, width = Inf)
