# Compare continuous sonde data (RRI) with monthly discrete grab samples (EMP station P8).
# Independent check on sensor calibration: sensor-vs-sensor checks cannot see a bias
# shared by all three depths.
#
# Needs data/DWQ_DWR_EDI_2025.csv (EDI package edi.458.15, downloaded manually; gitignored)
# and the QC output from rri_qc_final() in R/rri_multidepth_qc.R (loaded below if not already in memory).
library(tidyverse)
library(lubridate)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "rri_multidepth_qc.R"))
source(here("R", "mussel_analysis.R"))
if (!exists("wq_qc_final")) wq_qc_final <- rri_qc_final()

grab_station <- "P8"
tz <- "Etc/GMT+8"                 # both records are in Pacific Standard Time (assumed for grab data)
window_hrs <- 2                   # sonde median taken within +/- this many hours of the grab sample
min_readings <- 4                 # require at least this many sonde readings in the window

# Grab samples -----------------------------------------------------------------
# Drops QC-flagged results, non-detects and text entries (weather observations, etc.)
grab <- read_csv(here("data", "DWQ_DWR_EDI_2025.csv"), show_col_types = FALSE,
                 col_types = cols(.default = col_character())) |>
  mutate(Date = as_date(Date), val = suppressWarnings(as.numeric(Result_Value))) |>
  filter(Station == grab_station, Date >= as.Date("2007-12-25"),
         is.na(QC_Flag), !is.na(val), Detection_Status == "Detected") |>
  mutate(t = ymd_hms(paste(Date, Time), tz = tz)) |>
  filter(!is.na(t), !(Analyte == "Turbidity" & Result_Unit != "FNU"))

# Which sonde column corresponds to each grab analyte, and which flag marks a doubtful sonde reading
analyte_map <- tribble(
  ~Analyte,                                 ~param,                                  ~grab_type, ~sonde_col,        ~flag_col,
  "Specific Conductance",                   "Conductivity (uS/cm)",                  "surface",  "spcond_uscm_qc",  "cond_drift_flag",
  "(Bottom) Specific Conductance",          "Conductivity (uS/cm)",                  "bottom",   "spcond_uscm_qc",  "cond_drift_flag",
  "Water Temperature",                      "Temperature (C)",                       "surface",  "watertemp_c_qc",  "temp_drift_flag",
  "(Bottom) Water Temperature",             "Temperature (C)",                       "bottom",   "watertemp_c_qc",  "temp_drift_flag",
  "Dissolved Oxygen mg/L",                  "DO (mg/L)",                             "surface",  "do_mgl_qc",       "do_drift_flag",
  "(Bottom) Dissolved Oxygen mg/L",         "DO (mg/L)",                             "bottom",   "do_mgl_qc",       "do_drift_flag",
  "pH",                                     "pH",                                    "surface",  "ph_qc",           "ph_drift_flag",
  "(Bottom) pH",                            "pH",                                    "bottom",   "ph_qc",           "ph_drift_flag",
  "Turbidity",                              "Turbidity (FNU)",                       "surface",  "turb_fnu_qc",     "turb_flag",
  "(Bottom) Turbidity",                     "Turbidity (FNU)",                       "bottom",   "turb_fnu_qc",     "turb_flag",
  "Chlorophyll Fluorescence ug/L",          "Chlorophyll: grab fluorescence (ug/L)", "surface",  "chlor_ugl_qc",    "chl_drift_flag",
  "(Bottom) Chlorophyll Fluorescence ug/L", "Chlorophyll: grab fluorescence (ug/L)", "bottom",   "chlor_ugl_qc",    "chl_drift_flag",
  "Chlorophyll a",                          "Chlorophyll: lab chl-a (ug/L)",         "surface",  "chlor_ugl_qc",    "chl_drift_flag"
)
grab <- inner_join(grab, analyte_map, by = "Analyte")

# Sonde series, one list element per depth, with a logical "doubtful reading" column per parameter
sonde <- wq_qc_final |>
  mutate(turb_flag = !is.na(analysis_exclusion),
         across(ends_with("_drift_flag"), ~ !is.na(.x) & .x == "Fail: Depth Drift")) |>
  select(depth_m, date_time, ends_with("_qc"), ends_with("_drift_flag"), turb_flag) |>
  arrange(depth_m, date_time)
sonde_by_depth <- split(sonde, sonde$depth_m)

pair_one <- function(t, sonde_col, flag_col, depth) {
  x <- sonde_by_depth[[as.character(depth)]]
  tt <- as.numeric(x$date_time)
  i0 <- findInterval(as.numeric(t) - window_hrs * 3600 - 1, tt) + 1
  i1 <- findInterval(as.numeric(t) + window_hrs * 3600, tt)
  if (i1 < i0) return(list(sonde = NA_real_, flagged = NA))
  v <- x[[sonde_col]][i0:i1]
  list(sonde = if (sum(!is.na(v)) >= min_readings) median(v, na.rm = TRUE) else NA_real_,
       flagged = any(x[[flag_col]][i0:i1]))
}

# Pair every grab value with each sonde depth ---------------------------------------
pairs <- map_dfr(c(1, 3, 6), function(d) {
  r <- pmap(list(grab$t, grab$sonde_col, grab$flag_col), pair_one, depth = d)
  grab |> mutate(sonde_depth = d, sonde = map_dbl(r, "sonde"), flagged = map_lgl(r, "flagged"))
}) |>
  filter(!is.na(sonde),
         (grab_type == "surface" & sonde_depth %in% c(1, 3)) | (grab_type == "bottom" & sonde_depth == 6)) |>
  transmute(param, grab_type, sonde_depth, t, year = year(t), grab = val, sonde, diff = sonde - grab, flagged)

# Agreement statistics ------------------------------------------------------------------
# Surface grab vs 1m and 3m sondes; bottom grab (about 11 m at P8) vs the 6m sonde.
agreement <- pairs |>
  group_by(param, grab_type, sonde_depth) |>
  summarise(n = n(), grab_median = median(grab), sonde_median = median(sonde),
            median_diff = median(diff), mad_diff = mad(diff),
            median_ratio = median(sonde / grab, na.rm = TRUE),
            spearman = cor(sonde, grab, method = "spearman"), .groups = "drop")
print(agreement, n = Inf, width = Inf)

# Does the sonde disagree more with the grab samples when our QC flagged that period?
# (Monthly samples rarely land inside a flagged period, so this is a small-sample check.)
pairs |>
  group_by(param, flagged) |>
  summarise(n = n(), median_abs_diff = median(abs(diff)), .groups = "drop") |>
  print(n = Inf)

# Plots ---------------------------------------------------------------------------------
param_levels <- c("Conductivity (uS/cm)", "Temperature (C)", "DO (mg/L)", "pH",
                  "Turbidity (FNU)", "Chlorophyll: grab fluorescence (ug/L)")
plot_pairs <- pairs |>
  filter(param %in% param_levels) |>
  mutate(param = factor(param, levels = param_levels),
         pair = if_else(grab_type == "surface", paste0("surface grab vs ", sonde_depth, "m sonde"), "bottom grab vs 6m sonde"))

p_scatter <- ggplot(plot_pairs, aes(grab, sonde, color = flagged)) +
  geom_abline(color = "grey50") +
  geom_point(alpha = 0.5, size = 0.9) +
  facet_grid(param ~ pair, scales = "free") +
  scale_color_manual(values = c("FALSE" = "#2c7fb8", "TRUE" = "#d95f0e"), name = "sonde inside a\nflagged period") +
  theme_bw(base_size = 9) +
  labs(title = paste("Sonde vs", grab_station, "grab samples"), x = "Grab sample", y = "Sonde (RRI)")

# pH, DO and temperature as a difference; the others as log2(sonde / grab)
p_time <- plot_pairs |>
  mutate(y = if_else(grepl("pH|DO|Temp", param), diff, log2(pmax(sonde, 0.05) / pmax(grab, 0.05)))) |>
  ggplot(aes(t, y, color = factor(sonde_depth))) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_point(alpha = 0.55, size = 0.9) +
  facet_wrap(~ param, ncol = 1, scales = "free_y") +
  scale_color_brewer(palette = "Set1", name = "sonde\ndepth (m)") +
  theme_bw(base_size = 9) +
  labs(title = "Sonde minus grab over time (positive = sonde higher)", x = NULL,
       y = "difference (pH, DO, temp) or log2 ratio (others)")

ggsave(here("plots", "grab_vs_sonde_scatter.png"), p_scatter, width = 9, height = 11, dpi = 130)
ggsave(here("plots", "grab_vs_sonde_over_time.png"), p_time, width = 9, height = 12, dpi = 130)
