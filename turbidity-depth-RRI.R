# Load required packages
library(tidyverse)
library(lubridate)
library(here)
library(slider)
library(oddwater)
library(imputeTS)

# Import sonde data ------------------------------------------------------------
# Data was QC's and provided by Andrew Tran in C-EMP 
df_1m <- read_csv(here("data","RRI_data_1m.csv")) |> mutate(depth_m = 1)
df_3m <- read_csv(here("data","RRI_data_3m.csv")) |> mutate(depth_m = 3) 
df_6m <- read_csv(here("data","RRI_data_6m.csv")) |> mutate(depth_m = 6)

# Combine data from all three depths
wq_data <- bind_rows(df_1m, df_3m, df_6m)

# Add temporal and water year type data
wq_data <- wq_data |>
  mutate(
    # all data is timestamped to Pacific Standard Time (fixed -08:00, no DST).
    # read_csv() already parses the "-08:00" offset in the raw strings correctly
    # into a UTC-labeled POSIXct, so we only need to relabel the tzone here, not
    # reparse the clock digits (ymd_hms() would shift every timestamp by 8 hours).
    date_time = with_tz(date_time, tzone = "Etc/GMT+8"),
    year = year(date_time),
    month = month(date_time),
    day = day(date_time),
    water_year = if_else(month >= 10, year + 1, year)
    # water year begins October 1st
  )

# Add table of WY types downloaded from CDEC and select SJR version
sj_wy_types <- read_csv(here("data","WY-type.csv")) |> 
  select("Year","Yr-type-SJR") |> 
  rename("wy_type" = "Yr-type-SJR") |> 
  rename("water_year" = "Year")

# Join WY metadata with main data frame
wq_data <- wq_data |>
  left_join(sj_wy_types, by = "water_year") |>
  # Set factors so graphs plot in logical wetness order rather than alphabetical
  mutate(wy_type = factor(wy_type, levels = c("W", "AN", "BN", "D", "C")))

# Basic data QC
# Remove any erroneous negative sonde turbidity readings
wq_clean <- wq_data |> filter(turb_fnu >= 0 | is.na(turb_fnu))

# Secondary (automated) QC -------------------------------------------------
# Field-based primary QC (Wagner et al. 2006 method) has already been applied
# upstream. These steps catch residual statistical outliers. Each step adds
# columns without dropping rows, so you can inspect wq_qc1, wq_qc2, etc.
# individually before trusting the final result. See R/qc_functions.R for
# the method behind each step and its tunable thresholds.
source(here("R", "qc_functions.R"))

# Step 1: gross/physical range check
wq_qc1 <- check_gross_range(wq_clean, min_val = 0, max_val = 4000)
wq_qc1 |> count(gross_range_flag)

# Step 1b: range screen for temperature, pH, conductivity, DO and chlorophyll
# (impossible values -> NA in <param>_qc columns; limits in param_limits)
wq_qc1 <- check_param_ranges(wq_qc1)
wq_qc1 |>
  pivot_longer(starts_with("range_flag_"), names_to = "param", values_to = "flag") |>
  count(param, depth_m, flag) |> filter(flag == "Fail: Range")

# Step 2: robust local outlier detection (modified z-score on rolling median/MAD)
wq_qc2 <- flag_robust_outliers(wq_qc1, window = 12, z_thresh = 3.5, noise_floor = 5)
wq_qc2 |> count(spike_flag)

# Step 2b: flag erratic periods (sustained disturbances longer than step 2's
# rolling window, which it can't see). Flagged for review, not removed.
wq_qc2b <- flag_erratic_periods(wq_qc2, window = 48, jump_thresh = 5)
wq_qc2b |> count(depth_m, erratic_flag)

# Step 2c: flag cross-depth drift (one sensor climbing away from the other two,
# e.g. biofouling that resets at a service visit). Flagged for review, not removed.
wq_qc2c <- flag_depth_drift(wq_qc2b, diff_thresh = 8, ratio_thresh = 3)
wq_qc2c |> count(depth_m, drift_flag)

# Step 2d: flag flat lines (12+ identical readings in a row = stuck sensor).
# Flagged for review, not removed.
wq_qc2d <- flag_flat_lines(wq_qc2c, min_run = 12, pad_readings = 24)
wq_qc2d |> count(depth_m, flat_flag)

# Step 2e: flag each of the other parameters drifting away from the other two
# depths (24 hr median vs the median of all three depths, sustained 48+ hrs).
# Adds <prefix>_drift_flag columns for any analysis of those parameters; this
# does not change turbidity decisions. Thresholds sit beyond the 99th percentile
# of normal depth-to-depth differences (see flag_param_depth_drift()).
drift_settings <- list(
  ph   = list(param = "ph_qc",          diff_thresh = 0.5, ratio_thresh = 1),
  temp = list(param = "watertemp_c_qc", diff_thresh = 1.0, ratio_thresh = 1),
  cond = list(param = "spcond_uscm_qc", diff_thresh = 100, ratio_thresh = 1.25),
  do   = list(param = "do_mgl_qc",      diff_thresh = 1.5, ratio_thresh = 1),
  chl  = list(param = "chlor_ugl_qc",   diff_thresh = 5,   ratio_thresh = 2)
)
wq_qc2e <- wq_qc2d
for (p in names(drift_settings)) {
  s <- drift_settings[[p]]
  wq_qc2e <- flag_param_depth_drift(wq_qc2e, param = s$param, prefix = p,
                                    diff_thresh = s$diff_thresh, ratio_thresh = s$ratio_thresh)
}
wq_qc2e |>
  pivot_longer(ends_with("_drift_flag") & !starts_with("drift"), names_to = "param", values_to = "flag") |>
  count(param, depth_m, flag) |> filter(flag == "Fail: Depth Drift")

# Step 3: classify flagged points as isolated spikes vs sustained excursions
wq_qc3 <- classify_run_duration(wq_qc2e, max_spike_len = 2)
wq_qc3 |> count(duration_class)

# Step 4: check sustained excursions for corroborating signal in other channels
wq_qc4 <- corroborate_sustained_excursions(wq_qc3, window = 12, corrob_z_thresh = 2)
wq_qc4 |> filter(duration_class == "Sustained Excursion") |> count(run_corroborated)

# Step 4b: check sustained excursions for a coincident signal at the other
# two depths (same housing, same site -- a real event should show up in
# turbidity at more than one depth even if conductivity/DO don't move)
wq_qc4b <- corroborate_cross_depth(wq_qc4, z_thresh = 3.5, pad_readings = 1)
wq_qc4b |> filter(duration_class == "Sustained Excursion") |> count(run_corroborated, depth_corroborated)

# Step 5: finalize QC decision (removes gross-range fails and isolated spikes;
# flags sustained excursions corroborated by neither other channels nor other
# depths for manual review rather than removing them)
wq_qc_final <- finalize_qc_flags(wq_qc4b)
wq_qc_final |> count(qaqc_flag)

# Step 6: analysis-ready turbidity for the Golden Mussel clarity comparison.
# On top of turb_fnu_final, excludes the "Flagged: ... (review)" points
# (excursions, erratic periods, depth drift) and
# any summer episode above 100 FNU (+/- 6 hrs). These are most likely
# vessel-wake/dredging resuspension or instrument problems -- real but
# non-ambient conditions that confound a clarity comparison. Kept as a
# separate column (rather than overwriting turb_fnu_final) so a with/without
# sensitivity check is still possible later; analysis_exclusion records why.
wq_qc_final <- prepare_analysis_turbidity(wq_qc_final, ceiling = 100, pad_readings = 24)
wq_qc_final |> filter(month %in% 7:9) |> count(analysis_exclusion)

# 5. SUMMER TURBIDITY ISOLATION & INVASION FLAG
summer_turb <- wq_qc_final %>%
  # Isolate July, August, September
  filter(month %in% c(7, 8, 9))

# 6. VISUALIZATION
# Uses turb_fnu_analysis (flagged excursions excluded) -- this is the actual
# Golden Mussel clarity comparison.
ggplot(summer_turb, aes(x = factor(water_year), y = turb_fnu_analysis, fill = wy_type)) +
  geom_boxplot(outlier.alpha = 0.3,
               outlier.shape = 21) +
  # space = "free_x" prevents years from artificially stretching in WY types with fewer occurrences
  facet_grid(depth_m ~ ., scales = "free_x", space = "free_x",
             labeller = labeller(depth_m = function(x) paste(x, "m Depth"))) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    title = "San Joaquin River Summer Turbidity (July - September)",
    subtitle = "Assessing water clarity impacts following the detection of Golden Mussels",
    x = "Water Year",
    y = "Turbidity (FNU)",
    fill = "Water Year Type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot continuous summer time series, faceted by year
# Extract the unique depth values from your dataset (1, 3, and 6)
depth_levels <- unique(summer_turb$depth_m)

# Loop through each depth
for (current_depth in depth_levels) {
  
  # Filter the data for the active depth in the loop
  plot_data <- summer_turb |> 
    filter(depth_m == current_depth)
  
  # Build the faceted time-series plot
  # Uses turb_fnu_final (flagged excursions still visible) rather than
  # turb_fnu_analysis, since the point of this plot is to visually screen
  # for fouling/sensor issues -- hiding flagged events would defeat that.
  p <- ggplot(plot_data, aes(x = date_time, y = turb_fnu_final)) +
    geom_line(linewidth = 0.3) +
    # Facet by year to isolate specific deployments and artifacts
    facet_wrap(~ water_year, scales = "free_x") +
    theme_bw() +
    labs(
      title = paste("Summer Turbidity Time Series -", current_depth, "m Depth"),
      subtitle = "Visual QC screen for continuous fouling or sensor failure",
      x = "Date",
      y = "Turbidity (FNU)"
    )
  
  # Force the plot to render
  print(p)
  
  # Optional: Save each plot automatically to your working directory
  ggsave(paste0("QC_Turbidity_",current_depth,"m.png"), 
         plot = p, 
         path = here("plots"), 
         width = 10, 
         height = 6)
}
