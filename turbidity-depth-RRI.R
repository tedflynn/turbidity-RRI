library(tidyverse)
library(lubridate)
library(here)
library(slider)
library(oddwater)

# Import sonde data ------------------------------------------------------------
# Data was QC's and provided by Andrew Tran in C-EMP 
df_1m <- read_csv(here("turbidity-depth-RRI","data","RRI_data_1m.csv")) %>% 
  mutate(depth_m = 1)
df_3m <- read_csv(here("turbidity-depth-RRI","data","RRI_data_3m.csv")) %>% 
  mutate(depth_m = 3) 
df_6m <- read_csv(here("turbidity-depth-RRI","data","RRI_data_6m.csv")) %>% 
  mutate(depth_m = 6)

# Combine data from all three depths
wq_data <- bind_rows(df_1m, df_3m, df_6m)

# Add temporal and water year type data
wq_data <- wq_data %>%
  mutate(
    # all data is timestamped to Pacific Standard Time
    date_time = ymd_hms(date_time, tz = "Etc/GMT+8", truncated = 3), 
    # truncated argument above prevents failure on 00:00:00
    year = year(date_time),
    month = month(date_time),
    day = day(date_time),
    water_year = if_else(month >= 10, year + 1, year)
    # water year begins October 1st
  )

# Add table of WY types downloaded from CDEC and select SJR version
sj_wy_types <- read_csv(here("turbidity-depth-RRI","data","WY-type.csv")) |> 
  select("Year","Yr-type-SJR") |> 
  rename("wy_type" = "Yr-type-SJR") |> 
  rename("water_year" = "Year")

# Join WY metadata with main data frame
wq_data <- wq_data %>%
  left_join(sj_wy_types, by = "water_year") %>%
  # Set factors so graphs plot in logical wetness order rather than alphabetical
  mutate(wy_type = factor(wy_type, levels = c("W", "AN", "BN", "D", "C")))

# Basic data QC
# Remove any erroneous negative sonde turbidity readings
wq_clean <- wq_data %>%
  filter(turb_fnu >= 0 | is.na(turb_fnu))

# Define a turbidity maximum based on RRI historical baseline (e.g., 100 FNU)
qaqc_threshold <- 100

wq_clean <- wq_clean |>
  mutate(
    # If a value exceeds the threshold, replace it with NA, otherwise keep it
    turb_fnu = if_else(turb_fnu > qaqc_threshold, NA_real_, turb_fnu),
    # create a flag column so you have a record of how many points were removed
    qaqc_flag = if_else(is.na(turb_fnu) & !is.na(do_mgl), "Flagged", "Valid")
  )

# Generate completeness matrix by Water Year and Depth
qc_summary <- wq_clean %>%
  group_by(water_year, depth_m) %>%
  summarise(
    total_expected_obs = n(),
    valid_turb_obs = sum(!is.na(turb_fnu)),
    completeness_pct = round((valid_turb_obs / total_expected_obs) * 100, 1),
    mean_turb = mean(turb_fnu, na.rm = TRUE),
    .groups = "drop"
  )

wq_qc_relative <- wq_clean |>
  arrange(depth_m, date_time) |>
  group_by(depth_m, water_year) |>
  mutate(
    # Rolling median using 3 observations before and 3 after (7 total observations)
    rolling_med = slide_dbl(turb_fnu, median, na.rm = TRUE, .before = 3, .after = 3),
    
    # Calculate the absolute deviation from the local median
    median_dev = abs(turb_fnu - rolling_med),
    
    # Flag if the deviation is greater than 20% of the local median 
    # AND the deviation is greater than a 5 FNU noise floor
    flag_spike = if_else(
      median_dev > (0.20 * rolling_med) & median_dev > 5, 
      "Fail: Spike", 
      "Pass",
      missing = "No Data"
    )
  ) |>
  ungroup()

# View the results
wq_qc_relative |>
  count(flag_spike)

wq_qc_relative |>
  count(flag_shift)

# View a summary of the flags to see how much data is being caught
wq_qc_advanced |>
  count(flag_spike, flag_shift)

# 5. SUMMER TURBIDITY ISOLATION & INVASION FLAG
summer_turb <- wq_clean %>%
  # Isolate July, August, September
  filter(month %in% c(7, 8, 9)) 

# 6. VISUALIZATION
ggplot(summer_turb, aes(x = factor(water_year), y = turb_fnu, fill = wy_type)) +
  geom_boxplot(outlier.alpha = 0.3,
               outlier.shape = 21,        # Changes to a fillable circle
               outlier.fill = ) +
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
  p <- ggplot(plot_data, aes(x = date_time, y = turb_fnu)) +
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
         path = here("turbidity-depth-RRI","plots"), 
         width = 10, 
         height = 6)
}
