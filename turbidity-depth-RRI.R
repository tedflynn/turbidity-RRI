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

# 5. SUMMER TURBIDITY ISOLATION & INVASION FLAG
summer_turb <- wq_clean %>%
  # Isolate July, August, September
  filter(month %in% c(7, 8, 9)) 

# 6. VISUALIZATION
ggplot(summer_turb, aes(x = factor(water_year), y = turb_fnu, fill = wy_type)) +
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
         path = here("plots"), 
         width = 10, 
         height = 6)
}
