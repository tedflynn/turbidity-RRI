# RRI summer turbidity by depth: the QC'd three-depth record and its visual screen.
#
# The QC itself is in R/rri_multidepth_qc.R (the general steps in R/qc_functions.R plus the cross-depth
# steps that only apply to RRI's 1, 3 and 6 m sondes). rri_qc_final() runs it, or reads the cached
# result in data/qc/; use rri_qc_final(refresh = TRUE) after changing a data file or a QC setting.
# Output column meanings: qaqc_flag is the QC decision for each reading, turb_fnu_final is turbidity with
# gross-range failures and isolated spikes removed, and turb_fnu_analysis additionally leaves out
# flagged-for-review readings and summer episodes above 100 FNU (analysis_exclusion records why).
library(tidyverse)
library(lubridate)
library(here)
library(slider)

source(here("R", "qc_functions.R"))
source(here("R", "rri_multidepth_qc.R"))
source(here("R", "mussel_analysis.R"))     # prepare_analysis_turbidity(), used inside rri_qc_final()

sj_wy_types <- read_sj_wy_types()
wq_qc_final <- rri_qc_final()
wq_qc_final |> count(qaqc_flag)
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
