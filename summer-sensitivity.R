# Is summer turbidity after the Golden Mussel detection unusual compared with earlier summers?
# Compares each post-detection summer (July-September, analysis column) with the baseline of
# earlier summers, with and without the early years 2008-2010 (field methods and equipment
# there are less certain), and with a simple adjustment for wet vs not-wet water-year type.
#
# To add a new summer, put its data in data/ and rerun; years from `first_post_year` onward
# are treated as post-detection and compared with the years before it.
library(tidyverse)
library(lubridate)
library(here)

if (!exists("wq_qc_final")) source(here("turbidity-depth-RRI.R"))

first_post_year <- 2025       # first summer after the detection (assumed Oct 2024; adjust if needed)
early_years <- 2008:2010      # dropped in the second baseline
min_coverage <- 70            # a summer needs at least this % of valid readings to enter a baseline
expected_readings <- 92 * 96  # July-September at 15-minute resolution

# Yearly summer summaries by depth -------------------------------------------------------
yearly <- wq_qc_final |>
  filter(month %in% 7:9) |>
  mutate(year = year(date_time)) |>
  group_by(depth_m, year) |>
  summarise(coverage = round(100 * sum(!is.na(turb_fnu_analysis)) / expected_readings),
            median = median(turb_fnu_analysis, na.rm = TRUE),
            p25 = quantile(turb_fnu_analysis, 0.25, na.rm = TRUE),
            p75 = quantile(turb_fnu_analysis, 0.75, na.rm = TRUE),
            .groups = "drop") |>
  left_join(transmute(sj_wy_types, year = water_year, wy_type = as.character(wy_type)), by = "year") |>
  mutate(wet = wy_type %in% c("W", "AN"), log_median = log(median))

cat("Summer (Jul-Sep) median turbidity (FNU) and coverage (%), analysis column:\n")
yearly |>
  mutate(median = round(median, 2)) |>
  select(year, wy_type, depth_m, coverage, median) |>
  pivot_wider(names_from = depth_m, values_from = c(coverage, median)) |>
  print(n = Inf, width = Inf)

# Compare each post-detection summer with each baseline -------------------------------------------
baselines <- list(
  "all earlier years" = function(y) y < first_post_year,
  "without 2008-2010" = function(y) y < first_post_year & !(y %in% early_years)
)

compare_year <- function(d, post_year, baseline_name) {
  base <- yearly |> filter(depth_m == d, baselines[[baseline_name]](year), coverage >= min_coverage)
  post <- yearly |> filter(depth_m == d, year == post_year)
  nonwet <- base |> filter(!wet)
  trend <- suppressWarnings(cor.test(base$year, base$log_median, method = "kendall"))
  out <- tibble(post_year, depth_m = d, baseline = baseline_name, baseline_years = nrow(base),
                post_median = round(post$median, 2),
                baseline_median = round(median(base$median), 2),
                pct_baseline_lower = round(100 * mean(base$median < post$median)),
                robust_z = round((post$log_median - median(base$log_median)) / (1.4826 * mad(base$log_median)), 2),
                nonwet_years = nrow(nonwet), pct_nonwet_lower = round(100 * mean(nonwet$median < post$median)),
                baseline_trend_tau = round(unname(trend$estimate), 2), baseline_trend_p = round(trend$p.value, 3))
  # Prediction for this summer from a wet / not-wet model on the baseline years (needs the year type)
  if (!is.na(post$wy_type)) {
    pred <- predict(lm(log_median ~ wet, data = base), newdata = data.frame(wet = post$wet), interval = "prediction")
    out <- out |> mutate(predicted = round(exp(pred[1]), 2), pi95_low = round(exp(pred[2]), 2),
                         pi95_high = round(exp(pred[3]), 2),
                         inside_pi = post$median >= exp(pred[2]) & post$median <= exp(pred[3]))
  }
  out
}

post_years <- yearly |> filter(year >= first_post_year) |> pull(year) |> unique()
results <- expand_grid(post_year = post_years, d = c(1, 3, 6), baseline_name = names(baselines)) |>
  pmap_dfr(~ compare_year(..2, ..1, ..3))

cat("\nPost-detection summer vs baseline summers:\n")
results |>
  select(post_year, depth_m, baseline, baseline_years, post_median, baseline_median, pct_baseline_lower, robust_z) |>
  print(n = Inf, width = Inf)

cat("\nAdjusted for wet vs not-wet year type (prediction interval from the baseline years):\n")
results |>
  select(post_year, depth_m, baseline, nonwet_years, pct_nonwet_lower, any_of(c("predicted", "pi95_low", "pi95_high", "inside_pi"))) |>
  print(n = Inf, width = Inf)

cat("\nTrend within the baseline years (Kendall tau of log summer median vs year):\n")
results |>
  distinct(depth_m, baseline, baseline_trend_tau, baseline_trend_p) |>
  print(n = Inf)

# Plot -------------------------------------------------------------------------------------------------
p <- yearly |>
  mutate(depth = factor(paste0(depth_m, " m"), levels = c("1 m", "3 m", "6 m")),
         group = case_when(year >= first_post_year ~ "post-detection summer",
                           year %in% early_years ~ "2008-2010",
                           TRUE ~ "other earlier years")) |>
  ggplot(aes(year, median, color = group)) +
  annotate("rect", xmin = min(early_years) - 0.5, xmax = max(early_years) + 0.5, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_linerange(aes(ymin = p25, ymax = p75), linewidth = 0.6) +
  geom_point(aes(shape = wy_type), size = 2.4) +
  facet_wrap(~ depth, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("2008-2010" = "#7f7f7f", "other earlier years" = "#2c7fb8", "post-detection summer" = "#d95f0e"), name = NULL) +
  scale_shape_manual(values = c(W = 17, AN = 15, BN = 16, D = 1, C = 4), name = "SJR year type") +
  scale_x_continuous(breaks = seq(2008, max(yearly$year), 2)) +
  theme_bw(base_size = 10) +
  labs(title = "Summer (Jul-Sep) turbidity: yearly median and interquartile range", x = NULL, y = "Turbidity (FNU)")
ggsave(here("plots", "summer_median_by_year.png"), p, width = 9, height = 9, dpi = 140)
