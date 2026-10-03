# Is there a shift in turbidity after the Golden Mussel detection (October 2024), using every month?
#
# The summer tests use three months a year. This uses all twelve, as an interrupted time series on
# monthly median turbidity (log10):
#   1. For each station, a model of the baseline months (before the detection) with a calendar-month
#      effect and the station's own flow gives an expected value for every month; the anomaly is the
#      observed minus the expected log10 turbidity.
#   2. Each station's anomalies (or the mean over a group) are then modelled as
#         anomaly ~ linear trend in years + post-detection step
#      with first-order autocorrelation in the errors (nlme::gls, AR(1)), so the step is the change from
#      where the baseline trend would have put the station, not from the baseline average.
# Percent effects are 10^step - 1. 2025 and 2026 are the only post-detection years, so the step is
# estimated from at most 24 months, and the station groups are correlated: read the group lines as one
# piece of evidence, not seven. Same groups as delta-pooled-test.R (chosen after the per-station results
# were seen). Vernalis (SJR) is the control; "minus Vernalis" subtracts its monthly anomaly.
#
# Sources flow-turbidity-network.R for `daily` (daily station turbidity), `gauges` (daily flow covariate x)
# and `stations_tbl` if they are not already in memory.
library(tidyverse)
library(lubridate)
library(here)
library(nlme)

if (!exists("gauges") || !exists("daily")) source(here("flow-turbidity-network.R"))

detect <- ymd("2024-10-01")        # first detection; the step starts with this month
alt_start <- ymd("2025-01-01")     # sensitivity: step starts with the first full year
min_days_per_month <- 20
min_baseline_months <- 36

groups <- list(
  "Delta (all 7)"       = c("PPT", "TWI", "BLP", "ANH", "TRN", "OSJ", "MOK"),
  "Tidal mainstem"      = c("PPT", "TWI", "BLP", "ANH"),
  "Off mainstem"        = c("TRN", "OSJ", "MOK"),
  "Upstream (MSD, RRI)" = c("MSD", "RRI")
)

# Monthly series and anomalies --------------------------------------------------------------------------------------
monthly_turb <- daily |>
  mutate(month_start = floor_date(date, "month")) |>
  group_by(station, month_start) |>
  summarise(days = n(), turb = median(turb), .groups = "drop") |>
  filter(days >= min_days_per_month) |>
  mutate(station = as.character(station))

monthly_flow <- function(g) {
  gauges[[g]] |>
    mutate(month_start = floor_date(date, "month")) |>
    group_by(month_start) |>
    summarise(days = n(), x = mean(x), .groups = "drop") |>
    filter(days >= min_days_per_month) |>
    select(month_start, x)
}

time_cols <- function(d) {
  d |> mutate(t = year(month_start) + (month(month_start) - 0.5) / 12 - 2017,
              idx = (year(month_start) - 2000) * 12 + month(month_start),
              t_detect = year(detect) + (month(detect) - 0.5) / 12 - 2017,
              post = as.numeric(month_start >= detect),
              post_alt = as.numeric(month_start >= alt_start),
              ramp = pmax(0, t - t_detect),
              season = case_when(month(month_start) %in% c(12, 1, 2) ~ "DJF", month(month_start) %in% 3:5 ~ "MAM",
                                 month(month_start) %in% 6:8 ~ "JJA", TRUE ~ "SON"))
}

anoms <- pmap_dfr(stations_tbl, function(station, gauge, role) {
  d <- monthly_turb |>
    filter(station == !!station) |>
    inner_join(monthly_flow(gauge), by = "month_start") |>
    mutate(ly = log10(turb), moy = factor(month(month_start), levels = 1:12))
  base <- filter(d, month_start < detect)
  if (nrow(base) < min_baseline_months) return(NULL)
  fit <- lm(ly ~ moy + x, base)
  d |>
    filter(as.character(moy) %in% as.character(base$moy)) |>
    mutate(anom = ly - predict(fit, newdata = pick(everything())), station = !!station) |>
    select(station, month_start, anom) |>
    time_cols()
})

cat("Months with data per station: baseline (before", format(detect), ") and since:\n")
anoms |> group_by(station) |>
  summarise(baseline_months = sum(post == 0), post_months = sum(post == 1),
            first = format(min(month_start), "%Y-%m"), .groups = "drop") |>
  mutate(station = factor(station, levels = station_order)) |> arrange(station) |> print(n = Inf)

# Interrupted time series ---------------------------------------------------------------------------------------------
its <- function(d, terms, keep = terms[grepl("post|ramp", terms)]) {
  frm <- reformulate(terms, "anom")
  fit <- tryCatch(gls(frm, data = d, correlation = corAR1(form = ~ idx), na.action = na.omit),
                  error = function(e) NULL)
  ar1 <- NA_real_
  if (is.null(fit)) {
    fit <- lm(frm, d); tt <- summary(fit)$coefficients; pcol <- 4
  } else {
    tt <- summary(fit)$tTable; pcol <- 4
    ar1 <- unname(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  }
  map_dfr(keep, function(k) tibble(
    term = k, effect = tt[k, 1], se = tt[k, 2], p = tt[k, pcol], ar1 = ar1,
    trend_per_year = if ("t" %in% terms) tt["t", 1] else NA_real_,
    months = nrow(d), post_months = sum(d$post == 1)))
}
pct <- function(x) round(100 * (10^x - 1))
show_effect <- function(r) r |> mutate(effect_pct = pct(effect), lo_pct = pct(effect - 1.96 * se),
                                       hi_pct = pct(effect + 1.96 * se), p = round(p, 3),
                                       ar1 = round(ar1, 2), trend_pct_per_year = pct(trend_per_year))

# 1. Each station
by_station <- anoms |> group_by(station) |> group_modify(~ its(.x, c("t", "post"))) |> ungroup() |>
  mutate(station = factor(station, levels = station_order)) |> arrange(station)
cat("\nStep after the detection, each station (all months, baseline trend allowed for; % change in turbidity):\n")
by_station |> show_effect() |>
  select(station, effect_pct, lo_pct, hi_pct, p, trend_pct_per_year, ar1, post_months) |> print(n = Inf)

# 2. Groups
group_series <- function(members) {
  anoms |> filter(station %in% members) |> group_by(month_start) |>
    summarise(stations = n(), anom = mean(anom), .groups = "drop") |>
    filter(stations >= ceiling(length(members) / 2)) |> time_cols()
}
sjr <- anoms |> filter(station == "SJR") |> select(month_start, sjr = anom)
group_tbl <- map_dfr(names(groups), function(g) {
  d <- group_series(groups[[g]])
  rbind(mutate(its(d, c("t", "post")), group = g, version = "group"),
        mutate(its(inner_join(d, sjr, by = "month_start") |> mutate(anom = anom - sjr), c("t", "post")),
               group = g, version = "group minus Vernalis"))
})
cat("\nStep after the detection, station groups (mean anomaly across stations):\n")
group_tbl |> show_effect() |>
  select(group, version, effect_pct, lo_pct, hi_pct, p, trend_pct_per_year, ar1, post_months) |>
  arrange(version, group) |> print(n = Inf, width = Inf)

# 3. Is the Delta result sensitive to how the step is defined?
delta <- group_series(groups[["Delta (all 7)"]])
variants <- bind_rows(
  mutate(its(delta, c("t", "post")), variant = "step from Oct 2024, with trend"),
  mutate(its(delta, c("t", "post_alt")), variant = "step from Jan 2025, with trend"),
  mutate(its(delta, "post"), variant = "step from Oct 2024, no trend"),
  mutate(its(delta, c("t", "ramp")), variant = "ramp from Oct 2024, with trend (effect per year)")
)
cat("\nDelta (all 7): other definitions of the change:\n")
variants |> show_effect() |> select(variant, effect_pct, lo_pct, hi_pct, p) |> print(n = Inf, width = Inf)

# The monthly pattern behind it, and how the answer depends on when the change is assumed to start.
# The start-date scan is exploratory: the dates were looked at after seeing the monthly series, and the
# p-values are not adjusted for trying several, so it shows sensitivity, not a test.
cat("\nDelta (all 7): monthly mean anomaly since mid-2024 (stations contributing; blank months lacked data):\n")
delta |> filter(month_start >= ymd("2024-07-01")) |>
  transmute(month = format(month_start, "%Y-%m"), stations, pct_vs_expected = pct(anom)) |> print(n = Inf)

scan <- map_dfr(seq(detect, ymd("2025-07-01"), by = "3 months"), function(s) {
  its(mutate(delta, post_scan = as.numeric(month_start >= s)), c("t", "post_scan")) |>
    mutate(start = format(s, "%Y-%m"), months_after = sum(delta$month_start >= s))
})
cat("\nDelta (all 7): step with trend, by assumed start month (exploratory):\n")
scan |> show_effect() |> select(start, months_after, effect_pct, lo_pct, hi_pct, p) |> print(n = Inf)

# 4. Does the change show up in every season?
delta_s <- delta |> mutate(post_DJF = post * (season == "DJF"), post_MAM = post * (season == "MAM"),
                           post_JJA = post * (season == "JJA"), post_SON = post * (season == "SON"))
cat("\nDelta (all 7): step by season (months in that season since the detection):\n")
its(delta_s, c("t", "post_DJF", "post_MAM", "post_JJA", "post_SON")) |> show_effect() |>
  mutate(season = sub("post_", "", term),
         months_in_season = map_int(season, ~ sum(delta_s$post == 1 & delta_s$season == .x))) |>
  select(season, months_in_season, effect_pct, lo_pct, hi_pct, p) |> print()

# Figure ------------------------------------------------------------------------------------------------------------------
fits <- map_dfr(names(groups), function(g) {
  d <- group_series(groups[[g]])
  fit <- gls(anom ~ t + post, data = d, correlation = corAR1(form = ~ idx), na.action = na.omit)
  d |> mutate(group = g, trend_line = coef(fit)[["(Intercept)"]] + coef(fit)[["t"]] * t)
})
p <- fits |>
  mutate(group = factor(group, levels = names(groups))) |>
  ggplot(aes(month_start, anom)) +
  geom_hline(yintercept = 0, color = "grey70") +
  annotate("rect", xmin = detect, xmax = max(fits$month_start), ymin = -Inf, ymax = Inf, fill = "#d95f0e", alpha = 0.08) +
  geom_line(color = "grey60") +
  geom_point(aes(color = month_start >= detect), size = 1.1) +
  geom_line(aes(y = trend_line), color = "#2c7fb8", linewidth = 0.7) +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#d95f0e"), guide = "none") +
  theme_bw(base_size = 10) +
  labs(title = "Monthly turbidity relative to flow expectation, mean across stations",
       subtitle = "Negative = lower than the baseline seasonal + flow model. Line: fitted trend. Shaded: after Oct 2024",
       x = NULL, y = "Anomaly (log10 turbidity)")
ggsave(here("plots", "monthly_anomaly_its.png"), p, width = 8, height = 9, dpi = 140)
