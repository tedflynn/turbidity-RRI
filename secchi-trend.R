# Long-term water clarity from the EMP discrete record (monthly Secchi depth, 1975-2025), as context
# for the continuous turbidity results.
#
# 1. Trend by era: is the estuary-wide clearing (higher Secchi) visible, and when did it happen?
# 2. Was summer 2025 (the first summer after the Golden Mussel detection) unusually clear once the
#    era trend and wet years are allowed for? Station by station, then pooled.
#
# Secchi depth rises when the water clears, so this is the opposite sign to turbidity. The EMP file
# ends in December 2025, so 2026 is not covered. Summer = July-September; a station-year needs
# at least two summer samples. Wet years are W and AN in the San Joaquin year type (data/WY-type.csv).
#
# Needs data/DWQ_DWR_EDI_2025.csv. Station roles are assigned from EMP station location (San Joaquin
# corridor vs Sacramento side and Suisun), so they are worth checking against the station list:
#   P8 San Joaquin at Buckley Cove (near RRI), D26 Potato Point, D16 Twitchell Island, D12 Antioch,
#   D19 Franks Tract, D28A Old River at Rancho del Rio.
library(tidyverse)
library(lubridate)
library(here)

first_post_year <- 2025
era_start <- 2000            # a step decline in suspended sediment is documented around 1999
sj_corridor <- c("P8", "D26", "D16", "D12", "D19", "D28A")
context <- c("D4", "D22", "D10", "D6", "D7", "D8")
min_summers <- 15            # summers in the era needed to model a station

wy <- read_csv(here("data", "WY-type.csv"), show_col_types = FALSE, name_repair = "minimal") |>
  select(1, last_col()) |> setNames(c("year", "wy_type")) |>
  mutate(wet = wy_type %in% c("W", "AN"))

secchi <- read_csv(here("data", "DWQ_DWR_EDI_2025.csv"), show_col_types = FALSE,
                   col_types = cols(.default = col_character())) |>
  filter(Analyte == "Secchi Depth", is.na(QC_Flag), Station %in% c(sj_corridor, context)) |>
  transmute(station = Station, date = ymd(Date), cm = suppressWarnings(as.numeric(Result_Value))) |>
  filter(!is.na(cm), cm > 0)

summers <- secchi |>
  filter(month(date) %in% 7:9) |>
  group_by(station, year = year(date)) |>
  summarise(n = n(), cm = median(cm), .groups = "drop") |>
  filter(n >= 2) |>
  mutate(log_cm = log10(cm)) |>
  left_join(wy, by = "year") |>
  mutate(group = if_else(station %in% sj_corridor, "San Joaquin corridor", "Sacramento side / Suisun"))

cat("Summers with data per station (n >= 2 samples), first and last year:\n")
summers |> group_by(station, group) |>
  summarise(summers = n(), first = min(year), last = max(year), median_cm = round(median(cm)), .groups = "drop") |>
  arrange(group, station) |> print(n = Inf)

# 1. Trend by era ----------------------------------------------------------------------------------------------------
era_trend <- function(d, from, to) {
  x <- filter(d, year >= from, year <= to)
  if (nrow(x) < 8) return(tibble(summers = nrow(x), pct_per_year = NA_real_, tau = NA_real_, p = NA_real_))
  fit <- lm(log_cm ~ year, x)
  k <- suppressWarnings(cor.test(x$year, x$log_cm, method = "kendall"))
  tibble(summers = nrow(x), pct_per_year = round(100 * (10^coef(fit)[["year"]] - 1), 1),
         tau = round(unname(k$estimate), 2), p = round(k$p.value, 3))
}
eras <- list("1975-1999" = c(1975, 1999), "2000-2024" = c(2000, 2024), "2010-2024" = c(2010, 2024))
trend_tbl <- summers |>
  group_by(station, group) |>
  group_modify(function(d, key) imap_dfr(eras, ~ mutate(era_trend(d, .x[1], .x[2]), era = .y))) |>
  ungroup()
cat("\nChange in summer Secchi depth by era (% per year; positive = clearer water; Kendall tau and p):\n")
trend_tbl |> select(group, station, era, summers, pct_per_year, tau, p) |>
  pivot_wider(names_from = era, values_from = c(summers, pct_per_year, tau, p), names_vary = "slowest") |>
  arrange(group, station) |> print(n = Inf, width = Inf)

# Decade means of each station's anomaly (log10 Secchi minus its own mean), averaged over stations
cat("\nMean anomaly in log10 Secchi (vs each station's own 1975-2024 mean), by 5-year block and group:\n")
summers |> filter(year <= 2025) |>
  group_by(station) |> mutate(anom = log_cm - mean(log_cm[year <= 2024])) |> ungroup() |>
  mutate(block = cut(year, c(1974, 1979, 1984, 1989, 1994, 1999, 2004, 2009, 2014, 2019, 2024, 2025),
                     labels = c("75-79", "80-84", "85-89", "90-94", "95-99", "00-04", "05-09", "10-14", "15-19", "20-24", "2025"))) |>
  group_by(group, block) |> summarise(anomaly = round(mean(anom), 3), .groups = "drop") |>
  pivot_wider(names_from = block, values_from = anomaly) |> print(width = Inf)

# 2. Summer 2025 against the era model ---------------------------------------------------------------------------------
# log10(Secchi) ~ year + wet, fit to the era's baseline summers (before first_post_year), leaving each year
# out of its own fit. Residuals are standardized against the other baseline years' residuals.
station_resid <- function(d) {
  d <- filter(d, year >= era_start, !is.na(wet))
  base <- filter(d, year < first_post_year)
  if (nrow(base) < min_summers) return(NULL)
  map_dfr(d$year, function(y) {
    train <- filter(base, year != y)
    fit <- lm(log_cm ~ year + wet, train)
    test <- filter(d, year == y)
    tibble(year = y, resid = test$log_cm - predict(fit, test))
  }) |>
    mutate(z = map2_dbl(year, resid, function(y, r) {
      ref <- resid[year != y & year < first_post_year]
      (r - mean(ref)) / sd(ref)
    }))
}
resid_tbl <- summers |>
  group_by(station, group) |>
  group_modify(~ { r <- station_resid(.x); if (is.null(r)) tibble() else r }) |>
  ungroup()

cat("\nSummer", first_post_year, "Secchi against the 2000-2024 trend + wet-year model (z > 0 = clearer than expected):\n")
resid_tbl |>
  filter(year == first_post_year) |>
  left_join(filter(summers, year == first_post_year) |> select(station, cm), by = "station") |>
  mutate(pct_vs_expected = round(100 * (10^resid - 1)), z = round(z, 2), cm = round(cm)) |>
  select(group, station, secchi_cm = cm, pct_vs_expected, z) |>
  arrange(group, station) |> print(n = Inf)

# Pooled over each group; the spread of the baseline years' pooled value includes how correlated stations are
pool_test <- function(g) {
  Tdf <- resid_tbl |> filter(group == g, !is.na(z)) |> group_by(year) |>
    summarise(stations = n(), T = mean(z), .groups = "drop") |> filter(stations >= ceiling(n_distinct(resid_tbl$station[resid_tbl$group == g]) / 2))
  base <- filter(Tdf, year < first_post_year); post <- filter(Tdf, year >= first_post_year)
  if (nrow(post) == 0) return(NULL)
  n <- nrow(base)
  t <- (post$T - mean(base$T)) / (sd(base$T) * sqrt(1 + 1 / n))
  tibble(group = g, post_year = post$year, stations = post$stations, baseline_years = n, T = round(post$T, 2),
         baseline_sd = round(sd(base$T), 2), t = round(t, 2), p_clearer = round(1 - pt(t, n - 1), 3),
         pct_baseline_lower = round(100 * mean(base$T < post$T)))
}
cat("\nPooled across stations (T = mean z; p_clearer = one-sided probability of a value this high or higher):\n")
map_dfr(unique(resid_tbl$group), pool_test) |> print(n = Inf, width = Inf)

# Figure -----------------------------------------------------------------------------------------------------------------
p <- summers |>
  group_by(station) |> mutate(anom = log_cm - mean(log_cm[year <= 2024])) |> ungroup() |>
  group_by(group, year) |> summarise(anomaly = mean(anom), stations = n(), .groups = "drop") |>
  filter(stations >= 3) |>
  ggplot(aes(year, anomaly)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 1999.5, linetype = "dashed", color = "grey60") +
  geom_line(color = "grey55") +
  geom_point(aes(color = year >= first_post_year), size = 1.8) +
  geom_smooth(data = ~ filter(.x, year >= era_start, year < first_post_year), method = "lm", se = FALSE,
              color = "#2c7fb8", linewidth = 0.6) +
  facet_wrap(~ group, ncol = 1) +
  scale_color_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#d95f0e"), guide = "none") +
  theme_bw(base_size = 10) +
  labs(title = "Summer Secchi depth, 1975-2025: mean anomaly across stations",
       subtitle = "log10 Secchi minus each station's own mean; higher = clearer. Dashed line: 1999. Orange: 2025",
       x = NULL, y = "Anomaly in log10 Secchi")
ggsave(here("plots", "secchi_long_term.png"), p, width = 8, height = 6, dpi = 140)
