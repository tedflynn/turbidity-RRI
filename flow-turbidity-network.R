# How stable is the turbidity-flow relationship through time and along the San Joaquin River,
# and is 2025 unusual at any station once flow is accounted for?
#
# Ten stations with continuous turbidity, each paired with its nearest CDEC flow gauge (table
# below). Each station's turbidity comes from the stored QC (run-qc.R; RRI uses the three-depth
# pipeline in R/rri_multidepth_qc.R); daily medians are compared with daily flow.
#
# Flow gauges (data/CDEC_*): VNS and MSD are daily means. SJJ, PRI, TRN, OSJ and MOK are
# tidal 15-minute flow, Godin-filtered to a daily net flow (see R/flow_functions.R). The
# RRI gauge itself ends in 2018, so RRI turbidity is paired with Vernalis flow instead. That
# was chosen from a check (printed below) of how well Vernalis and Brandt Bridge (BDT) track
# the RRI gauge: Vernalis does better, and BDT is tidal and often at or below zero in summer.
#
# Covariate x: log10(flow) for the daily river gauges (positive-only, so days at or below zero
# are dropped); flow in 1000 cfs for the tidally filtered gauges, whose net flow in summer is
# near zero and often negative.
#
# Provisional 2026 (use_cdec_2026): until the full 2026 files arrive, turbidity for 2026 comes from
# CDEC downloads (data/CDEC_<station>_turb_2026.csv), which are real-time and not QC'd, are turbidity
# only, and are in NTU at SJR, PPT, MSD, RRI, TWI and ANH but FNU elsewhere (no conversion applied).
# The turbidity-only QC was checked on 2024-2025: leaving out conductivity and DO changed no station's
# summer median and excluded at most 0.1% more readings. CDEC has a single (1 m) sonde at RRI, so
# with this switch on RRI is the 1 m series throughout, not the median across three depths.
#
# Vernalis (SJR) is the only station known to be free of mussels; Mossdale and everything
# downstream may not be. To add 2026, extend the data files and rerun; years from
# `first_post_year` onward are treated as post-detection.
library(tidyverse)
library(lubridate)
library(slider)
library(here)

source(here("R", "qc_functions.R"))
source(here("R", "station_qc.R"))
source(here("R", "rri_multidepth_qc.R"))
source(here("R", "mussel_analysis.R"))
source(here("R", "flow_functions.R"))

first_post_year <- 2025
use_cdec_2026 <- TRUE          # set FALSE once the full 2026 data files are in place
early_years <- 2008:2010
min_valid_readings_per_day <- 48
min_days_per_year <- 200      # a calendar year needs this many days to give a correlation
min_days_per_summer <- 60

# Stations, roughly upstream to downstream along the river, then the off-mainstem sites
stations_tbl <- tribble(
  ~station, ~gauge, ~role,
  "SJR", "VNS", "control, upstream",
  "MSD", "MSD", "mainstem",
  "RRI", "VNS", "mainstem",
  "PPT", "PRI", "mainstem",
  "TWI", "SJJ", "mainstem, tidal",
  "BLP", "SJJ", "mainstem, tidal",
  "ANH", "SJJ", "mainstem, tidal",
  "TRN", "TRN", "off mainstem",
  "OSJ", "OSJ", "off mainstem",
  "MOK", "MOK", "off mainstem"
)
station_order <- stations_tbl$station

# Turbidity ----------------------------------------------------------------------------------------
# The QC is separate: run-qc.R stores the QC'd single-sonde stations in data/qc/ and rri_qc_final() the RRI
# three-depth record (the one special case). build_daily_turbidity() (R/mussel_analysis.R) leaves out the
# flagged readings and takes daily medians. Rerun run-qc.R first if a data file or QC setting changed.
if (!exists("daily")) {
  daily <- build_daily_turbidity(station_order, use_cdec_2026, min_valid_readings_per_day)
  cat("QC outcome for the single-sonde stations:\n")
  print(attr(daily, "qc_summary"))
}

# Flow ------------------------------------------------------------------------------------------------
daily_gauge <- function(id) read_cdec_daily(here("data", paste0("CDEC_", id, "_flow_daily.csv")))
tidal_gauge <- function(id) godin_daily_flow(here("data", paste0("CDEC_", id, "_flow_15min.csv")))
as_log <- function(x) x |> filter(flow > 0) |> transmute(date, x = log10(flow))
as_net <- function(x) x |> transmute(date, x = flow / 1000)

raw_flow <- list(VNS = daily_gauge("VNS"), MSD = daily_gauge("MSD"), BDT = daily_gauge("BDT"),
                 SJJ = tidal_gauge("SJJ"), PRI = tidal_gauge("PRI"), TRN = tidal_gauge("TRN"),
                 OSJ = tidal_gauge("OSJ"), MOK = tidal_gauge("MOK"))
gauges <- c(map(raw_flow[c("VNS", "MSD")], as_log), map(raw_flow[c("SJJ", "PRI", "TRN", "OSJ", "MOK")], as_net))

cat("\nFlow gauges: days with usable flow, record span, and summer days at or below zero:\n")
imap_dfr(raw_flow, ~ tibble(gauge = .y, days = nrow(.x), first = min(.x$date), last = max(.x$date),
                            summer_days_le_0 = sum(month(.x$date) %in% 7:9 & .x$flow <= 0))) |> print()

cat("\nDays with usable flow in July-September, last three years:\n")
imap_dfr(raw_flow, ~ .x |> filter(month(date) %in% 7:9, year(date) >= 2024) |> count(year = year(date)) |> mutate(gauge = .y)) |>
  pivot_wider(names_from = year, values_from = n) |> print()

# How well does Brandt Bridge (BDT) stand in for the RRI gauge, and how does Vernalis compare? ---------
rri_gauge <- tidal_gauge("RRI")
cat("\nRRI gauge (to 2018) against BDT and Vernalis, Spearman correlation of daily flow:\n")
standin <- inner_join(rename(rri_gauge, rri = flow), rename(raw_flow$BDT, bdt = flow), by = "date") |>
  inner_join(rename(raw_flow$VNS, vns = flow), by = "date")
standin |>
  mutate(period = if_else(month(date) %in% 7:9, "Jul-Sep", "all days")) |>
  bind_rows(mutate(standin, period = "all days")) |>
  group_by(period) |>
  summarise(days = n(), rri_vs_bdt = cor(rri, bdt, method = "spearman"),
            rri_vs_vernalis = cor(rri, vns, method = "spearman"), bdt_vs_vernalis = cor(bdt, vns, method = "spearman")) |>
  mutate(across(where(is.numeric) & !days, ~ round(.x, 2))) |>
  distinct() |>
  print()

# Summer level by station and year ------------------------------------------------------------------------------------
summer_medians <- daily |>
  filter(month(date) %in% 7:9) |>
  group_by(station, year = year(date)) |>
  summarise(days = n(), median = median(turb), .groups = "drop") |>
  filter(days >= min_days_per_summer)
cat("\nSummer (Jul-Sep) median turbidity by station and year (FNU/NTU):\n")
summer_medians |>
  select(-days) |>
  mutate(median = round(median, 1)) |>
  pivot_wider(names_from = station, values_from = median) |>
  select(year, any_of(station_order)) |>
  arrange(year) |>
  print(n = Inf)

# 1. Turbidity-flow correlation by year -----------------------------------------------------------------------
joined <- function(s, gauge) {
  daily |> filter(station == s) |> inner_join(gauges[[gauge]], by = "date")
}
yearly_cor <- function(s, gauge) {
  x <- joined(s, gauge) |> mutate(year = year(date), summer = month(date) %in% 7:9)
  bind_rows(
    x |> group_by(year) |> filter(n() >= min_days_per_year) |>
      summarise(rho = cor(turb, x, method = "spearman"), days = n(), season = "all days", .groups = "drop"),
    x |> filter(summer) |> group_by(year) |> filter(n() >= min_days_per_summer) |>
      summarise(rho = cor(turb, x, method = "spearman"), days = n(), season = "Jul-Sep", .groups = "drop")
  ) |> mutate(station = s, gauge = gauge)
}
by_year <- bind_rows(
  pmap_dfr(stations_tbl, function(station, gauge, role) mutate(yearly_cor(station, gauge), driver = "own gauge")),
  map_dfr(setdiff(station_order, "SJR"), function(s) mutate(yearly_cor(s, "VNS"), driver = "Vernalis flow"))
) |> mutate(station = factor(station, levels = station_order))

show_by_year <- function(season_, driver_) {
  by_year |>
    filter(season == season_, driver == driver_) |>
    mutate(rho = round(rho, 2)) |>
    select(year, station, rho) |>
    pivot_wider(names_from = station, values_from = rho) |>
    select(year, any_of(station_order)) |>
    arrange(year) |>
    print(n = Inf)
}
cat("\nSpearman correlation of daily turbidity with the station's own flow, by year, all days:\n")
show_by_year("all days", "own gauge")
cat("\n... July-September only:\n")
show_by_year("Jul-Sep", "own gauge")
cat("\nSpearman correlation with Vernalis flow (same day), July-September:\n")
show_by_year("Jul-Sep", "Vernalis flow")

stability <- by_year |>
  group_by(station, driver, season) |>
  summarise(
    years = sum(year < first_post_year),
    median_rho = median(rho[year < first_post_year]),
    min_rho = min(rho[year < first_post_year]),
    max_rho = max(rho[year < first_post_year]),
    tau_vs_year = suppressWarnings(cor(year[year < first_post_year], rho[year < first_post_year], method = "kendall")),
    post_rho = list(setNames(round(rho[year >= first_post_year], 2), year[year >= first_post_year])),
    pct_baseline_lower = list(map_dbl(rho[year >= first_post_year], ~ round(100 * mean(rho[year < first_post_year] < .x)))),
    .groups = "drop")
cat("\nHow stable is the correlation? Baseline years (before the first post-detection year):\n")
stability |>
  filter(driver == "own gauge") |>
  mutate(across(c(median_rho, min_rho, max_rho, tau_vs_year), ~ round(.x, 2)),
         post_rho = map_chr(post_rho, ~ paste(.x, collapse = "/")),
         pct_baseline_lower = map_chr(pct_baseline_lower, ~ paste(.x, collapse = "/"))) |>
  arrange(season, station) |>
  print(n = Inf, width = Inf)

# 2. Does 2025 (or a later summer) depart from what flow predicts? ---------------------------------------------
# Leave-one-year-out on July-September days: each baseline year is predicted from a model fit to
# the other baseline years, so baseline residuals are as honest as the post-detection ones.
#   flow only:        log10(turbidity) ~ own flow
#   flow + Vernalis:  log10(turbidity) ~ own flow + log10(Vernalis turbidity, 3 days earlier) + Vernalis flow
vern_turb <- daily |> filter(station == "SJR") |> transmute(date = date + 3, vernalis_turb = turb)
fit_data <- function(s, gauge, with_vernalis) {
  d <- joined(s, gauge) |> filter(month(date) %in% 7:9)
  if (with_vernalis) {
    d <- d |>
      inner_join(vern_turb, by = "date") |>
      inner_join(rename(gauges$VNS, vern_x = x), by = "date") |>
      mutate(lv = log10(vernalis_turb))
  }
  d |> mutate(year = year(date), ly = log10(turb)) |> filter(is.finite(ly))
}

loyo <- function(s, gauge, with_vernalis, baseline_years, post_year) {
  d <- fit_data(s, gauge, with_vernalis)
  frm <- if (with_vernalis) ly ~ x + lv + vern_x else ly ~ x
  base_years <- intersect(baseline_years, unique(d$year))
  base_resid <- map_dfr(base_years, function(y) {
    test <- filter(d, year == y)
    if (nrow(test) < 40 || n_distinct(setdiff(base_years, y)) < 4) return(NULL)
    fit <- lm(frm, filter(d, year %in% setdiff(base_years, y)))
    tibble(year = y, mean_resid = mean(test$ly - predict(fit, test)))
  })
  test <- filter(d, year == post_year)
  if (nrow(test) < 40 || nrow(base_resid) < 5) return(NULL)
  fit <- lm(frm, filter(d, year %in% base_years))
  post_resid <- mean(test$ly - predict(fit, test))
  tibble(station = s, model = if (with_vernalis) "flow + Vernalis" else "flow only", post_year,
         baseline_years = nrow(base_resid), sd_resid_log10 = round(sd(base_resid$mean_resid), 3),
         pct_vs_expected = round(100 * (10^post_resid - 1)),
         pct_baseline_lower = round(100 * mean(base_resid$mean_resid < post_resid)),
         z = round((post_resid - mean(base_resid$mean_resid)) / sd(base_resid$mean_resid), 2),
         model_r2 = round(summary(fit)$r.squared, 2))
}

post_years <- sort(unique(year(daily$date)[year(daily$date) >= first_post_year]))
all_years <- sort(unique(year(daily$date)))
baselines <- list("all earlier years" = setdiff(all_years, post_years),
                  "without 2008-2010" = setdiff(all_years, c(post_years, early_years)))
prediction <- expand_grid(stations_tbl, vern = c(FALSE, TRUE), baseline = names(baselines), post_year = post_years) |>
  filter(!(station == "SJR" & vern)) |>
  pmap_dfr(function(station, gauge, role, vern, baseline, post_year) {
    r <- loyo(station, gauge, vern, baselines[[baseline]], post_year)
    if (is.null(r)) NULL else mutate(r, baseline = baseline, role = role)
  }) |>
  mutate(station = factor(station, levels = station_order)) |>
  arrange(model, baseline, post_year, station)

cat("\nPost-detection summer vs expectation from flow (z = standard deviations from the baseline years;\n",
    "pct_baseline_lower = share of baseline years with a lower residual):\n", sep = "")
print(prediction, n = Inf, width = Inf)

# 3. Figure -----------------------------------------------------------------------------------------------------
p <- by_year |>
  filter(driver == "own gauge") |>
  ggplot(aes(year, rho, color = season, group = season)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_line(alpha = 0.5, na.rm = TRUE) +
  geom_point(aes(shape = year >= first_post_year), size = 1.8) +
  facet_wrap(~ station, ncol = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_color_manual(values = c("all days" = "#2c7fb8", "Jul-Sep" = "#d95f0e"), name = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(title = "Turbidity-flow correlation by year, each station with its nearest flow gauge",
       x = NULL, y = "Spearman correlation (daily turbidity, daily flow)")
ggsave(here("plots", "turbidity_flow_by_year.png"), p, width = 8, height = 10, dpi = 140)
