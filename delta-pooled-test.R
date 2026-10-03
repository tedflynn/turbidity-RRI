# Is summer turbidity in the Delta stations lower than flow predicts, taken together?
#
# flow-turbidity-network.R tests each station separately, but the stations share water, weather and
# operations, so ten separate z-scores overcount the evidence. This pools them.
#
# For each station and summer (Jul-Sep), the residual is the mean of log10(turbidity) minus the
# prediction from the flow model fit to the other baseline years (leave-one-year-out, as in
# flow-turbidity-network.R). Each station's residuals are standardized against its baseline years
# (also leaving the year out), and a group's statistic T for a year is the mean of its stations' z.
# The post-detection T is compared with the baseline years' T. Because that spread comes from real
# years, it already contains how strongly the stations move together; the mean pairwise correlation
# and the implied effective number of independent stations are printed too.
#
# Groups were fixed from station type, but after the per-station results had been seen, so treat the
# choice as post hoc:
#   Delta (all 7)     PPT, TWI, BLP, ANH, TRN, OSJ, MOK
#   Tidal mainstem    PPT, TWI, BLP, ANH
#   Off mainstem      TRN, OSJ, MOK
#   Upstream (others) MSD, RRI   (Vernalis, SJR, is the control and is kept separate)
# "Minus Vernalis" subtracts Vernalis's own z for the same year, which removes anything regional that
# also moved the control. The 2026 turbidity is provisional (see provisional-vs-qc-2025.R); an
# adjustment line rescales it by each station's 2025 CDEC / QC'd summer ratio.
#
# Sources flow-turbidity-network.R if its objects are not in memory, and uses `provisional_ratio`
# from provisional-vs-qc-2025.R if present.
library(tidyverse)
library(lubridate)
library(slider)
library(here)

if (!exists("prediction")) source(here("flow-turbidity-network.R"))

groups <- list(
  "Delta (all 7)"     = c("PPT", "TWI", "BLP", "ANH", "TRN", "OSJ", "MOK"),
  "Tidal mainstem"    = c("PPT", "TWI", "BLP", "ANH"),
  "Off mainstem"      = c("TRN", "OSJ", "MOK"),
  "Upstream (MSD, RRI)" = c("MSD", "RRI")
)

# Year-level residuals -------------------------------------------------------------------------------------------
year_resid <- function(s, gauge, with_vernalis, min_days = 40) {
  d <- fit_data(s, gauge, with_vernalis)
  frm <- if (with_vernalis) ly ~ x + lv + vern_x else ly ~ x
  years <- d |> count(year) |> filter(n >= min_days) |> pull(year)
  base <- years[years < first_post_year]
  map_dfr(years, function(y) {
    train <- setdiff(base, y)
    if (length(train) < 4) return(NULL)
    fit <- lm(frm, filter(d, year %in% train))
    test <- filter(d, year == y)
    tibble(station = s, year = y, resid = mean(test$ly - predict(fit, test)))
  })
}

resid_all <- bind_rows(
  pmap_dfr(stations_tbl, function(station, gauge, role) mutate(year_resid(station, gauge, FALSE), model = "flow only")),
  pmap_dfr(filter(stations_tbl, station != "SJR"), function(station, gauge, role) mutate(year_resid(station, gauge, TRUE), model = "flow + Vernalis"))
)

# Optional adjustment for the provisional 2026 data: if CDEC read lower than the QC'd record by `ratio`
# in 2025, the 2026 numbers are scaled up by 1 / ratio (the same shift applied to every 2026 reading).
if (exists("provisional_ratio")) {
  resid_adjusted <- resid_all |>
    left_join(provisional_ratio, by = "station") |>
    mutate(resid = if_else(year >= 2026 & !is.na(ratio), resid - log10(ratio), resid)) |>
    select(-ratio)
}

# Standardize against the baseline years, leaving the year itself out of its own reference -----------------
standardize <- function(res) {
  res |>
    group_by(model, station) |>
    group_modify(function(g, key) {
      base <- filter(g, year < first_post_year)
      g |> mutate(z = map2_dbl(year, resid, function(y, r) {
        ref <- base$resid[base$year != y]
        if (length(ref) < 5) NA_real_ else (r - mean(ref)) / sd(ref)
      }))
    }) |>
    ungroup()
}

# Group statistic and its test against the baseline years ----------------------------------------------------------
group_T <- function(z, members) {
  need <- ceiling(length(members) / 2)
  z |> filter(station %in% members, !is.na(z)) |>
    group_by(year) |> summarise(stations = n(), T = mean(z), .groups = "drop") |>
    filter(stations >= need)
}
test_T <- function(Tdf) {
  base <- filter(Tdf, year < first_post_year)
  post <- filter(Tdf, year >= first_post_year)
  if (nrow(base) < 6 || nrow(post) == 0) return(NULL)
  n <- nrow(base)
  pmap_dfr(post, function(year, stations, T) {
    t <- (T - mean(base$T)) / (sd(base$T) * sqrt(1 + 1 / n))
    tibble(post_year = year, stations = stations, baseline_years = n, T = round(T, 2),
           baseline_mean = round(mean(base$T), 2), baseline_sd = round(sd(base$T), 2), t = round(t, 2),
           p_lower = round(pt(t, n - 1), 3), pct_baseline_lower = round(100 * mean(base$T < T)))
  })
}
mean_pairwise_cor <- function(z, members) {
  w <- z |> filter(station %in% members, year < first_post_year) |>
    select(year, station, z) |> pivot_wider(names_from = station, values_from = z) |> select(-year)
  r <- suppressWarnings(cor(w, use = "pairwise.complete.obs"))
  mean(r[upper.tri(r)], na.rm = TRUE)
}

run_groups <- function(z, label) {
  map_dfr(names(groups), function(g) {
    members <- groups[[g]]
    res <- test_T(group_T(z, members))
    if (is.null(res)) return(NULL)
    rbar <- mean_pairwise_cor(z, members)
    n_st <- length(members)
    mutate(res, group = g, set = label, mean_pairwise_r = round(rbar, 2),
           effective_stations = round(n_st / (1 + (n_st - 1) * rbar), 1))
  })
}

z_flow <- standardize(filter(resid_all, model == "flow only"))
z_vern <- standardize(filter(resid_all, model == "flow + Vernalis"))

show <- function(x) x |> select(set, group, post_year, stations, T, baseline_mean, baseline_sd, t, p_lower, pct_baseline_lower,
                                 mean_pairwise_r, effective_stations) |> arrange(set, post_year, group)

cat("Pooled z (negative = lower than flow predicts), post-detection summers vs baseline years.\n",
    "T = mean z across the group's stations; p_lower = one-sided probability of a value this low or lower\n",
    "from the baseline years' T distribution (t distribution, accounts for estimating its spread):\n\n", sep = "")
pooled <- bind_rows(run_groups(z_flow, "flow only"), run_groups(z_vern, "flow + Vernalis"))
print(show(pooled), n = Inf, width = Inf)

# Minus Vernalis (difference from the control) ------------------------------------------------------------------------
z_sjr <- z_flow |> filter(station == "SJR") |> select(year, z_sjr = z)
did <- map_dfr(names(groups), function(g) {
  Tdf <- group_T(z_flow, groups[[g]]) |> inner_join(z_sjr, by = "year") |> mutate(T = T - z_sjr)
  res <- test_T(select(Tdf, year, stations, T))
  if (is.null(res)) NULL else mutate(res, group = g, set = "flow only, minus Vernalis")
})
cat("\nSame, as a difference from Vernalis (the control) in each year:\n")
did |>
  select(set, group, post_year, stations, T, baseline_mean, baseline_sd, t, p_lower, pct_baseline_lower) |>
  arrange(post_year, group) |>
  print(n = Inf, width = Inf)

# Allowing for the estuary-wide decline in turbidity and for wet years ------------------------------------------------
# Turbidity across the estuary has been falling for decades, so a plain "is this year below the baseline
# average" test can flag recent years even with no new cause. This fits the group's baseline T against
# year (a linear trend) and against wet / not-wet year type, then asks whether 2025 and 2026 fall below
# what that model expects (prediction interval, so the model's own uncertainty is included).
# WY-type.csv ends in 2025; 2026 is entered below as reported (Above Normal) and should be checked,
# including which index (San Joaquin or Sacramento) the classification refers to.
wy_type_2026 <- "AN"
wy <- read_csv(here("data", "WY-type.csv"), show_col_types = FALSE, name_repair = "minimal") |>
  select(1, last_col()) |> setNames(c("year", "wy_type")) |>
  filter(year != 2026) |>
  bind_rows(tibble(year = 2026, wy_type = wy_type_2026))

test_T_adjusted <- function(Tdf, trend, wet) {
  d <- Tdf |> left_join(wy, by = "year") |> mutate(wet = wy_type %in% c("W", "AN"))
  base <- filter(d, year < first_post_year)
  post <- filter(d, year >= first_post_year)
  terms <- c(if (trend) "year", if (wet) "wet")
  fit <- lm(reformulate(if (length(terms)) terms else "1", "T"), base)
  pr <- predict(fit, post, se.fit = TRUE)
  t <- (post$T - pr$fit) / sqrt(pr$se.fit^2 + sigma(fit)^2)
  tibble(post_year = post$year, wy_type = post$wy_type, T = round(post$T, 2), expected = round(pr$fit, 2),
         t = round(t, 2), p_lower = round(pt(t, df.residual(fit)), 3),
         trend_per_year = if (trend) round(coef(fit)[["year"]], 3) else NA_real_,
         trend_p = if (trend) round(summary(fit)$coefficients["year", 4], 3) else NA_real_,
         wet_effect = if (wet) round(coef(fit)[["wetTRUE"]], 2) else NA_real_)
}
adjusted <- expand_grid(group = names(groups),
                        model = c("mean only", "year trend", "wet year", "year trend + wet year")) |>
  pmap_dfr(function(group, model) {
    Tdf <- group_T(z_flow, groups[[group]])
    res <- test_T_adjusted(Tdf, trend = grepl("trend", model), wet = grepl("wet", model))
    mutate(res, group = group, model = model)
  })
cat("\nPooled T against what the baseline trend and year type expect (flow-only residuals):\n")
adjusted |>
  select(group, model, post_year, wy_type, T, expected, t, p_lower, trend_per_year, trend_p, wet_effect) |>
  arrange(group, post_year, model) |>
  print(n = Inf, width = Inf)

# Dropping one station at a time from the Delta group ---------------------------------------------------------------
cat("\nDelta group with one station left out (flow only):\n")
loo <- map_dfr(groups[["Delta (all 7)"]], function(drop) {
  res <- test_T(group_T(z_flow, setdiff(groups[["Delta (all 7)"]], drop)))
  if (is.null(res)) NULL else mutate(res, left_out = drop)
})
loo |> select(left_out, post_year, T, baseline_mean, baseline_sd, t, p_lower, pct_baseline_lower) |>
  arrange(post_year, left_out) |> print(n = Inf, width = Inf)

# Provisional-data adjustment ----------------------------------------------------------------------------------------
if (exists("resid_adjusted")) {
  cat("\n2026 rescaled by each station's 2025 CDEC / QC'd summer ratio (flow only); 2025 is unchanged:\n")
  print(provisional_ratio |> mutate(ratio = round(ratio, 3)), n = Inf)
  adj <- run_groups(standardize(filter(resid_adjusted, model == "flow only")), "flow only, 2026 rescaled")
  adj |> filter(post_year == 2026) |> show() |> print(n = Inf, width = Inf)
}

# Figure ----------------------------------------------------------------------------------------------------------------
p <- map_dfr(names(groups), function(g) mutate(group_T(z_flow, groups[[g]]), group = g)) |>
  mutate(group = factor(group, levels = names(groups)), post = year >= first_post_year) |>
  ggplot(aes(year, T)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_line(color = "grey55") +
  geom_point(aes(color = post), size = 2) +
  facet_wrap(~ group, ncol = 2) +
  scale_color_manual(values = c(`FALSE` = "#2c7fb8", `TRUE` = "#d95f0e"), guide = "none") +
  theme_bw(base_size = 10) +
  labs(title = "Summer turbidity relative to flow expectation, pooled across stations",
       subtitle = "Mean standardized residual (negative = lower than expected); orange = post-detection summers",
       x = NULL, y = "Pooled z")
ggsave(here("plots", "delta_pooled_by_year.png"), p, width = 8, height = 6, dpi = 140)
