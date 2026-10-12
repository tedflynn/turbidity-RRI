# Which continuous water quality stations belong in a region-wide turbidity analysis of the Golden Mussel
# introduction? Screening from the station list (data/CWQ-Stations.csv: code, name, lat/long, parameters,
# record start and end), by location, length of pre-detection record and available parameters.
#
# Computed here: distance from the introduction site, years of record before detection, parameter set,
# and distance to the nearest other station (redundancy). Assigned by judgement from the Delta's channel
# layout and flow paths, and recorded in `judgement` below so it can be argued with: the region of each
# station, its tier, and the reason. Nothing in this script knows where mussels are actually established;
# that needs CDFW/DWR detection data (see the note in the output).
#
# Output: data/station_screening.csv and plots/station_selection_map.png.
library(tidyverse)
library(lubridate)
library(here)

detection <- as_date("2024-10-01")        # first detection, Port of Stockton (October 2024)
site <- c(lat = 37.963, lon = -121.365)   # RRI, next to the Port of Stockton: stand-in for the introduction site

stn <- read_csv(here("data", "CWQ-Stations.csv"), show_col_types = FALSE) |>
  transmute(station = Station_Code, name = str_to_title(Station_Name), lat = Latitude, lon = Longitude,
            parameters = Parameters,
            begin = mdy(Dataset_Begin_Date), end = mdy(Dataset_End_Date),
            has_chl = str_detect(Parameters, "Chlorophyll"), has_do = str_detect(Parameters, "DO"),
            has_ph = str_detect(Parameters, "pH"))

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371; p <- pi / 180
  a <- sin((lat2 - lat1) * p / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * r * asin(sqrt(a))
}
stn <- stn |>
  mutate(km_from_site = haversine_km(lat, lon, site["lat"], site["lon"]),
         years_before_detection = round(as.numeric(detection - begin) / 365.25, 1),
         nearest_other_km = map_dbl(seq_len(n()), function(i)
           min(haversine_km(lat[i], lon[i], lat[-i], lon[-i]))))

# Judgement: region, tier and reason -----------------------------------------------------------------------
in_analysis <- c("SJR", "MSD", "RRI", "PPT", "TWI", "BLP", "ANH", "TRN", "OSJ", "MOK")
judgement <- tribble(
  ~station, ~region,                ~tier,        ~reason,
  # south Delta, Old and Middle River: pumped from the SJR, near the export intakes
  "TWA", "South Delta (Old R.)",    "Add: core",   "Longest mid-Old River record (1999); on the pumped pathway from the SJR to the export intakes",
  "MHO", "South Delta (Middle R.)", "Add: core",   "Longest Middle River record (1999), upper end of the river",
  "GLC", "South Delta (Grant Line)", "Add: core",  "Export end of the south Delta by Clifton Court (2007); the place managers care most about",
  "OH1", "South Delta (Old R.)",    "Secondary",   "Head of Old River, 3 km from Mossdale: redundant with MSD, but a record from 2001",
  "MRX", "South Delta (Middle R.)", "Secondary",   "Lower Middle River; same river as MHO, so strongly correlated with it",
  "MRU", "South Delta (Middle R.)", "Secondary",   "Middle River again; redundant with MHO and MRX",
  "ORX", "South Delta (Old R.)",    "Secondary",   "Old River above Doughty Cut; redundant with TWA and OH1",
  "ODM", "South Delta (Old R.)",    "Secondary",   "Below the DMC barrier; redundant with GLC and the intake stations",
  "GLE", "South Delta (Grant Line)", "Secondary",  "Grant Line Canal east; redundant with GLC",
  "ORM", "South Delta (Old R.)",    "Secondary",   "Old River above Mountain House Creek; redundant",
  "VCU", "South Delta (Victoria)",  "Secondary",   "Victoria Canal, west of the pumps; redundant",
  "ORI", "Export intakes",          "Context only", "Clifton Court intake (2018): under 7 years before detection, too few for a baseline, but the best place to watch for a mussel effect at the pumps",
  "WCI", "Export intakes",          "Context only", "West Canal at the Clifton Court intake (2018); same limit as ORI",
  # central Delta
  "OBI", "Central Delta",           "Add: core",   "Old River near Bacon Island (2008), all parameters; interior central Delta between the SJR channel and Franks Tract",
  "HOL", "Central Delta",           "Secondary",   "Holland Cut (2005), long, but no pH or chlorophyll",
  "FRK", "Central Delta (Franks Tract)", "Secondary", "Open, shallow, SAV-dominated (2015): useful for separating SAV from mussels, but only 9 years before detection",
  "BET", "Central Delta",           "Secondary",   "Bethel Island (2014); redundant with FRK and HOL",
  "HLT", "Central Delta (Middle R.)", "Secondary", "Middle River near Holt (2009), no pH or DO",
  "ORQ", "Central Delta",           "Secondary",   "Old River at Quimby Island (2005); turbidity and conductivity only",
  # west Delta, confluence
  "FAL", "West Delta",              "Add: core",   "False River (2005): link between the SJR and Franks Tract; long, with chlorophyll",
  "TSL", "West Delta",              "Secondary",   "Threemile Slough (2008): connects the SJR and Sacramento; no pH or chlorophyll",
  "FCT", "West Delta",              "Secondary",   "Fisherman's Cut (2014); no pH or chlorophyll, redundant with FAL and BLP",
  "SSI", "West Delta",              "Secondary",   "Sherman Island (2015); between BLP and the Sacramento side",
  # Sacramento River: separate source water, candidate controls
  "SRH", "Sacramento River",        "Add: core",   "Sacramento at Hood (2008): upstream, separate watershed and sediment supply; candidate control",
  "RVB", "Sacramento River",        "Add: core",   "Rio Vista (2008): lower Sacramento River; candidate control, tidal",
  "MAL", "Delta gateway",           "Excluded",    "Chipps Island: EDI specific conductance for 2008-2024 is a sorted ramp (about 130 to 17,500 uS/cm) and turbidity noise changes in 2011 and 2025, so the record cannot be QC'd",
  # Suisun Bay: salinity limits freshwater mussels
  "MRZ", "Suisun Bay",              "Secondary",   "Martinez (2007): brackish to saline, so no freshwater mussel effect; an estuary-wide reference with different turbidity drivers",
  "GZL", "Suisun Bay",              "Secondary",   "Grizzly Bay (2015): brackish; reference only",
  "HON", "Suisun Bay",              "Secondary",   "Honker Bay (2015): brackish; reference only",
  "RYC", "Suisun Bay",              "Secondary",   "Ryer / Suisun Bay (2015): brackish; reference only"
) |>
  bind_rows(tibble(station = in_analysis,
                   region = c("SJR corridor", "SJR corridor", "SJR corridor", "Central Delta", "West Delta", "West Delta",
                              "West Delta", "Central Delta", "Central Delta", "Central Delta"),
                   tier = "In analysis now",
                   reason = c("Control, upstream of the introduction", "Mossdale: CDFW reports mussels", "Next to the Port of Stockton",
                              "Prisoners Point", "Twitchell Island", "Blind Point", "Antioch Ship Canal",
                              "Turner Cut", "Old River at Franks Tract near Terminous", "Mokelumne River")))

screen <- stn |>
  left_join(judgement, by = "station") |>
  mutate(tier = factor(tier, levels = c("In analysis now", "Add: core", "Secondary", "Context only", "Excluded"))) |>
  arrange(tier, region, station)

stopifnot(!anyNA(screen$tier))
write_csv(screen, here("data", "station_screening.csv"))

cat("Stations by tier:\n"); print(count(screen, tier))
cat("\nCore additions, with what makes each one useful:\n")
screen |> filter(tier == "Add: core") |>
  select(station, region, begin, years_before_detection, km_from_site, nearest_other_km, parameters) |> print(width = Inf)
cat("\nSecondary and context stations, by region (years of record before detection):\n")
screen |> filter(tier %in% c("Secondary", "Context only")) |>
  select(station, region, tier, years_before_detection, km_from_site, has_chl, has_ph, has_do) |> print(n = Inf, width = Inf)

# Map: stations by tier. No basemap (none installed); the Delta's shape is clear from the stations alone.
tier_col <- c("In analysis now" = "#1B4F72", "Add: core" = "#D55E00", "Secondary" = "#7FB3D5", "Context only" = "#999999", "Excluded" = "#B03A2E")
nudge <- tribble(~station, ~dx, ~dy,       # hand-placed labels where stations crowd together
                 "OH1", 0.020, 0.010, "TWA", 0.020, -0.012, "MRU", -0.026, -0.004, "MHO", -0.026, 0.010, "MRX", -0.026, 0.008,
                 "ODM", 0.020, -0.012, "ORI", 0.020, 0.008, "WCI", 0.020, -0.004, "GLC", 0.020, -0.012, "ORM", 0.022, -0.004,
                 "ORX", -0.026, -0.010, "GLE", 0.020, -0.012, "VCU", 0.020, 0.004, "FRK", -0.026, 0.004, "BET", 0.020, 0.004,
                 "HOL", 0.020, -0.012, "OBI", 0.020, -0.012, "HLT", -0.026, 0.004, "ORQ", -0.026, 0.004, "TRN", 0.020, 0.006,
                 "OSJ", 0.020, 0.008, "PPT", 0.020, -0.014, "MOK", 0.020, 0.004, "FAL", -0.026, 0.004, "TWI", -0.026, 0.004,
                 "TSL", -0.026, 0.004, "FCT", -0.026, -0.010, "BLP", -0.026, -0.010, "SSI", 0.020, 0.008, "ANH", -0.026, 0.004,
                 "RRI", 0.020, 0.008, "MSD", 0.020, 0.004, "SJR", 0.020, 0.004, "RVB", 0.020, 0.004, "MAL", 0.020, 0.008)
lab <- screen |> left_join(nudge, by = "station") |> mutate(dx = replace_na(dx, 0.02), dy = replace_na(dy, 0.006))
p <- ggplot(lab, aes(lon, lat)) +
  geom_point(aes(colour = tier, size = tier == "Add: core"), alpha = 0.95) +
  geom_text(aes(lon + dx, lat + dy, label = station, colour = tier), size = 2.9, hjust = 0.5, show.legend = FALSE) +
  annotate("point", x = site["lon"], y = site["lat"], shape = 1, size = 7, stroke = 1) +
  annotate("text", x = site["lon"] - 0.03, y = site["lat"] - 0.045, label = "introduction site\n(Port of Stockton, approx.)", size = 3, lineheight = 0.9) +
  scale_colour_manual(values = tier_col, name = NULL) +
  scale_size_manual(values = c(`FALSE` = 2.4, `TRUE` = 3.8), guide = "none") +
  coord_fixed(ratio = 1 / cos(38 * pi / 180)) +
  labs(x = "Longitude", y = "Latitude", title = "Continuous water quality stations by role in the region-wide analysis") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
dir.create(here("plots"), showWarnings = FALSE)
ggsave(here("plots", "station_selection_map.png"), p, width = 9, height = 8, dpi = 130)
