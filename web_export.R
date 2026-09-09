##########################################
## Web export — JSON payloads for kristofmoeller.com hero widgets
## Run AFTER Time_series_analysis.R has populated all_data.txt and filtered_data.txt.
## Produces web_export/{stations.json, yearly_probability.json, meta.json}.
##########################################

# 1. Setup -----------------------------------------------------------------
script_dir <- here::here()
source(file.path(script_dir, "Time_series_analysis_custom_functions.R"))
install_packages()

if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(jsonlite)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

out_dir <- file.path(script_dir, "web_export")
dir.create(out_dir, showWarnings = FALSE)

# 2. Load the two canonical outputs of the main pipeline -------------------
all_data <- read_delim(
  file.path(script_dir, "all_data.txt"),
  delim = "\t", col_names = TRUE, show_col_types = FALSE
)

filtered_data <- read_delim(
  file.path(script_dir, "filtered_data.txt"),
  delim = "\t", col_names = TRUE, show_col_types = FALSE
)

# Re-derive the station_map exactly as in Time_series_analysis.R so we can
# attach short station codes (B1..B9, D1..D18, S1..S14, NW1..NW4, L1, L2, N1, N2, S12).
# Keep this block in sync with the master script if it ever changes.
station_map <- c(
  "Kjempebakken"                   = "NW4",
  "Møkland"                        = "NW3",
  "RA2"                            = "B1",
  "A13"                            = "B2",
  "B7"                             = "B3",
  "Korsfjorden"                    = "NW2",
  "GA1"                            = "B4",
  "C3"                             = "B5",
  "Bjørnafjorden"                  = "NW1",
  "Indre Oslofjord"                = "S12",
  "Håøyfjorden"                    = "S13",
  "Kosterfjorden (NR16)"           = "S11",
  "SLV Bottnefjorden"              = "S10",
  "Arendal"                        = "S14",
  "Stretudden"                     = "S9",
  "Havstensfjord"                  = "S8",
  "SLV Saltöfjorden"               = "S7",
  "Å17"                            = "S6",
  "SLÄGGÖ"                         = "S4",
  "SLV Havstensfjorden-Ljungskile" = "S5",
  "Koljöfjord"                     = "S3",
  "SLV Lyresund-Stigfjorden"       = "S2",
  "Åstol"                          = "S1",
  "DANAFJORD"                      = "D18",
  "BY15 GOTLANDSDJ"                = "B6",
  "N7 OST Nidingen"                = "D17",
  "VIB3708"                        = "L2",
  "N14 Falkenberg"                 = "D15",
  "NOR409"                         = "D16",
  "ANHOLT E"                       = "D13",
  "NOR5503"                        = "D14",
  "VIB3727"                        = "L1",
  "L9  LAHOLMSBUKTEN"              = "D12",
  "9E REF M1V1"                    = "B7",
  "ARH170006"                      = "D10",
  "VSJ20925"                       = "D11",
  "RKB1"                           = "N2",
  "KBH431"                         = "D6",
  "ROS60"                          = "D7",
  "VEJ0006870"                     = "D9",
  "FYN6900017"                     = "D8",
  "RIB1510007"                     = "N1",
  "BRKBMPK2"                       = "B8",
  "FYN6300043"                     = "D5",
  "BY2 ARKONA"                     = "B9",
  "TF0360"                         = "D4",
  "TF0046"                         = "D1",
  "MECKLENBURGER BUCHT"            = "D3",
  "Heiligendamm"                   = "D2"
)

# Helper: classify a short station code into its sub-basin group.
basin_from_code <- function(code) {
  case_when(
    is.na(code)                              ~ NA_character_,
    startsWith(code, "B")                    ~ "Baltic Sea",
    startsWith(code, "D")                    ~ "Danish Straits",
    startsWith(code, "NW")                   ~ "Norwegian Sea",
    startsWith(code, "N") & !startsWith(code, "NW") ~ "North Sea",
    startsWith(code, "S")                    ~ "Skagerrak",
    startsWith(code, "L")                    ~ "Limfjord",
    TRUE                                     ~ "Other"
  )
}

# 3. stations.json --------------------------------------------------------
# One record per analysis station, sorted by latitude descending.
# Cell density and observation counts come from all_data so this captures every
# confirmed presence, not just the filtered subset.

stations_payload <- all_data %>%
  ungroup() %>%
  mutate(
    code         = unname(station_map[combined_station]),
    present_flag = ifelse(probability == "present", 1L, 0L)
  ) %>%
  filter(!is.na(code)) %>%
  group_by(combined_station, code) %>%
  summarise(
    lat                = round(mean(lat, na.rm = TRUE), 4),
    lon                = round(mean(lon, na.rm = TRUE), 4),
    n_observations     = n(),
    n_present          = sum(present_flag, na.rm = TRUE),
    first_year         = suppressWarnings(min(year, na.rm = TRUE)),
    last_year          = suppressWarnings(max(year, na.rm = TRUE)),
    first_present_year = suppressWarnings(min(year[present_flag == 1], na.rm = TRUE)),
    max_cells_L        = suppressWarnings(max(cells_L[present_flag == 1], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Replace +/-Inf (no presences ever) with NA so JSON encodes as null.
    first_present_year = ifelse(is.finite(first_present_year), first_present_year, NA_integer_),
    max_cells_L        = ifelse(is.finite(max_cells_L), round(max_cells_L), NA_real_),
    basin              = basin_from_code(code),
    name               = combined_station
  ) %>%
  arrange(desc(lat)) %>%
  select(code, name, basin, lat, lon,
         n_observations, n_present, first_year, last_year,
         first_present_year, max_cells_L)

write_json(
  stations_payload,
  file.path(out_dir, "stations.json"),
  auto_unbox = TRUE, na = "null", pretty = FALSE
)

# 4. yearly_probability.json ----------------------------------------------
# Long-format: one row per station x year.
# Wilson 95% CI so the JSON carries everything needed to draw the ribbon
# without re-fitting on the client.

yearly_payload <- all_data %>%
  ungroup() %>%
  mutate(code = unname(station_map[combined_station])) %>%
  filter(!is.na(code), !is.na(probability), !is.na(year)) %>%
  group_by(code, year) %>%
  summarise(
    n         = n(),
    n_present = sum(probability == "present", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p     = n_present / n,
    # Wilson score 95% CI
    z     = 1.96,
    denom = 1 + z^2 / n,
    centre = (p + z^2 / (2 * n)) / denom,
    half   = (z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / denom,
    lo     = pmax(0, centre - half),
    hi     = pmin(1, centre + half),
    p      = round(p, 4),
    lo     = round(lo, 4),
    hi     = round(hi, 4)
  ) %>%
  select(code, year, n, n_present, p, lo, hi) %>%
  arrange(code, year)

write_json(
  yearly_payload,
  file.path(out_dir, "yearly_probability.json"),
  auto_unbox = TRUE, na = "null", pretty = FALSE
)

# 5. meta.json ------------------------------------------------------------
meta <- list(
  generated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  data_through   = suppressWarnings(max(all_data$year, na.rm = TRUE)),
  data_from      = suppressWarnings(min(all_data$year, na.rm = TRUE)),
  n_stations     = nrow(stations_payload),
  n_observations = sum(stations_payload$n_observations),
  n_present      = sum(stations_payload$n_present),
  citation       = "Moeller, K., Carstensen, J., Jakobsen, H., Engesmo, A., Karlson, B. (2026) Time series analysis of the toxic dinoflagellate Alexandrium pseudogonyaulax across Northern European waters. Harmful Algae 153, 103060."
)

write_json(
  meta,
  file.path(out_dir, "meta.json"),
  auto_unbox = TRUE, na = "null", pretty = TRUE
)

message("web_export/ written: ",
        paste(list.files(out_dir), collapse = ", "))
