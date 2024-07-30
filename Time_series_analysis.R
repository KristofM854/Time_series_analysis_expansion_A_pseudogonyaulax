##########################################
## Time series analysis of the expansion of Alexandrium pseudogonyaulax in northern European waters as part of Kristof's Phd thesis
## by Kristof Möller (Alfred Wegener Institut), Jacob Carstensen and Hans Jakobsen (each Aarhus University), Annette Engesmo (Norway) and Bengt Karlson (Sweden)
## Published in Limnology & Oceanography: 
## All raw-data available on PANGAEA: 
## Questions to: kristof-moeller@outlook.de
## Kristof Möller 06/24
## Alfred-Wegener-Institute Bremerhaven
##########################################
# Load custom functions #####
source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

# Install needed packages
install_packages()

####### NORWAY Data #########
# Load data #####
# Phytoplankton counts data
Nor_counts <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Norway\\phyto_feb24.txt",
    delim = "\t",
    col_names = TRUE
  )
Nor_CTD <- # including water profiles
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Norway\\norway_ctd_full.txt",
    delim = "\t",
    col_names = TRUE
  )
# Physical and chemical parameters of the water column
Nor_waterquality <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Norway\\nutrients.txt",
    delim = "\t",
    col_names = TRUE
  )

# General data transformations : CTD-data #####
# change long into wide format of all Norwegian data sets and remove unwanted columns
Nor_CTD <-
  pivot_wider(Nor_CTD[,c(9, 36:46)], values_from = verdinum, names_from = Parameter_navn, values_fn = mean) 
Nor_CTD <- Nor_CTD %>% dplyr::select(-Oksygenmetning, -Oksygen, -Tetthet)

# change colnames of CTD data
colnames(Nor_CTD) <- c(
  "station",
  "depth_start",
  "depth_end",
  "date",
  "year",
  "month",
  "day",
  "doy",
  "lat",
  "lon",
  "temp",
  "sal"
) 

# calculate the density column
Nor_CTD <- Nor_CTD %>% drop_na(temp) %>%
  mutate(density = 999.842594 + 6.793952E-2 * temp - 9.095290E-3 * temp^2 + 1.001685E-4 * temp^3 - 1.120083E-6 * temp^4 + 6.536332E-9 * temp^5 +
           (8.24493E-1 - 4.0899E-3 * temp + 7.6438E-5 * temp^2 - 8.2467E-7 * temp^3 + 5.3875E-9 * temp^4) * sal +
           (-5.72466E-3 + 1.0227E-4 * temp - 1.6546E-6 * temp^2) * sal^1.5 + 4.8314E-4 * sal^2
  )

# combine both water depth in one (One sometimes contains NAs)
Nor_CTD <-
  Nor_CTD %>%  mutate(depth = rowMeans(dplyr::select(., one_of(
    "depth_start", "depth_end"
  )), na.rm = TRUE)) %>% dplyr::select(-depth_start, -depth_end)

# Introduce stratification key if density difference of first two metres vs last two metres exceeds 1 g/cm^3
Nor_CTD <- Nor_CTD %>% 
  group_by(station, date) %>% 
  dplyr::mutate(strat = as.factor(ifelse(
    abs(mean(c(density[depth >= max(depth) - 2], density[max(depth)]), na.rm = TRUE) -
          mean(c(density[depth <= min(depth) + 2], density[min(depth)]), na.rm = TRUE)) >= 1,
    "stratified",
    "not stratified"
  ))) %>%
  ungroup()

# Average physical parameters of CTD data over first 10m of the  water column 
Nor_CTD <- Nor_CTD %>%
  group_by(date, year, doy, day, month, strat, station) %>%
  filter(depth < 10) %>%
  dplyr::summarise(sal = mean(sal, na.rm = TRUE),
                   temp = mean(temp, na.rm = TRUE),
                   density = mean(density, na.rm = TRUE)) %>%
  ungroup()

# General data transformations : microalgae counts data #####
# keep relevant columns and change to English names
Nor_counts <- Nor_counts[, c(1:2, 6:10, 15, 17:18)]

colnames(Nor_counts) <-
  c("station",
    "date",
    "year",
    "month",
    "day",
    "species",
    "cells_L",
    "doy",
    "lat",
    "lon")

# Perform left join and select only the lat and lon columns from Nor_counts
Nor_waterquality <- left_join(Nor_waterquality,
                              Nor_counts %>% dplyr::select(station, lat, lon) %>% unique(),
                              by = c("Vannlokalitet" = "station"))

# Introduce probability key and change all non-Alexandrium pseudogonyaulax entries to "No Alexandrium"
Nor_counts <- Nor_counts %>%
  mutate(
    probability = case_when(
      species == "Alexandrium pseudogonyaulax" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AO = case_when(
      species == "Alexandrium ostenfeldii" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AT = case_when(
      species == "Alexandrium tamarense" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AM = case_when(
      species == "Alexandrium minutum" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    species = case_when(
      species == "Alexandrium pseudogonyaulax" ~ "Alexandrium pseudogonyaulax",
      species == "Alexandrium ostenfeldii" ~ "Alexandrium ostenfeldii",
      species == "Alexandrium tamarense" ~ "Alexandrium tamarense",
      species == "Alexandrium minutum" ~ "Alexandrium minutum",
      str_detect(species, "Alexandrium") ~ "Alexandrium spp.",
      is.na(species) ~ NA_character_,
      TRUE ~ "No Alexandrium"
    ))

Nor_counts <- Nor_counts %>%
  mutate(probability_Aspp = case_when(
    str_detect(species, "Alexandrium spp.") & !is.na(cells_L) ~ "present",
    is.na(species) | is.na(cells_L) ~ NA_character_,
    TRUE ~ "absent"
  ))

# remove any absent entries on dates where A. pseudogonyaulax was actually present
# stems from the operation before that changed all other species to "No Alexandrium"
# thus on all present days "No Alexandrium" exists as well and needs to be removed
Alex_intermediate <- Nor_counts %>% filter(!species %in% c("No Alexandrium", "Alexandrium pseudogonyaulax"))

Nor_counts <- Nor_counts %>%
  filter(species %in% c("No Alexandrium", NA, "Alexandrium pseudogonyaulax")) %>%
  group_by(station, date) %>%
  filter(!(any(probability == "present") & probability == "absent")) %>%
  slice(ifelse(any(probability == "present"), which(probability == "present"), 1)) %>%
  ungroup()

Nor_counts <- full_join(Nor_counts, Alex_intermediate) 

# General data transformations : chemical and physical parameters #####
# change column names
colnames(Nor_waterquality) <- c(
  "station",
  "year",
  "month",
  "day",
  "doy",
  "PO4",
  "NO3",
  "silicate",
  "nprat",
  "siprat",
  "sinrat",
  "date",
  "lat",
  "lon"
)

# Merge all Norwegian data files #####
# Only keep first station of combined stations so that they match each other
Nor_counts$station <- sapply(strsplit(Nor_counts$station, ", "), function(x) x[1])
Nor_CTD$station <- sapply(strsplit(Nor_CTD$station, ", "), function(x) x[1])
Nor_waterquality$station <- sapply(strsplit(Nor_waterquality$station, ", "), function(x) x[1])

NOR <- Nor_counts %>%
  full_join(Nor_CTD) %>%
  full_join(Nor_waterquality)

NOR <- NOR %>% mutate(month = month(date))

## NO MATCHES BETWEEN NUTRIENT DATA AND CELL COUNTS: CONSIDER USING ADJACENT NUTRIENT SAMPLING DAYS??

# introduce logistic column (0 = absent and 1 = present) for logistic regression
NOR <-
  NOR %>%
  mutate(probability = ifelse(species == "Alexandrium pseudogonyaulax", "present", "absent")) %>%
  mutate(logistic = ifelse(NOR$probability == "absent", 0, 1)) %>% convert_as_factor(probability)

# Average numeric data from the same stations and dates 
NOR <-
  NOR %>% group_by(
    station,
    date,
    day,
    month,
    year,
    doy,
    lat,
    lon,
    probability,
    probability_AT,
    probability_AO,
    probability_Aspp,
    probability_AM,
    strat,
    species
  ) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) 

# Nutrient data often exists one day before or after phytoplankton data. Combine these 
NOR <- average_adjacent_days(NOR, "station") %>% unique()

# # Save NOR as a new txt.file #####
write.table(NOR, file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/NOR.txt", sep = "\t", row.names = FALSE)

####### DANISH Data #########
# Read in data files #####
# Phytoplankton counts data
DK_counts <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\DK_counts.txt",
    delim = "\t",
    col_names = TRUE
  ) 
DK_counts <- DK_counts %>% mutate(
  date = as.Date(paste(day, month, year, sep = "."), format = "%d.%m.%Y"))
# CTD data
DK_CTD <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\DK_CTD.txt",
    delim = "\t",
    col_names = TRUE
  )
DK_CTD <- DK_CTD %>% mutate(date = as.Date(date, format = "%d. %b %y")) %>% mutate(day = day(date))

# remove the -000X from the Limfjord stations as it does not match the synthax of DK_counts
DK_CTD$station <- sapply(strsplit(DK_CTD$station, "-"), function(x) x[1])

# Chemical and physical water properties
DK_waterquality <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\DK_water_quality.txt",
    delim = "\t",
    col_names = TRUE
  )
DK_waterquality <-
  DK_waterquality %>% mutate(date = as.Date(date, format = "%d. %b %y")) %>% mutate(day = day(date))

# remove the -000X from the Limfjord stations as it does not match the synthax of DK_counts
DK_waterquality$station <-
  sapply(strsplit(DK_waterquality$station, "-"), function(x)
    x[1])

# Secci depth and Kd values
DK_secci_kd <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\DK_secci_kd.txt",
    delim = "\t",
    col_names = TRUE
  )

DK_secci_kd <-
  DK_secci_kd %>% mutate(date = as.Date(date, format = "%d. %b %y")) %>% mutate(day = day(date))
DK_secci_kd$station <-
  sapply(strsplit(DK_secci_kd$station, "-"), function(x)
    x[1])

# General data transformations: CTD-Data #####
# change column names
DK_CTD <- DK_CTD %>%
  dplyr::rename(
    depth = `depth (m)`,
    sal = `Salinity (ppt)`,
    temp = `°C`,
    Chl = `Chlorophyl (flourometer,µg Chl L-1)`
  )

DK_counts <-
  DK_counts %>% dplyr::rename(
    lon = Longitude,
    lat = Latitude,
    species = species_name,
    cells_L = cells_per_L,
  )

# Replace lat/lon of the stations in DK_CTD and DK_waterquality with lat/lon of the stations in DK_counts
# get lat lon of all stations in DK_counts
all_stations <- DK_counts %>% dplyr::select(station, lon, lat) %>% unique() %>%
  separate_rows(station, sep = ", ") %>%
  distinct() %>%
  group_by(station) %>%
  dplyr::reframe(lat = mean(lat), lon = mean(lon))

# Join with DK_CTD
DK_CTD <- left_join(DK_CTD,
                    all_stations,
                    by = "station")

# Perform left join and select only the lat and lon columns from DK_counts
DK_waterquality <- left_join(DK_waterquality,
                             (all_stations),
                             by = "station")

# Perform left join and select only the lat and lon columns from DK_counts
DK_secci_kd <- left_join(DK_secci_kd,
                         (all_stations),
                         by = "station")
# calculate density 
DK_CTD <- DK_CTD %>%
  mutate(density = 999.842594 + 6.793952E-2 * temp - 9.095290E-3 * temp^2 + 1.001685E-4 * temp^3 - 1.120083E-6 * temp^4 + 6.536332E-9 * temp^5 +
           (8.24493E-1 - 4.0899E-3 * temp + 7.6438E-5 * temp^2 - 8.2467E-7 * temp^3 + 5.3875E-9 * temp^4) * sal +
           (-5.72466E-3 + 1.0227E-4 * temp - 1.6546E-6 * temp^2) * sal^1.5 + 4.8314E-4 * sal^2
  )

# Introduce stratification key if density difference exceeds 1 g/cm^3 
DK_CTD <- DK_CTD %>% 
  group_by(station, date) %>% 
  dplyr::mutate(strat = as.factor(ifelse(
    abs(mean(c(density[depth >= max(depth) - 2], density[max(depth)]), na.rm = TRUE) -
          mean(c(density[depth <= min(depth) + 2], density[min(depth)]), na.rm = TRUE)) >= 1,
    "stratified",
    "not stratified"
  ))) %>%
  ungroup()

# Average physical parameters of CTD data over the whole water column (0-10m)
DK_CTD <- DK_CTD %>%
  group_by(date, year, month, day, strat, station, lat, lon) %>%
  filter(depth < 10) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup()

# General data transformations: Microalgae counts data #####
# Replace occurrences of "Alexandrium_pseudogoniaulax" with "Alexandrium_pseudogonyaulax" and introduce probability key 
DK_counts$species <-
  gsub(
    "_",
    " ",
    DK_counts$species
  )

DK_counts$species <-
  gsub(
    "Alexandrium pseudogoniaulax",
    "Alexandrium pseudogonyaulax",
    DK_counts$species
  )

# Replace all other phytoplankton occurences with "No Alexandrium"
DK_counts <- DK_counts %>%
  mutate(
    probability = case_when(
      species == "Alexandrium pseudogonyaulax" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AO = case_when(
      species == "Alexandrium ostenfeldii" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AT = case_when(
      species == "Alexandrium tamarense" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AM = case_when(
      species == "Alexandrium minutum" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    species = case_when(
      species == "Alexandrium pseudogonyaulax" ~ "Alexandrium pseudogonyaulax",
      species == "Alexandrium ostenfeldii" ~ "Alexandrium ostenfeldii",
      species == "Alexandrium tamarense" ~ "Alexandrium tamarense",
      species == "Alexandrium minutum" ~ "Alexandrium minutum",
      str_detect(species, "Alexandrium") ~ "Alexandrium spp.",
      is.na(species) ~ NA_character_,
      TRUE ~ "No Alexandrium"
    ))

DK_counts <- DK_counts %>%
  mutate(probability_Aspp = case_when(
  str_detect(species, "Alexandrium spp.") & !is.na(cells_L) ~ "present",
  is.na(species) | is.na(cells_L) ~ NA_character_,
  TRUE ~ "absent"
))

# Only keep one "No Alexandrium" entry per date
Alex_intermediate <- DK_counts %>% filter(!species %in% c("No Alexandrium", NA, "Alexandrium pseudogonyaulax"))

DK_counts <- DK_counts %>%
  filter(species %in% c("No Alexandrium", NA, "Alexandrium pseudogonyaulax")) %>%
  group_by(station, date, year, month, day) %>%
  filter(!(any(probability == "present") & probability == "absent")) %>%
  slice(ifelse(any(probability == "present"), which(probability == "present"), 1)) %>%
  ungroup()

DK_counts <- full_join(DK_counts, Alex_intermediate)

# Group by 'station' and 'date' and summarize all numeric columns
DK_CTD <- DK_CTD %>%
  group_by(station, date, strat) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

# Group by 'station' and 'date' and summarize all numeric columns
DK_waterquality <- DK_waterquality %>%
  group_by(station, date, lat, lon) %>%
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

# Merge all data files #####
DK <- DK_counts %>%
  full_join(DK_waterquality) %>%
  full_join(DK_CTD) %>%
  full_join(DK_secci_kd) 

# combine plankton counts and nutrients data if thye only differ in one day
DK <- average_adjacent_days(DK, "station") %>% unique()

# Convert date, introduce weeks and day of the year (doy)
DK <- DK %>%
  mutate(
    week = week(date),
    doy = yday(date)
  )

# Read wind data files #####
DK_wind <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\wind_data.txt",
    delim = "\t",
    col_names = TRUE
  )
DK_wind_stations <-
  fread(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Denmark\\stationen_wind.txt",
    fill = TRUE
  )

# General data transformations wind data #####
# change date format, introduce doy, and change comma to points and change to numeric
DK_wind <- DK_wind %>% mutate(
  date =
    as.Date(
      paste(day, month, year, sep = "-"),
      format = "%d-%m-%Y"
    ),
  doy = yday(date),
  `wind_m/s` = as.numeric(gsub(",", ".", `wind_m/s`))
) %>% dplyr::select(-release_date)

# remove last two digits of the station name in the DK_wind dataset to match the format in the station file
DK_wind <- DK_wind %>%
  mutate_at(vars(station), ~ as.numeric(str_sub(., end = -3)))

# select columns of interest
DK_wind_stations <- DK_wind_stations %>%
  dplyr::select(station, Lat, Lon)

# Merge both wind dataframes
DK_wind <- full_join(DK_wind, DK_wind_stations) %>%
  drop_na(Lat)

# Find the two closest stations between danish monitoring stations and meterological "wind" stations
closest_stations <-
  find_closest_stations(
    DK_wind %>% distinct(station, .keep_all = TRUE),
    DK_counts %>% distinct(station, .keep_all = TRUE)
  )

# Merge DK with DK_wind and the closest stations
# adjust column names to match DK
colnames(closest_stations) <- c("wind_match", "station", "distance")

DK <- merge(DK, closest_stations, by = "station", all.x = T)

colnames(DK_wind)[1] <- "wind_match"

# Merge DK and DK wind and finalize #####
DK <-
  merge(
    DK,
    DK_wind %>% dplyr::select(wind_match, date, year, month, day, `wind_m/s`),
    by = c("wind_match", "date", "year", "month", "day"),
    all.x = T,
    all.y = T
  )

# Introduce probability key
DK$probability <-
  as.factor(ifelse(
    DK$species == "Alexandrium pseudogonyaulax",
    "present",
    ifelse(is.na(DK$species), NA , "absent")
  ))

# Rename columns
DK <- DK %>% dplyr::rename(C = "C_(µgC_L-1)",
                           DIN = "DIN_(µgL-1)",
                           NH4 = "NH4_(µgL-1)",
                           NO3 = "nitrat+nitrit_(µgL-1)",
                           TN = "Total_nitrogen_(µgL-1)",
                           DIP = `DIP_(µgL-1)`,
                           PO4 = `Total_phosphate_(µgL-1)`,
                           chl = `Chlorophyll_a_(µgr_L-1)`,
                           silicate = `Silicate_(µgL-1)`,
                           chl_2 = "Chl",
                           kd = "Kd",
                           secci = "Secci",
                           wind_ms = `wind_m/s`) %>%
  dplyr::select(-Phylum, -Class, -'Cell_vol(µm3)', -distance, -wind_match)

# Convert ug/L to umol/L  to match other datasets
DK <-
  DK %>% mutate(
    C = C / 12.0107,
    DIN = DIN / 14.006720,
    NH4 = NH4 / 18.0383,
    NO3 = NO3 / 62.0049,
    TN = TN / 14.006720,
    DIP = DIP / 30.973762 ,
    PO4 = PO4 / 94.9712
  ) %>% drop_na(station, date) # remove wind data without any matches

# # Save DK as a new txt.file #####
write.table(DK, file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/DK.txt", sep = "\t", row.names = FALSE)

######### SWEDISH DATA #########
# Load Swedish data #####
SW_phys2 <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Sweden\\sharkweb_phys2.txt",
    delim = "\t",
    col_names = TRUE
  )
SW_phyto <-
  read_delim(
    "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Sweden\\sharkweb_phyto_new.txt",
    delim = "\t",
    col_names = TRUE
  )

# # General Data Preparation for SW_phyto #####
SW_phyto <-
  SW_phyto %>% dplyr::select(
    reported_station_name,
    sample_date,
    sample_latitude_dm,
    sample_longitude_dm,
    wind_speed_ms,
    secchi_depth_m,
    sample_min_depth_m,
    sample_max_depth_m,
    scientific_name,
    parameter,
    value,
    unit,
    coefficient,
    sedimentation_volume_ml
  )

colnames(SW_phyto) <-
  c("station",
    "date",
    "lat",
    "lon",
    "wind_ms",
    "secci",
    "min_depth",
    "max_depth",
    "name",
    "parameter",
    "counts",
    "unit",
    "counts_coefficient",
    "sedimentation_volume")

# Replace all commas with dots in all columns
SW_phyto[,c(5:14)] <-
  sapply(SW_phyto[,c(5:14)], function(col)
    gsub(",", ".", col))

# Convert counts and counts_coefficient to numeric
SW_phyto <- SW_phyto %>%
  mutate(
    counts = as.numeric(counts),
    counts_coefficient = as.numeric(counts_coefficient)
  )

# Perform the multiplication conditionally
SW_phyto <- SW_phyto %>%
  mutate(
    counts = ifelse(grepl("counted", parameter, ignore.case = TRUE) & !is.na(counts_coefficient), 
                    counts * counts_coefficient, 
                    counts)
  )

# Convert LAT and LON to decimal degrees
SW_phyto <- SW_phyto %>%
  mutate(lat = sapply(lat, dmm_to_dd),
         lon = sapply(lon, dmm_to_dd),
         date = as.Date(date))

# General Data transformations for SW_phys #####
# Keep columns of interest and change to english names
SW_phys2 <- SW_phys2[, c(9, 15, 20:21, 26:28, 34, 38:39, 41, 43, 45, 47, 69, 71, 73, 75, 77, 79, 81, 83, 87, 89, 91, 93, 95, 97, 99, 109)]

SW_phys2 <- SW_phys2 %>%
  mutate(lat = sapply(`Provets latitud (DM)`, dmm_to_dd),
         lon = sapply(`Provets longitud (DM)`, dmm_to_dd),
         Provtagningsdatum = as.Date(Provtagningsdatum, format = "%Y-%m-%d")) %>%
  dplyr::select(-`Provets latitud (DM)`, -`Provets longitud (DM)`)

colnames(SW_phys2) <-
  c(
    "station",
    "date",
    "wind_direction",
    "wind_speed",
    "air_temp",
    "secci",
    "sampling_depth",
    "pressure",
    "temp",
    "temp2",
    "sal",
    "sal2",
    "PO4",
    "TP",
    "NO2",
    "NO3",
    "NO2+NO3",
    "NH4",
    "TN",
    "silicate",
    "chl",
    "DOC",
    "POC",
    "TOC",
    "PON",
    "cur",
    "cur_dir",
    "CDOM",
    "lat",
    "lon"
  )

# Convert columns with commas to numeric and change , to .
SW_phys2[, c(3:28)] <-
  sapply(SW_phys2[, c(3:28)], function(col)
    as.numeric(gsub(",", ".", col)))


# calculate the mean of both salinity and temperature columns excluding NAs and remove the redundant ones 
SW_phys2 <- SW_phys2 %>%
  ungroup() %>%
  mutate(
    temp = rowMeans(dplyr::select(., one_of("temp", "temp2")), na.rm = TRUE),
    sal = sal / 1000,
    `NO2+NO3` = coalesce(NO3, 0) + coalesce(NO2, 0)
  ) %>%
  dplyr::select(-temp2, -sal2) %>%
  mutate(`NO2+NO3` = replace(`NO2+NO3`, is.na(NO3) &
                               is.na(NO2), NA)) %>%
  mutate(`NO2+NO3` = replace(`NO2+NO3`, NO2 == 0 & is.na(NO3), NA)) %>%
  mutate(DIN = NH4 + `NO2+NO3`) %>% # Calculate the DIN (Dissolved Inorganic Nitrogen) concentration
  dplyr::select(-NO2, -NO3) %>%
  dplyr::rename(NO3 = 'NO2+NO3')

# calculate the density column
SW_phys2 <- SW_phys2 %>%
  mutate(density = 999.842594 + 6.793952E-2 * temp - 9.095290E-3 * temp^2 + 1.001685E-4 * temp^3 - 1.120083E-6 * temp^4 + 6.536332E-9 * temp^5 +
           (8.24493E-1 - 4.0899E-3 * temp + 7.6438E-5 * temp^2 - 8.2467E-7 * temp^3 + 5.3875E-9 * temp^4) * sal +
           (-5.72466E-3 + 1.0227E-4 * temp - 1.6546E-6 * temp^2) * sal^1.5 + 4.8314E-4 * sal^2
  )

# Introduce stratification key if density difference exceeds 1 g/cm^3 (?!)
SW_phys2 <- SW_phys2 %>% 
  group_by(station, date) %>% 
  dplyr::mutate(strat = as.factor(ifelse(
    abs(mean(c(density[sampling_depth >= max(sampling_depth) - 2], density[max(sampling_depth)]), na.rm = TRUE) -
          mean(c(density[sampling_depth <= min(sampling_depth) + 2], density[min(sampling_depth)]), na.rm = TRUE)) >= 1,
    "stratified",
    "not stratified"
  ))) %>%
  ungroup()

# # General Data Preparation for SW_phyto #####
# Average over the first 10m of the water column
SW_phyto <- SW_phyto %>% dplyr::rename(cells_L = counts)

SW_phyto[, c(5:8)] <-
  sapply(SW_phyto[, c(5:8)], function(col)
    as.numeric(col))

SW_phyto <-
  SW_phyto %>% filter(max_depth <= 10) %>% 
  group_by(station, date, name, lat, lon, parameter) %>% 
  reframe(wind_ms = mean(wind_ms, na.rm = T),
          cells_L = mean(cells_L, na.rm = T),
          secci = mean(secci, na.rm = T)) %>%
  filter(str_detect(parameter, c("Abundance|counted"))) %>%
  filter(parameter != "Abundance class")

# Introduce probability key and change all non-Alexandrium pseudogonyaulax entries to "No Alexandrium"
SW_phyto <- SW_phyto %>%
  mutate(
    probability = case_when(
      name == "Alexandrium pseudogonyaulax" & !is.na(cells_L) ~ "present",
      is.na(name) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AO = case_when(
      name == "Alexandrium ostenfeldii" & !is.na(cells_L) ~ "present",
      is.na(name) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AT = case_when(
      name == "Alexandrium tamarense" & !is.na(cells_L) ~ "present",
      is.na(name) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AM = case_when(
      name == "Alexandrium minutum" & !is.na(cells_L) ~ "present",
      is.na(name) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    name = case_when(
      name == "Alexandrium pseudogonyaulax" ~ "Alexandrium pseudogonyaulax",
      name == "Alexandrium ostenfeldii" ~ "Alexandrium ostenfeldii",
      name == "Alexandrium tamarense" ~ "Alexandrium tamarense",
      name == "Alexandrium minutum" ~ "Alexandrium minutum",
      str_detect(name, "Alexandrium") ~ "Alexandrium spp.",
      is.na(name) ~ NA_character_,
      TRUE ~ "No Alexandrium"
    ))

SW_phyto <- SW_phyto %>%
  mutate(probability_Aspp = case_when(
    str_detect(name, "Alexandrium spp.") & !is.na(cells_L) ~ "present",
    is.na(name) | is.na(cells_L) ~ NA_character_,
    TRUE ~ "absent"
  ))

# remove any absent entries on dates where A. pseudogonyaulax was actually present
# Before that filter out all other Alexandrium species and join them after
Alex_intermediate <- SW_phyto %>% filter(!name %in% c("No Alexandrium", NA, "Alexandrium pseudogonyaulax"))

SW_phyto <- SW_phyto %>%
  filter(name %in% c("No Alexandrium", NA, "Alexandrium pseudogonyaulax")) %>%
  group_by(station, date) %>%
  filter(!(any(probability == "present") & probability == "absent")) %>%
  slice(ifelse(any(probability == "present"), which(probability == "present"), 1)) %>%
  ungroup()

SW_phyto <- full_join(SW_phyto, Alex_intermediate)

# Average the physical parameters of interest over a sampling depth of 10m to mimic danish data set
# in chunks of the dataframe because my laptop is ******* slow
# Calculate number of chunks
chunk_size <- 25000
num_chunks <- ceiling(nrow(SW_phys2) / chunk_size)

# Initialize an empty list to store the summarized dataframes
summarized_dfs <- list()

# Loop through each chunk, summarize it, and store the result
for (i in 1:num_chunks) {
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, nrow(SW_phys2))
  
  chunk <- SW_phys2[start_row:end_row, ]
  
  summarized_chunk <- chunk %>%
    drop_na(station) %>%
    group_by(station, date, strat) %>%
    filter(sampling_depth < 10) %>%
    dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
  
  # Store the summarized chunk in the list
  summarized_dfs[[i]] <- summarized_chunk
}

# Combine all summarized dataframes 
SW_phys2 <- do.call(rbind, summarized_dfs)

# Combine microalgae and physical parameters datasets
SW_work <-
  full_join(SW_phyto,
            SW_phys2) 

SW_work <-
  SW_work %>% mutate(wind_ms2 = rowMeans(dplyr::select(., starts_with("wind")), na.rm = T))

# Introduce additional time columns
SW_work <- SW_work %>%
  mutate(
    month = month(date),
    year = year(date),
    doy = yday(date),
    week = week(date),
    day = day(date)
  ) %>%
  as.data.frame() %>%
  dplyr::select(-wind_ms, -wind_speed)

# Rename species column to match other datasets
SW_work <-
  SW_work %>% dplyr::rename(species = name, wind_ms = wind_ms2)

# combine plankton counts and nutrients data if they only differ in one day
SW_work <- average_adjacent_days(SW_work, "station") %>% unique()

# # Save SW_work as a new txt.file #####
write.table(SW_work, file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/SW_work.txt", sep = "\t", row.names = FALSE)

########## IOW-Odin data #########
# Load data #####
HD <- read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Odin-Baltic Sea\\odin2_2024-01-31_095801.txt",
  delim = "\t",
  skip = 2,  # Skip the first line
  col_names = FALSE,
  col_types = cols(.default = "c"))

# Read the first two lines
first_two_lines <- read_lines(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Odin-Baltic Sea\\odin2_2024-01-31_095801.txt",
  n_max = 2
)

# Split the lines into vectors based on tabs
line1_vector <- strsplit(first_two_lines[1], "\t")[[1]]
line2_vector <- strsplit(first_two_lines[2], "\t")[[1]]

# Concatenate corresponding values from line1 and line2
combined_header <- paste(line1_vector, line2_vector, sep = "_")

# Set the combined header to the dataframe
colnames(HD) <- combined_header

# Load additional physical parameters (wind data and secci depth) #####
station_details = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\Odin-Baltic Sea\\Station_details.txt",
  delim = "\t",
  col_names = T
)

# General data transformations #####
# remove all columns containing biomass or carbon in the name as we are only interested in cell counts
HD <- HD %>%
  dplyr::select(-starts_with("NA_Carbon_"), -starts_with("NA_Biomass_"))

# Rename columns
# Function to remove "NA_" prefix
remove_na_prefix <- function(col_name) {
  gsub("^NA_", "", col_name)
}
HD <- HD %>%
  rename_with(remove_na_prefix, starts_with("NA_"))

# Replace , with .
HD[, ] <- sapply(HD[,], gsub, pattern = ",", replacement = ".")

# Apply transformation to numeric for columns consisting of numbers
HD[,3:604] <- sapply(HD[,3:604], as.numeric)

# Change date format
HD <- HD %>% mutate(`Time_(start)` = as.Date(HD$`Time_(start)`, format = "%d.%m.%Y"))

# Introduce month, year, and doy columns and convert all numeric data to numeric
HD <- HD %>%
  mutate(month = month(`Time_(start)`),
         year = year(`Time_(start)`),
         doy = yday(`Time_(start)`)) %>% 
  dplyr::rename(station = Name, lon = `Longitude_(start)_[°]`, lat = `Latitude_(start)_[°]`)

# Calculate means of temperature and salinity, excluding NAs
HD <- HD %>%
  mutate(
    temp = rowMeans(dplyr::select(., starts_with("Temp")), na.rm = TRUE),
    sal = rowMeans(dplyr::select(., starts_with("Sal")), na.rm = TRUE)) %>%
  dplyr::select(-temp1, -sal1, -temp2, -sal2, -temp3, -sal3)

HD <- HD %>% mutate_at(c("temp", "sal"), as.numeric)

# Pivot all cell count columns into a single column with the id in the species column
HD <- HD %>%
  pivot_longer(cols = contains("1/m**3"),
               names_to = "species",
               values_to = "cells_L")

HD <- HD %>%
  mutate(`NO3+NO2` = rowMeans(dplyr::select(., starts_with("NO")), na.rm = TRUE))

# Split station names into shorter strings
HD$station <- sapply(strsplit(HD$station, "_"), function(x) x[1])

# Combine species (cf and non-cf) counts and shorten names
HD <- HD %>%
  mutate(species = gsub(
    "Alexandrium_pseudogonyaulax.*",
    "Alexandrium pseudogonyaulax",
    species
  ))
HD <- HD %>%
  mutate(species = gsub(
    "Alexandrium_ostenfeldii.*",
    "Alexandrium ostenfeldii",
    species
  ))

# Introduce probability key and change all non-Alexandrium pseudogonyaulax entries to "No Alexandrium"
HD <- HD %>%
  mutate(
    probability = case_when(
      species == "Alexandrium pseudogonyaulax" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AO = case_when(
      species == "Alexandrium ostenfeldii" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AT = case_when(
      species == "Alexandrium tamarense" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    probability_AM = case_when(
      species == "Alexandrium minutum" & !is.na(cells_L) ~ "present",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "absent"
    ),
    species = case_when(
      species == "Alexandrium pseudogonyaulax" ~ "Alexandrium pseudogonyaulax",
      species == "Alexandrium ostenfeldii" ~ "Alexandrium ostenfeldii",
      species == "Alexandrium tamarense" ~ "Alexandrium tamarense",
      species == "Alexandrium minutum" ~ "Alexandrium minutum",
      str_detect(species, "Alexandrium") ~ "Alexandrium spp.",
      is.na(species) | is.na(cells_L) ~ NA_character_,
      TRUE ~ "No Alexandrium"
    ))

HD <- HD %>%
  mutate(probability_Aspp = case_when(
    str_detect(species, "Alexandrium spp.") & !is.na(cells_L) ~ "present",
    is.na(species) | is.na(cells_L) ~ NA_character_,
    TRUE ~ "absent"
  ))

# change station and date colnames
HD <- HD %>% dplyr::rename(date = `Time_(start)`)

# Only keep distinct rows as in the IOW dataset each species contains a row even though it was not detected or counted

HD <- HD %>% distinct()

# keep all present entries of Alexandrium
Alex_intermediate <-
  HD %>% filter(
    probability == "present" |
      probability_AO == "present" |
      probability_AT == "present" | 
      probability_Aspp == "present" |
      probability_AM == "present" 
  )

HD <- HD %>%
  filter(species %in% c("No Alexandrium", NA)) %>%
  group_by(station, date) %>% 
  arrange(cells_L) %>%
  slice_head() %>%
  ungroup()

HD <- full_join(HD, Alex_intermediate)

# Include secci depth and wind speed #####
# select columns of interest and change columns to english and to match HD dataframe
station_details <-
  station_details %>% dplyr::select(c(1, 2, 4, 6, 12, 19))

colnames(station_details) <-
  c("station", "date", "lon", "lat", "wind_ms", "secci")

# extract date from time column and convert to date
station_details$date <-
  str_sub(station_details$date,
          start = 1,
          end = 10) %>% as.Date(format = "%d.%m.%Y")

# Average all numeric columns of the same station and date 
HD <-
  HD %>% ungroup() %>% group_by(station,
                                date,
                                probability,
                                probability_AO,
                                probability_AT,
                                probability_AM,
                                probability_Aspp,
                                species,
                                `Depth_(start)_[m]`,
                                `Depth_(end)_[m]`) %>% dplyr::summarise_if(is.numeric, mean, na.rm = T)

# # join algae (HD) and station_details dataframe
HD <- left_join(HD, station_details, by = c("station", "date", "lat", "lon"))

# calculate the density column
HD <- HD %>%
  mutate(density = 999.842594 + 6.793952E-2 * temp - 9.095290E-3 * temp^2 + 1.001685E-4 * temp^3 - 1.120083E-6 * temp^4 + 6.536332E-9 * temp^5 +
           (8.24493E-1 - 4.0899E-3 * temp + 7.6438E-5 * temp^2 - 8.2467E-7 * temp^3 + 5.3875E-9 * temp^4) * sal +
           (-5.72466E-3 + 1.0227E-4 * temp - 1.6546E-6 * temp^2) * sal^1.5 + 4.8314E-4 * sal^2
  )

# Introduce stratification key if density difference exceeds 1 g/cm^3 (?!)
HD <- HD %>% dplyr::rename(depth = `Depth_(start)_[m]`)

# Introduce DIN (dissolved inorganic nitrogen) as dissolved nitrogen (DN) - dissolved organic nitrogen (DON)
HD$DIN <- HD$DN - HD$DON

HD <- HD %>%
  dplyr::group_by(station, date) %>%
  dplyr::mutate(strat = as.factor(ifelse(
    mean(density[depth >= (max(depth) - 2)], na.rm = TRUE) - mean(density[depth <= 2], na.rm = TRUE) > 1,
    "stratified",
    "not stratified"
  )))

# average the parameters of interest over a sampling depth of 10m to be equal to Danish and Swedish dataset
# divide cell counts / 1000 as they are per m^3
HD$cells_L <- HD$cells_L / 1000

HD <- HD[,] %>%
  group_by(
    station,
    date,
    probability,
    probability_AO,
    probability_AT,
    probability_AM,
    probability_Aspp,
    strat,
    species
  ) %>%
  filter(as.numeric(`Depth_(end)_[m]`) <= 10) %>%
  group_modify( ~ dplyr::summarise(.x, across(where(is.numeric), mean, na.rm = TRUE)))

# deselect columns 
HD <-
  HD %>% dplyr::select(-c(
    `Time_(end)`,
    `Longitude_(end)_[°]`,
    `Latitude_(end)_[°]`,
    `Depth_(end)_[m]`,
    dens,
    NO3
  ))

# Rename columns
HD <- HD %>% dplyr::rename(NO3 = `NO3+NO2`, chl = Chl)

# combine plankton counts and nutrients data if they only differ in one day
HD$day <- day(HD$date)

# If an Alexandrium entry on a given date/station combination exists only keep these and remove No Alexandrium
HD <- HD %>%
  group_by(station, date) %>%
  filter(
    any(species != "No Alexandrium") & species != "No Alexandrium" |
      all(species == "No Alexandrium")
  ) %>%
  ungroup()

# export data
write.table(HD,
            file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/HD.txt",
            sep = "\t",
            row.names = FALSE)
# 
# # # Merge dataframes of all monitoring stations #####
# # Read in previously modified and saved dataframes
HD = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\HD.txt",
  delim = "\t",
  col_names = T
)
SW_work = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\SW_work.txt",
  delim = "\t",
  col_names = T
)
NOR = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\NOR.txt",
  delim = "\t",
  col_names = T
)
DK = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\DK.txt",
  delim = "\t",
  col_names = T
)

# update day, month, year and doy from date of all datasets
update_dates <- function(df) {
  df <- df %>% mutate(date = as.Date(date)) 
  df %>%
    mutate(
      day = day(date),
      doy = yday(date),
      month = month(date),
      year = year(date)
    )
}

# Convert data frames to data tables
setDT(DK)
setDT(SW_work)
setDT(HD)

data_frames <- list(NOR, DK, SW_work, HD)

# Update dates for all data frames in the list
updated_data_frames <- lapply(data_frames, update_dates)

all_data <-  as.data.frame(reduce(updated_data_frames, full_join))

# combine stations in close proximity
# First average lat lon of each station, as they have small differences
# Loop over each station, average lat lon, store in list
unique_stations <-
  all_data %>% ungroup() %>% dplyr::select(station, lat, lon) %>% unique() %>% drop_na(lat) %>%
  mutate(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T)) %>% distinct()

station_list <- list()

for (i in 1:nrow(unique_stations)) {
  station_data <- all_data %>% filter(station == unique_stations$station[i])
  station_data <-
    station_data %>% mutate(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
  station_list[[unique_stations$station[i]]] <- station_data
}

# combine each dataframe from the list
all_data <- do.call(rbind, station_list)

# find stations in close proximity; threshold = +/- 0.025 lat/lon
# dataframe contains combined stations, old stations and the new averaged lat / lon of combined stations
close_stations <- combine_close_stations(all_data)

# rename station to combined station or the join after does not work properly
close_stations <-
  close_stations %>% dplyr::rename(combined_station = station)

# # # save close stations to not compute them each time
# write.table(close_stations,
#             file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/close_stations.txt",
#             sep = "\t",
#             row.names = FALSE)

close_stations <- read.delim(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/close_stations.txt",
         sep = "\t")

# Update all_data with new stations, lat and lons 
all_data <-
  full_join(
    all_data %>% dplyr::select(-lat, -lon),
    close_stations,
    by = c("station" = "old_station")
  )

# overwrite "limiting_conditions"
all_data <- all_data %>% mutate(
  limiting_conditions = as.factor(case_when(
    is.na(DIN) | is.na(PO4) ~  NA,
    DIN < 2 & PO4 < 0.2 ~ "yes",
    DIN > 2 | PO4 > 0.2 ~ "no"
  ))
)
# remove rows with no plankton data
# all_data <- all_data %>% drop_na(species)

all_data <- all_data %>%
  group_by(
    combined_station,
    date,
    probability,
    probability_AM,
    probability_AO,
    probability_Aspp,
    probability_AT,
    logistic,
    species
  ) %>%
  dplyr::summarize(
    strat = combine_strat(strat),
    limiting_conditions = combine_limiting_conditions(limiting_conditions),
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = 'drop'
  )

all_data <- average_adjacent_days(all_data, 'combined_station') %>% unique()

# Create nutrient and wind lagged columns #####
# Define parameters of interest
parameters_of_interest <- c("NH4", "NO3", "sal", "temp", "chl", "wind_ms", "TN", "PO4", "silicate")

# Loop through parameters and lag weeks
n_lags <- 1:3
for (param in parameters_of_interest) {
  for (num_weeks in n_lags) {
    all_data <- create_lagged_columns(all_data, param, num_weeks)
  }
}
# 
# # Introduce geographical regions #####
# region_boxes <- data.frame(
#   region = c(
#     "eastern_baltic",
#     "western_baltic",
#     "danish_southsea",
#     "kateggat",
#     "skagerrak",
#     "limfjord",
#     "norwegian_coast",
#     "norwegian_coast"
#   ),
#   min_lat = c(54, 54, 55.5, 56.38, 57.49, 56.4, 57, 62),
#   max_lat = c(57.5, 55.5, 56.38, 57.49, 60, 57, 65, 75),
#   min_lon = c(13.4, 9.5, 9.5, 9.5, 8, 8.01, 4, 8.01),
#   max_lon = c(21, 13.4, 13.3, 13.3, 15, 9.7, 8, 30),
#   lon = c(13.4, 9.5, 9.5, 6, 8.7, 8.01, 4, 8.01),
#   lat = c(54, 54, 55.5, 56.38, 57.49, 56.4, 57, 65)
# )
# 
# # Assign regions to each station
# all_data <- assign_regions_to_stations(all_data, region_boxes)

# introduce logistic column (0 = absent and 1 = present) for logistic regression
all_data <-
  all_data %>%
  mutate(logistic = ifelse(all_data$probability == "absent", 0, 1)) %>% convert_as_factor(probability, logistic)

# # Save all_data as a new txt.file #####
write.table(all_data,
            file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/all_data_lag.txt",
            sep = "\t",
            row.names = FALSE)

# Garbage collection: call after large objects have been removed #####
gc()

# Delete workspace, clean environment/variables and restarte R if used memory piles up
# dev.off()
rm(list = ls())

# restart R
.rs.restartR()

# Load custom functions #####
source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

install_packages()
# 
# read in all_data to start here: #####
all_data = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\all_data_lag.txt",
  delim = "\t",
  col_names = T
)

# Introduce 5 year and 10 year column #####
# Define the breaks and labels for the desired year segments
breaks <- c(1997, 2007, 2012, 2017, 2022)  # Note the last value goes one year beyond the last segment

labels <- c("1997-2006", "2007-2011", "2012-2016", "2017-2021")

# Mutate the dataframe
all_data <- all_data %>%
  mutate(year2 = cut(year, breaks = breaks, labels = labels, right = F))

# Find stations that meet the criteria (sampling in the years of 2010:2020 with maximum two years missing; in the month 4-10, and minimum 20 observations per year) ######
stations_to_keep <- all_data %>%
  filter(year %in% 2010:2020, month %in% 5:9) %>%
  drop_na(probability) %>% # Filter data for the specified years and months
  group_by(combined_station) %>%
  dplyr::summarise(unique_years = n_distinct(year), unique_month = n_distinct(month)) %>%  # Count the number of unique years for each station
  filter(unique_years >= (2020 - 2010 + 1 - 3) & unique_month >= 3) %>%  # Keep stations with less than 3 missing years
  pull(combined_station)

stations_to_keep2 <- all_data %>%
  filter(year %in% 2010:2020) %>% 
  drop_na(probability) %>% # Filter data for the specified years and months
  group_by(combined_station, year) %>%
  dplyr::summarise( unique_doy = n_distinct(doy)) %>%  # Count the number of unique years for each station
  filter(unique_doy >= 10) %>%  # Filter for at least 10 observations per year
  ungroup() %>%
  pull(combined_station) %>%
  unique()

# Find stations that meet the criteria
common_stations <- intersect(stations_to_keep, stations_to_keep2)

# Remove stations from all_data that do not meet the criteria above
filtered_data <- all_data %>%
  filter(combined_station %in% common_stations)

# convert chlorophyll to mol/L
filtered_data$chl <- filtered_data$chl / 893.51

# Introduce nutrient ratios 
filtered_data <-
  filtered_data %>% mutate('NO3_PO4' = DIN / PO4,
                           'Si_N' = silicate / DIN,
                           DON = TN - DIN,
                           'C_chl' = C / chl) %>% 
                    mutate('DON_DIN' = DON / DIN)

# replace combined station name with station
filtered_data <-
  filtered_data %>% dplyr::rename(station = combined_station)

unique_stations <-
  filtered_data %>% ungroup() %>% dplyr::select(station, lat, lon) %>% unique()
station_list <- list()

for (i in 1:nrow(unique_stations)) {
  station_data <- filtered_data %>% filter(station == unique_stations$station[i])
  station_data <-
    station_data %>% mutate(lat = mean(lat, na.rm = T), lon = mean(lon, na.rm = T))
  station_list[[unique_stations$station[i]]] <- station_data
}

# combine each dataframe from the list
filtered_data <- do.call(rbind, station_list)

# Change absent/present in probability columns to 0/1, which the following models accept
filtered_data <- filtered_data %>% mutate(probability = as.factor(ifelse(probability == "present", 1, ifelse(is.na(cells_L), NA, 0))),
                                          probability_AO = as.factor(ifelse(probability_AO == "present", 1, ifelse(is.na(cells_L), NA, 0))),
                                          probability_AT = as.factor(ifelse(probability_AT == "present", 1, ifelse(is.na(cells_L), NA, 0))),
                                          probability_AM = as.factor(ifelse(probability_AM == "present", 1, ifelse(is.na(cells_L), NA, 0))),
                                          probability_Aspp = as.factor(ifelse(probability_Aspp == "present", 1, ifelse(is.na(cells_L), NA, 0))))

# Group by station and year, then count the number of "present" observations in probability column
counts_overview <- all_data %>%
  filter(probability == "present" & combined_station %in% common_stations) %>%
  group_by(combined_station) %>%
  dplyr::summarise(
    present_count = n(),
    lat = as.numeric(sprintf("%.2f", first(lat))),
    lon = as.numeric(sprintf("%.2f", first(lon))),
    .groups = "drop" # Override grouped output
  )

counts_overview2 <- all_data %>%
  filter(probability == "present" & combined_station %in% common_stations) %>%
  group_by(combined_station) %>%
  dplyr::summarise(
    max_cells_L =  sprintf("%.2e", max(cells_L, na.rm = TRUE)),
    temp_range = paste0(ifelse(all(is.na(
      temp
    )), NA, round(
      min(temp, na.rm = TRUE), digits = 1
    )), "-", ifelse(all(is.na(
      temp
    )), NA, round(
      max(temp, na.rm = TRUE), digits = 1
    ))),
    sal_range = paste0(ifelse(all(is.na(
      sal
    )), NA, round(
      min(sal, na.rm = TRUE), digits = 1
    )), "-", ifelse(all(is.na(
      sal
    )), NA, round(
      max(sal, na.rm = TRUE), digits = 1
    ))),
    year_range = paste0(min(year), "-", max(year)),
    month_range = paste0(min(month), "-", max(month)))

counts_overview3 <- all_data %>%
  filter(combined_station %in% common_stations) %>%
  group_by(combined_station) %>%
  dplyr::summarise(
    sampling_range = paste0(min(year), "-", max(year)),
    lat = as.numeric(sprintf("%.2f", first(lat))),
    lon = as.numeric(sprintf("%.2f", first(lon))))

# Update dates for all data frames in the list
station_characteristics <-
  as.data.frame(reduce(
    list(counts_overview, counts_overview2, counts_overview3),
    full_join
  )) %>% arrange(desc(lat))

# reorder column order
station_characteristics <- station_characteristics %>% relocate(combined_station, lat, lon, sampling_range)

# export station characteristics
tab <-
  station_characteristics %>%
  flextable() %>%
  autofit() %>%
  save_as_docx(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/station_characteristics.docx")

# # plot stations that passed the filtering criteria and export ####
# Only keep the first station name of combined stations
filtered_data$station <- sapply(strsplit(filtered_data$station, ", "), function(x) x[1])

# Generate google map
register_google(key = "AIzaSyCFXkNzarqALqL3qD3Jt5sv-LvPWjLc1mA")

# Plot different filtering steps for supplemental
# Step 1 #####
# Get dataframe of unique combination of station / lat / lon; arrange by latitude and introduce numbering
unique_stations <-
  all_data %>% dplyr::select(combined_station, lat, lon, probability) %>% unique()

# Create the additional column
unique_stations <- unique_stations %>%
  group_by(combined_station, lat, lon) %>%
  mutate(has_both = ifelse(any(probability == "present") & any(probability == "absent"), "both", "absent_only")) %>%
  ungroup()

# Combine the data from both plots
combined_data <- all_data %>%
  dplyr::select(lat, lon, combined_station, species, cells_L) %>%
  distinct() %>%
  drop_na(cells_L) %>%
  group_by(combined_station, species) %>%
  dplyr::summarize(lat = mean(lat), lon = mean(lon), .groups = 'drop') 


pacman::p_load(leaflet, htmlwidgets, leaflet.extras)

# Define a custom rainbow color palette for species
custom_palette <- rainbow(n = length(unique(combined_data$species)))

# Create a color factor based on the rainbow palette
species_palette <- colorFactor(custom_palette, combined_data$species)

# Create a leaflet map
map <- leaflet() %>%
  addProviderTiles(providers$Esri.OceanBasemap)

# Add circle markers for each species
unique_species <- unique(combined_data$species)
for (species in unique_species) {
  species_data <- combined_data %>% filter(species == !!species)
  map <- map %>%
    addCircleMarkers(
      data = species_data,
      ~ lon,
      ~ lat,
      color = ~ species_palette(species),
      radius = 3,  # Increase the radius for better visibility
      fill = TRUE,
      fillOpacity = 0.7,
      group = species
    )
}

# Add layer control to toggle species
map <- map %>%
  addLayersControl(
    overlayGroups = unique_species,
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = species_palette,
    values = combined_data$species,
    title = "Species"
  ) %>%
  setMaxBounds(-2, 50, 30, 73)

# Save the leaflet map as an HTML file to the specified location
saveWidget(map, file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/interactive_map.html")

# Get dataframe of unique combination of station / lat / lon; arrange by latitude and introduce numbering
unique_stations <-
  filtered_data  %>% dplyr::select(station, lat, lon, probability) %>% unique() 

unique_stations <- unique_stations %>%
  arrange(station, desc(probability)) %>%
  group_by(station) %>%
  slice_head(n = 1) %>%
  ungroup() %>% 
  arrange(desc(lat))
  
unique_stations <-
  unique_stations %>% arrange(desc(lat)) %>%
  mutate(station_number = rep(1:nrow(unique_stations)))

pacman::p_load(ggrepel, ggOceanMaps, ggspatial)

coordinates <- data.frame(lon = c(3.5, 3.5, 23, 23), lat = c(54, 73.5, 73.5, 54))

stations_meeting_criteria <- basemap(data = coordinates, bathymetry = F, legends = F, land.col = "grey90", rotate = T) +
 geom_spatial_point(
    data = unique_stations,
    aes(x = lon, y = lat, col = probability),
    size = 1.5)  +
  labs(
    title = "Yearly analysed stations",
    x = "Longitude (°)", 
    y = "Latitude (°)"
  ) +
  scale_color_manual(
    values = c("0" = "#E69F00", "1" = "#009E73"),  # Colorblind-friendly orange and green
    labels = c("0" = "absent", "1" = "present")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12, face = "bold", colour = "black"),  # Increase legend text size
    legend.position = c(0.40, 0.85),  # Move legend to top right corner
    legend.justification = c(1, 1)  # Ensure legend is placed correctly within the corner
  ) + 
  # geom_spatial_text_repel(data = unique_stations %>% filter(station_number %in% c(rep(1:11), 23)), aes(x = lon, y = lat, label = station_number),
  #                 point.padding = 0.25,
  #                 color = "darkred",
  #                 size = 3,
  #                 force = 20,
  #                 segment.size = 0.2,
  #                 max.overlaps = 20,
  #                 min.segment.length = 0) +
  ggspatial::annotation_north_arrow(
    location = "tl",  # top left corner
    pad_x = unit(0.1, "in"), 
    pad_y = unit(0.1, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  ggspatial::annotation_scale(
    location = "tr",  # top right corner
    pad_x = unit(0.1, "in"), 
    pad_y = unit(0.1, "in")
  )

coordinates2 <- data.frame(lon = c(7.5, 7.5, 17, 17), lat = c(54, 60, 60, 54))

stations_meeting_criteria2 <- basemap(data = coordinates2, bathymetry = F, legends = F, land.col = "grey85", rotate = T) +
  ggspatial::geom_spatial_point(
    data = unique_stations,
    aes(x = lon, y = lat, col = probability),
    size = 1.5) +
  labs(
    title = "Yearly analysed stations",
    x = "Longitude (°)", 
    y = "Latitude (°)"
  ) +
  scale_color_manual(
    values = c("0" = "#E69F00", "1" = "#009E73"),  # Colorblind-friendly orange and green
    labels = c("0" = "absent", "1" = "present")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12, face = "bold", colour = "white"),  # Increase legend text size
    legend.position = "none",  # Move legend to top right corner
    legend.justification = c(1, 1)  # Ensure legend is placed correctly within the corner
  ) 
# + 
#   geom_spatial_text_repel(data = unique_stations, aes(x = lon, y = lat, label = station_number),
#                    point.padding = 0.25,
#                   color = "darkred",
#                    size = 3,
#                    force = 20,
#                    segment.size = 0.2,
#                    max.overlaps = 20,
#                    min.segment.length = 0)  

P_combined <- ggarrange(stations_meeting_criteria, stations_meeting_criteria2, ncol = 2)

ggsave("stations_meeting_criteria_both_without_label.png",
       P_combined,
       path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
       dpi = 300,
       width = 8,
       height = 5,
       units = "in"
)

# construct 10 year plot with mean cell densities of each species including the amount of present observations #######
mean_cell_dens <- filtered_data %>%
  filter(species != "No Alexandrium" & species != "Alexandrium minutum") %>%
  drop_na(cells_L) %>%
  group_by(station, year2, species, lat, lon) %>%
  drop_na(cells_L) %>%
  dplyr::summarise(
    cells_L = mean(cells_L, na.rm = TRUE), n = n()) %>%
  drop_na(year2)

mean_cell_dens <- full_join(mean_cell_dens, unique_stations) %>% mutate(log_cells_L = log(cells_L))

pacman::p_load(colorspace)

scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Define a list to store the plots
plot_list <- list()

# Specific breaks and their log-transformed values
breaks <- c(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7)
log_breaks <- log(breaks)

# Custom labels for the breaks
custom_labels <- scales::scientific_format()(breaks)

coordinates2 <- data.frame(lon = c(6, 6, 20, 20), lat = c(54, 60, 60, 54))

pacman::p_load(colorspace, viridisLite, ggOceanMaps)

# Create a custom color gradient based on viridis with added orange/red tones
viridis_colors <- viridis(256)
# Modify the upper end of the viridis palette to transition to orange/red
custom_colors <- c(viridis_colors[1:200], colorRampPalette(c("orange", "red"))(56))

# Define custom breaks and sizes
custom_breaks <- c(1, 5, 10, 20, 40)  # Define breaks, starting with 0
custom_sizes <- c(1, 5, 10, 20, 40)    # Corresponding sizes for the breaks

# Loop through each species
for (i in seq_along(c("Alexandrium pseudogonyaulax", "Alexandrium ostenfeldii", "Alexandrium tamarense", "Alexandrium spp."))) {
  species_name <- c("Alexandrium pseudogonyaulax", "Alexandrium ostenfeldii", "Alexandrium tamarense", "Alexandrium spp.")[i]
  
  # Filter the data for the current species
  species_data <- mean_cell_dens %>% filter(species == species_name) %>% drop_na(n)
  
  # Create the plot for the current species
  plot2 <-  basemap(data = coordinates2, bathymetry = F, legends = F, land.col = "grey85", rotate = T) +
    ggspatial::geom_spatial_point(
      data = species_data,
      aes(x = lon, y = lat, col = log_cells_L, size = n)) +
    labs(
      title = "Yearly analysed stations",
      x = "Longitude (°)", 
      y = "Latitude (°)"
    ) +
    facet_grid(
      . ~ year2) +
    scale_color_gradientn(
      colors = viridisLite::plasma(200),
      limits = log(c(10, 10^7)),
      breaks = log_breaks,
      labels = custom_labels
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    scale_size_area(
      limits = c(0, 40),     # Adjusting the limits
      breaks = custom_breaks, # Custom breaks
      max_size = 6,          # Max size for points
      labels = custom_sizes   # Custom sizes for breaks
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      # plot.background = element_rect(fill = "white"),
      legend.title = element_blank(),
      strip.background = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0, vjust = 1, colour = "white", face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) 
  # Condition to remove the legend for specified species
  if (species_name %in% c("Alexandrium pseudogonyaulax", "Alexandrium spp.", "Alexandrium ostenfeldii")) {
    plot2 <- plot2 + theme(legend.position = "none")
  }
  if (species_name %in% c("Alexandrium pseudogonyaulax")) {
    plot2 <- plot2 + theme(strip.text = element_text(size = 14, vjust = 1),
                           plot.margin = unit(c(0, 0, -4, 0), "cm"))
  }
  
  if (species_name %in% c("Alexandrium ostenfeldii", "Alexandrium tamarense")) {
    plot2 <- plot2 + theme(plot.margin = unit(c(-4, 0, -4, 0), "cm"))
  }
  
  # Additional customization for "Alexandrium spp."
  if (species_name == "Alexandrium spp.") {
    plot2 <- plot2 + 
      labs(x = "Longitude", y = "Latitude") +
      theme(
        axis.text.x = element_text(size = 12),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 12),
        axis.ticks.y = element_line(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(-4, 0, 0, 0), "cm")
      ) 
  }
  
  # Add the plot to the list
  plot_list[[species_name]] <- plot2
}

# Combine the plots using patchwork
pacman::p_load(patchwork)

Alex_dens <- wrap_plots(plot_list)  + plot_layout(axis_titles = "collect", ncol = 1)

species_data <- mean_cell_dens %>% filter(species == "Alexandrium pseudogonyaulax") %>% drop_na(n)

breaks <- c(10^1, 10^2, 10^3, 10^4, 10^5)
log_breaks <- log(breaks)
custom_labels <- scales::scientific_format()(breaks)

plot2 <-  basemap(data = coordinates2, bathymetry = F, legends = F, land.col = "grey85", rotate = T) +
  ggspatial::geom_spatial_point(
    data = species_data,
    aes(x = lon, y = lat, col = log_cells_L, size = n)) +
  labs(
    title = "Yearly analysed stations",
    x = "Longitude (°)", 
    y = "Latitude (°)"
  ) +
  # scale_color_manual(
  #   values = c("0" = "#E69F00", "1" = "#009E73"),  # Colorblind-friendly orange and green
  #   labels = c("0" = "absent", "1" = "present")
  # ) +
  facet_wrap(~year2, nrow = 2, ncol = 2) +
  scale_color_gradientn(
    colors = viridisLite::plasma(200),
    limits = log(c(10, 10^5)),
    breaks = log_breaks,
    labels = custom_labels
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  scale_size_area(
    limits = c(0, 40),     # Adjusting the limits
    breaks = custom_breaks, # Custom breaks
    max_size = 6,          # Max size for points
    labels = custom_sizes   # Custom sizes for breaks
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0, vjust = 1, colour = "white", face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 

# construct 10 year plot with mean cell densities of each species including the amount of present observations #######
mean_cell_dens2 <- filtered_data %>%
  filter(species != "No Alexandrium" & species != "Alexandrium minutum") %>%
  drop_na(cells_L) %>%
  group_by(station, year, species, lat, lon) %>%
  drop_na(cells_L) %>%
  dplyr::summarise(
    cells_L = mean(cells_L, na.rm = TRUE), n = n()) %>%
  drop_na(year)

mean_cell_dens3 <- filtered_data %>%
  filter(species != "No Alexandrium" & species != "Alexandrium minutum") %>%
  drop_na(cells_L) %>%
  group_by(station, year, species, lat, lon, month) %>%
  drop_na(cells_L) %>%
  dplyr::summarise(
    cells_L = mean(cells_L, na.rm = TRUE), n = n()) %>%
  drop_na(year)

mean_cell_dens2 <- full_join(mean_cell_dens2, unique_stations) %>% mutate(log_cells_L = log(cells_L))
mean_cell_dens3 <- full_join(mean_cell_dens3, unique_stations) %>% mutate(log_cells_L = log(cells_L))

species_data2 <-
  mean_cell_dens2 %>% filter(species == "Alexandrium pseudogonyaulax" | species == "Alexandrium ostenfeldii") %>% drop_na(n)
species_data3 <-
  mean_cell_dens3 %>% filter(species == "Alexandrium pseudogonyaulax" | species == "Alexandrium ostenfeldii" ) %>% drop_na(n)

custom_breaks2 <- c(2, 4, 6, 8)  # Define breaks, starting with 0
custom_sizes2 <- c(2, 4, 6, 8)    # Corresponding sizes for the breaks

plot3 <-  basemap(data = coordinates2, bathymetry = F, legends = F, land.col = "grey85", rotate = T) +
  ggspatial::geom_spatial_point(
    data = species_data2 %>% filter(year >= 2006 & year <= 2011 & species == "Alexandrium pseudogonyaulax"),
    aes(x = lon, y = lat, col = log_cells_L, size = n)) +
  labs(
    title = "",
    x = "Longitude (°)", 
    y = "Latitude (°)"
  ) +
  # scale_color_manual(
  #   values = c("0" = "#E69F00", "1" = "#009E73"),  # Colorblind-friendly orange and green
  #   labels = c("0" = "absent", "1" = "present")
  # ) +
  facet_wrap(~year) +
  scale_color_gradientn(
    colors = viridisLite::plasma(200),
    limits = log(c(10, 10^5)),
    breaks = log_breaks,
    labels = custom_labels
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  scale_size_area(
    limits = c(0, 7),     # Adjusting the limits
    breaks = custom_breaks2, # Custom breaks
    max_size = 4,          # Max size for points
    labels = custom_sizes2   # Custom sizes for breaks
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0, vjust = 1, colour = "white", face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 
plot3

ggsave(
  "AP_cell_densities2.png",
  plot3,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  dpi = 300,
  width = 6,
  height = 4,
  units = "in"
)

coordinates4 <- data.frame(lon = c(8, 8, 11, 11), lat = c(56.5, 57.5, 57.5, 56.5))
species_data4 <- species_data3 %>% filter(str_detect(station, "VIB") & species == "Alexandrium pseudogonyaulax")

plot4 <-  basemap(data = coordinates4, bathymetry = F, legends = F, land.col = "grey85", rotate = T) +
  ggspatial::geom_spatial_point(
    data = species_data4,
    aes(x = lon, y = lat, col = log_cells_L, size = n)) +
  labs(
    title = "Yearly analysed stations",
    x = "Longitude (°)", 
    y = "Latitude (°)"
  ) +
  # scale_color_manual(
  #   values = c("0" = "#E69F00", "1" = "#009E73"),  # Colorblind-friendly orange and green
  #   labels = c("0" = "absent", "1" = "present")
  # ) +
  facet_wrap(~month + year) +
  scale_color_gradientn(
    colors = viridisLite::plasma(200),
    limits = log(c(10, 10^5)),
    breaks = log_breaks,
    labels = custom_labels
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    # plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0, vjust = 1, colour = "white", face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size = 12),
    axis.ticks.y = element_line(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 14),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 
plot4

ggsave(
  "AP_cell_densities.png",
  plot2,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  dpi = 300,
  width = 6,
  height = 4,
  units = "in"
)

# Save filtered_data as a new txt.file ######
write.table(filtered_data,
            file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/filtered_data.txt",
            sep = "\t",
            row.names = FALSE)

# Load custom functions
source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

install_packages()

# read in filtered_data to start here: #####
filtered_data = read_delim(
  "C:\\Users\\krist\\OneDrive\\Dokumente\\AWI\\Promotion\\Data Analysis\\filtered_data.txt",
  delim = "\t",
  col_names = T
)

filtered_data <- filtered_data %>%
  mutate(limiting_conditions = as.factor(limiting_conditions),
         strat = as.factor(strat))

test <- filtered_data %>% filter(year >= 2000) %>% group_by(station, month) %>%  dplyr::summarise(TP = mean(TP, na.rm = T), TN = mean(TN, na.rm = T),
                                                                               NH4 = mean(NH4, na.rm = T), NO3 = mean(NO3, na.rm = T)) 

ggplot(test, aes(x = month, y = NH4))+
  geom_line() +
  geom_point() +
  facet_wrap(~station, scales = "free")

# apply check_and_fit function that checks for GLM requirements and 
# then performs GLMs to calculate the monthly and yearly probabilities of A. pseudogonyaulax presence
# uses all available data for each station
# increases calculation speed of the following for-loops as non-suiting stations are not even checked in the check_and_fit function
stations_matching_condition <-
  c()  # Initialize an empty vector to store matching stations
unique_stations2 <- unique(filtered_data$station)  # Get unique stations
for (i in seq_along(unique_stations2)) {
  each_station <- unique_stations2[i]
  if (nlevels(factor(
    filtered_data %>% filter(station == each_station) %>% pull(probability)
  )) >= 2) {
    # Append the station to the list of matching stations
    stations_matching_condition <-
      c(stations_matching_condition, each_station)
  }
}

# apply check_and_fit function on stations
# calculates monthly and yearly probability patterns including 95% CI intervals (glm)
# fits logistic regression for yearly probability patterns
# calculates probability patterns within a year (doyly) via GAM and a corresponding fit 

# rename combined station to station first. Old station names matching the combined stations are in "close_stations"
probability_columns <- c("probability", "probability_AO", "probability_AT", "probability_AM", "probability_Aspp")

for (each_station in stations_matching_condition) {
    check_and_fit(filtered_data,
                  each_station,
                  "station",
                  "probparm_station_monthly",
                  "probparm_station_yearly", 
                  "predicted_probs_station_yearly", 
                  "results_list_", 
                  "probparm_station_doyly", 
                  "probparm_station_doyly_fit")
}

# Filter the working directory to include only non-empty dataframes
non_empty_dfs <- lapply(ls(envir = .GlobalEnv, pattern = "^results_list_"), function(list_name) {
  results_list <- get(list_name, envir = .GlobalEnv)
  lapply(results_list, function(df) {
    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      return(df)
    }
    return(NULL)
  })
})

non_empty_dfs <- unlist(non_empty_dfs, recursive = FALSE)

# Export non-empty dataframes to the global environment
list2env(non_empty_dfs, envir = .GlobalEnv)

# Define parameters of interest
parameters_of_interest <-
  c(
    "NH4",
    "NO3",
    "sal",
    "temp",
    "TN",
    "PO4",
    "silicate",
    "wind_ms",
    "DIN",
    "DON",
    "DIP",
    "NO3_PO4",
    "DON_DIN",
    "chl",
    "C_chl",
    "strat",
    "limiting_conditions"
  )

# Initialize an empty vector to store the matching column names
matching_columns <- c()

# Loop through column names of all_data to include lagged parameter columns into parameters_of_interest
for (col_name in colnames(filtered_data)) {
  # Check if col_name contains any parameter in parameters_of_interest using grepl
  if (any(grepl(paste(parameters_of_interest, collapse = "|"), col_name))) {
    matching_columns <- c(matching_columns, col_name)
  }
}

# update parameters_of_interest
parameters_of_interest <- matching_columns

# remove NO3_NO2

parameters_of_interest <- setdiff(parameters_of_interest, "NO3+NO2")

# Prepare data frames to store data generated by the previous check_and_fit functions
Coef_stations <- data.frame() # for coefficients of check_and_fit
sample_size <- data.frame() # for sample size of glms (needed for back-transformation)

source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

# instead of using years_of_presence, rather use years_since_first_observation
# reasoning: after first observation, following years with only absent observations should be included

years_since_first_observation <- find_years_since_first_observation(filtered_data, probability_columns)

# Process each station data
station_parameters_coefficients <- lapply(unique(filtered_data$station), function(each_station) {
  station_data <- filtered_data %>% filter(station == each_station)
  years_station <- years_since_first_observation %>% filter(station == each_station) %>% ungroup() %>% dplyr::select(year, species)
  return(
    process_data(
      station_data,
      years_station,
      parameters_of_interest,
      "station",
      each_station,
      NULL,
      probability_columns
    )
  )
})

# Combine the results into data frames
Coef_stations <- do.call(rbind, lapply(station_parameters_coefficients, function(station_results) {
  do.call(rbind, lapply(station_results, function(probability_results) {
    do.call(rbind, lapply(probability_results, function(parameter_results) {
      parameter_results$coefficients_data
    }))
  }))
}))

sample_size <- do.call(rbind, lapply(station_parameters_coefficients, function(station_results) {
  do.call(rbind, lapply(station_results, function(probability_results) {
    do.call(rbind, lapply(probability_results, function(parameter_results) {
      parameter_results$sample_size_data
    }))
  }))
}))

sample_size <- sample_size %>%
  pivot_longer(
    cols = c(present, absent),
    names_to = "treat",
    values_to = "sample_size"
  )

# Join station coefficients and sample size dataframes for backtransformation of mean and sd
Coef_stations_back <-
  full_join(Coef_stations,
            sample_size) 

Coef_stations_back2 <-
  backtransformation_of_log_transformed_coefficients(unique(Coef_stations_back), parameters_of_interest, "each_station")

# Select columns and pivot wider to get parameter-wise columns
Coef_stations_back2 <- Coef_stations_back2 %>%
  dplyr::select(!n) %>%
  pivot_wider(values_from = data,
              names_from = parameter)

# General data transformation of the coefficient dataframe; introduce significance (< 0.05) column
Coef_stations2 <- Coef_stations_back2 %>%
  pivot_longer(
    cols = starts_with(parameters_of_interest),
    names_to = "parameter",
    values_to = "data"
  ) %>%
  pivot_wider(names_from = type,
              values_from = data) %>%
  mutate(significance = ifelse('p_value' < 0.05, "s", "ns"))

# Introduce confidence intervals for monthly/yearly probability distributions with 1/0 as upper and lower limit, respectively
# Get a list of all probparm dataframes
# Filter dataframes that start with "probparm_"
filtered_dataframes <- grep("^probparm_", ls(), value = TRUE) %>%
  keep(~ !grepl("doy", .))

# Apply the function to the filtered dataframes and update them
updated_dfs <- lapply(mget(filtered_dataframes), calculate_upr_lwr)

# Update the dataframes in the environment
for (i in seq_along(filtered_dataframes)) {
  assign(filtered_dataframes[i], updated_dfs[[i]])
}

# Create an empty list to store the results
probparm_station_doyly_fit_characteristics <- list()

# Find doys at which the probability exceeds 0.05 and drops below
probparm_station_doyly_fit_characteristics <- probparm_station_doyly_fit %>%
    group_by(station) %>%
    dplyr::reframe(
      "t1" = min(doy[data >= 0.1 & !is.na(data)], na.rm = TRUE),
      "t2" = max(doy[data >= 0.1 & !is.na(data)], na.rm = TRUE),
      "p_max" = doy[which.max(data[!is.na(data)])]
    ) %>%
    dplyr::ungroup() 

# Introduce water bodies to the 17 filtered stations
station_and_waterbody <-
  data.frame(station = unique(filtered_data$station)) %>%
  full_join(unique_stations) %>%
  arrange(station_number) %>%
  mutate(
    water_body = c(
      "estuary",
      "estuary",
      "open",
      "estuary",
      "estuary",
      "open",
      "estuary",
      "estuary",
      "coastal",
      "coastal",
      "estuary",
      "coastal",
      "estuary",
      "estuary",
      "estuary",
      "open",
      "estuary",
      "estuary",
      "estuary",
      "estuary",
      "coastal",
      "coastal",
      "open",
      "coastal",
      "estuary",
      "coastal",
      "open",
      "open",
      "estuary",
      "estuary",
      "coastal",
      "coastal",
      "coastal",
      "coastal",
      "coastal",
      "open",
      "estuary",
      "estuary",
      "coastal",
      "coastal",
      "estuary",
      "coastal",
      "coastal",
      "coastal",
      "estuary",
      "coastal",
      "open",
      "open",
      "open",
      "open",
      "open",
      "open",
      "open",
      "coastal"))

# merge with probparm_station_doyly_fit_characteristics_df containing t1, t2, pmax
probparm_station_doyly_fit_characteristics <-
  probparm_station_doyly_fit_characteristics %>%
  full_join(station_and_waterbody) %>%
  filter(!is.infinite(t1))

# Calculate seasonal means
# Load custom functions
source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

seasonal_means <-
  calculate_seasonal_mean(filtered_data, probability_columns)

seasonal_means <- bind_rows(seasonal_means)

seasonal_means2 <-
  pivot_wider(
    aggregate(seasonal_means, data ~ parameter + station, FUN = mean),
    names_from = "parameter",
    values_from = "data"
  )

seasonal_means2$probability <- ifelse(is.na(seasonal_means2$probability), 0, seasonal_means2$probability)

seasonal_means2 <- pivot_longer(seasonal_means2, starts_with("probability"), names_to = "species",
                                values_to = "probability",
                                values_drop_na = T)
seasonal_means2 <-  pivot_longer(
  seasonal_means2,
  cols = c(2:13), values_to = "data",
  names_to = "parameter"
)

custom_labeller <- function(variable) {
  # Map your parameters to desired facet titles
  titles <- c("NH4" = "NH<sub>4</sub><sup>+</sup> (&mu;g L<sup>-1</sup>)",
              "NO3" = "NO<sub>3</sub><sup>-</sup> (&mu;g L<sup>-1</sup>)",
              "sal" = "Salinity",
              "silicate" = "Silicate (&mu;mol L<sup>-1</sup>)",
              "temp" = "Temperature (°C)",
              "TN" = "Total nitrogen (&mu;mol L<sup>-1</sup>)",
              "PO4" = "PO<sub>4</sub><sup>3-</sup> (&mu;g L<sup>-1</sup>)",
              "limiting_conditions" = "probability of nutrient<br>limiting conditions",
              "strat" = "probability of<br>stratification",
              "chl" = "Chl-a (&mu;g L<sup>-1</sup>)",
              "wind_ms" = "wind speed (m s<sup>-1</sup>)",
              "C" = "Carbon (&mu;g L<sup>-1</sup>)",
              "DIP" = "DIP (&mu;g L<sup>-1</sup>)",
              "DIN" = "DIN (&mu;g L<sup>-1</sup>)",
              "NO3_PO4" = "DIN:PO<sub>4</sub><sup>3-</sup>",
              "Si_N" = "Si:N",
              "DON" = "DON (&mu;g L<sup>-1</sup>)",
              "DON_DIN" = "DON:DIN",
              "C_chl" = "C:Chl")
  
  # Return the corresponding title for each parameter
  return(titles[variable])
}

pacman::p_load(ggpointdensity)
sa_plot <- seasonal_means2 %>% filter(
  !species %in% c("probability_Aspp", "probability_AM", "probability_AT"))
#   &
#     !parameter %in% c("NO3", "NH4", "PO4", "DON_DIN", "DON", "C", "TN", "DIN", "chl", "DIP")
# )
# sa_plot$parameter <- factor(sa_plot$parameter, levels = c("temp", "sal", "NO3_PO4", "strat", "limiting_conditions"))

seasonal_means_plot <-
    ggplot(sa_plot, aes(x = data, y = probability, col = species)) + 
  geom_pointdensity(size = 1.5) +
  facet_wrap(
    ~ parameter,
    scales = "free",
    labeller = labeller(parameter = custom_labeller),
    strip.position = "bottom",
  ) +
  theme_classic() +
  theme(
    plot.title = element_blank(),  # Remove plot title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    legend.text = element_markdown(),
    axis.title.y = element_markdown(size = 14),
    strip.background = element_blank(),
    axis.text = element_markdown(size = 10),
    strip.text = element_markdown(),
    strip.placement = "outside",
    legend.position = "top"
  ) +
  scale_colour_colorblind(labels = c("<i>A. pseudogonyaulax</i>", "<i>A. ostenfeldii</i>", "<i>A. tamarense</i>")) +
  # stat_smooth(method = "loess", se = F, span = 3) +
  xlab("") +
  ylab("probability of observing<br><i>Alexandrium</i>")  # Split into two lines

seasonal_means_plot

ggsave(
  "fig3_manuscript.png",
  seasonal_means_plot,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  dpi = 300,
  width = 7,
  height = 5,
  units = "in"
)

aov_seasonality <-
  aov(
    data = probparm_station_doyly_fit_characteristics ,
    t1 ~ water_body
  )
summary(aov_seasonality)

TukeyHSD(aov_seasonality)

probparm_station_doyly <-
  full_join(
    probparm_station_doyly,
    probparm_station_doyly_fit_characteristics[, c(1, 5)]
  )
probparm_station_doyly_fit <-
  full_join(
    probparm_station_doyly_fit,
    probparm_station_doyly_fit_characteristics[, c(1, 5)]
  )

# export probparm_station_doyly_fit_characteristics as table
tab <-
  probparm_station_doyly_fit_characteristics %>%
  arrange(water_body) %>%
  flextable() %>%
  autofit() %>%
  save_as_docx(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/probparm_station_doyly_fit_characteristics.docx")

amount_of_observations  <- filtered_data %>%
  filter(
    probability ==  1
  ) %>%
  group_by(station, year) %>%
  dplyr::summarize(n_presence = n()) %>%
  arrange(station) %>%
  mutate(species = "A. pseudogonyaulax")

amount_of_observations2 <- filtered_data %>%
  filter(
    probability_AO ==  1
  ) %>%
  group_by(station, year) %>%
  dplyr::summarize(n_presence = n()) %>%
  arrange(station) %>%
  mutate(species = "A. ostenfeldii")

amount_of_observations_all <- rbind(amount_of_observations, amount_of_observations2)

# Join the datasets, retaining only unique combinations of station, lat, and lon
amount_of_observations_all <- amount_of_observations_all %>%
  left_join(unique_stations %>% dplyr::select(station, station_number), by = "station") 

amount_of_observations_all <-
  amount_of_observations_all %>% arrange(desc(station_number)) 

amount_of_observations_all$station_number <- factor(
  amount_of_observations_all$station_number, levels = unique(amount_of_observations_all$station_number))

amount_of_observations_all$species <-
  factor(
    amount_of_observations_all$species,
    levels = c("A. pseudogonyaulax", "A. ostenfeldii"))


# tab <-
#   amount_of_observations %>%
#   flextable() %>%
#   autofit() %>%
#   save_as_docx(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/amount_of_observations.docx")

custom_labeller <- function(variable) {
  # Map your parameters to desired facet titles
  titles <- c("A. ostenfeldii" = "<b>b)</b>",
              "A. pseudogonyaulax" = "<b>a)</b>")
  
  # Return the corresponding title for each parameter
  return(titles[variable])
}

p <-
  ggplot(
    amount_of_observations_all,
    aes(
      x = year,
      y = station_number,
      fill = n_presence
      # ,
      # colour = factor(probability, c(1, 0))
    )
  ) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 13),
                       breaks = seq(0, 13, by = 3),
                       labels = scales::label_number(accuracy = 1)) +  # Use your preferred color scale
  labs(title = "", x = "Year", y = "Station") +
  theme_classic() +
  facet_wrap(
    ~ species,
    ncol = 1,
    labeller = labeller(species = custom_labeller())) +
  theme(
    plot.title = element_markdown(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    axis.text.y = element_text(),
    strip.background = element_blank(),
    axis.text = element_text(),
    strip.placement = "outside",  # Move strip labels outside the plot area
    strip.text = element_markdown(hjust = 0, face = "bold") 
  ) +
  scale_y_discrete(guide = guide_axis(n.dodge=2))
p

# Save the heatmap as an image file
ggsave(filename = "amount_of_observations.png",
       plot = p,
       path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis",
       units = "in",
       height = 7, width = 4.5)

# plot monthly data of each station and export 
for (each_station in unique(probparm_station_monthly$station)) {
  p.subset <- probparm_station_monthly %>% filter(station == each_station)
  # Print each_station and filename for debugging
  filename <- paste0(each_station, "_CI", ".png")
  create_plot(
    p.subset,
    NULL,
    "station",
    "month",
    paste0(each_station, " monthly"),
    "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/station_monthly",
    filename  # Specify the complete filename here
  )
}

# add station_numbers and water_body to the dataset
filtered_data <-
  filtered_data %>% left_join(station_and_waterbody %>% dplyr::select(station, station_number, water_body))

source("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/All_Monitoring_stations_custom_functions.R")

# Call the function for each probability column
probability_columns <- c("probability", "probability_AO", "probability_AT")

result_all <- do.call(rbind, lapply(probability_columns, function(prob_col) {
  calculate_seasonal_probability(filtered_data %>% filter(station %in% unique_stations$station), prob_col)
}))

result_all$water_body <- factor(result_all$water_body, levels = c("estuary", "coastal", "open"),
                                labels = c("estuaries", "coastal waters", "open waters"))

all_plots_monthly <- list()

  for(i in unique(result_all$water_body)){
   P <- ggplot(result_all %>% filter(water_body == i & species != "probability_AT"), aes(x = month, y = data, group = species)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.75) + 
  ylab("probability of presence") +
  xlab("") +
  ggtitle(i) +
  theme(
    plot.title = element_text(hjust = 0,  vjust = -2, size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    legend.text = element_markdown(size = 8),
    strip.background = element_blank(),
    axis.text = element_text()
  ) + 
  scale_x_discrete(breaks = 1:12,
                          labels = month.abb[1:12],
                          limits = factor(1:12)) + 
  scale_y_continuous(breaks = seq(0.0, 0.3, 0.1),
                     labels = seq(0.0, 0.3, 0.1)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = species),
                alpha = 0.25) +
  coord_cartesian(ylim = c(0, 0.35)) +
     scale_fill_manual(
       values = c("#377eb8", "#ff7f00"),
       labels = c("<i>A. pseudogonyaulax</i>", "<i>A. ostenfeldii</i>"),
       guide = guide_legend(
         nrow = 1, 
         keyheight = unit(0.2, "cm"), 
         keywidth = unit(0.2, "cm")
       )
     )
   all_plots_monthly[[i]] <- P
  }
all_plots <-
  all_plots_monthly[[3]] + all_plots_monthly[[2]] + all_plots_monthly[[1]] + plot_layout(
    axes = "collect",
    ncol  = 1,
    guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "all_stations_monthly.png",
  plot = all_plots,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  units = "in",
  width = 3.0,
  height = 5
)

# plot doyly data of each station and export 
for (each_station in unique(probparm_station_doyly$station)) {
  p.subset <-
    probparm_station_doyly %>% filter(station == each_station) %>%
    mutate(doy = as.numeric(doy))
  p.subset2 <- probparm_station_doyly_fit %>% filter(station == each_station)
  # Print each_station and filename for debugging
  filename <- paste0(each_station, "_CI", ".png")
  create_plot(
    p.subset,
    p.subset2,
    "station",
    "doy",
    paste0(each_station, " doyly"),
    "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/station_doyly",
    paste0(each_station, ".png")  # Specify the complete filename here
  )
}

# plot yearly data of each station and export 
for (each_station in unique(probparm_station_yearly$station)) {
  p.subset <- probparm_station_yearly %>% filter(station == each_station) %>% distinct()
  p.subset2 <- predicted_probs_station_yearly %>% filter(station == each_station)
  filename <- paste0(each_station, "_CI", ".png")
  create_plot(
    p.subset,
    p.subset2,
    "station",
    "year",
    paste0(each_station),
    "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/station_yearly2",
    filename
  )
}

# Modelling the probability of presence as a function of temperature for each station

# Filter data and remove NAs for temperature
filtered_data_temp <- filtered_data %>% drop_na(temp) 

# Initialize an empty list to store the results for each station
results_list <- list()

# Loop through each station in the data
stations <- unique(filtered_data$station)
for (each_station in stations) {
  # Filter the data for the current station
  test <- filtered_data_temp %>% filter(station == each_station & month <= 10 & month >= 5)
  
  if(nrow(test) >= 50) {

  # Fit the GAM model
  M1a_gam <- gam(
    data = test,
    probability ~ s(temp, bs = "cp", k = 3),
    family = binomial
  )
  if(summary(M1a_gam)$s.pv < 0.1){
  # Create a new data frame with a range of temp values
  new_data <- data.frame(temp = seq(10, 22.5, length.out = 1000))

  # Use predict function to get the fitted values (modelled probabilities) and confidence intervals
  pred <- predict(M1a_gam, new_data, type = "response", se.fit = TRUE)
  
  # Add the predictions and confidence intervals to the new_data dataframe
  new_data$probability <- pred$fit
  new_data$se <- pred$se.fit
  new_data$lower <- new_data$probability - 1.96 * new_data$se
  new_data$upper <- new_data$probability + 1.96 * new_data$se
  new_data$station <- each_station
  
  # Add the result to the list
  results_list[[each_station]] <- new_data
}}}

# Combine the results into a single dataframe
results_df <- bind_rows(results_list)

# Plot the results for each station
P1 <- ggplot(results_df, aes(x = temp, y = probability, color = station)) +
  geom_line() +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(
    title = "",
    x = "temperature (°C)",
    y = "probability"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

results_df2 <- results_df 

results_df2 <- full_join(results_df2, station_and_waterbody[,c(1:3, 5:6)])

P2 <- ggplot(results_df2, aes(x = temp, y = probability, color = water_body, group = station)) +
  geom_line() +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(
    title = "",
    x = "temperature (°C)",
    y = "probability"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

P3 <- ggplot(
  results_df2 %>% group_by(water_body, temp) %>% dplyr::summarise(
    probability = mean(probability, na.rm = T)
  ),
  aes(x = temp, y = probability, color = water_body)
) +
  geom_line() +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(title = "",
       x = "temperature (°C)",
       y = "probability") +
  theme_minimal()+
  theme(legend.position = "none")

temp_combined <- P1 + P2 + P3 + plot_layout(ncol = 3)

ggsave(
  filename = "temp_combined.png",
  plot = temp_combined,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  units = "in",
  width = 7,
  height = 5
)

P3 <- ggplot(
  results_df2 %>% group_by(water_body, temp) %>% dplyr::summarise(
    probability = mean(probability, na.rm = T)
  ),
  aes(x = temp, y = probability, color = water_body)
) +
  geom_line(linewidth = 1) +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(title = "",
       x = "temperature (°C)",
       y = "probability of observing<br><i>A. pseudogonyaulax</i>") +
  theme_classic() +
  theme(
    plot.title = element_blank(),  # Remove plot title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    legend.text = element_markdown(),
    axis.title.y = element_markdown(size = 12),
    strip.background = element_blank(),
    axis.text = element_markdown(size = 10),
    strip.text = element_markdown(),
    strip.placement = "outside",
    legend.position = "top") +
  scale_color_colorblind(labels = c("coastal waters", "estuaries", "open waters")) +
  scale_x_continuous(breaks = seq(12, 22, 2))
  

P3
ggsave(
  filename = "temp_Ap.png",
  plot = P3,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  units = "in",
  width = 4,
  height = 3.5
)

Coef_stations2 <-
  full_join(Coef_stations2, unique_stations)

# ##### Plot means of parameters when A.p. is absent or present --> Coef_stations #####
dodge <- position_dodge(width = 0.9)
limits <- aes(ymin = mean, ymax = mean + sd)

# # Save each facet as a separate image file ####
facet_names <- unique(Coef_stations2 %>% dplyr::select(parameter) %>% drop_na(parameter))

Coef_stations_sub <-
  Coef_stations2 %>% filter(species %in% c("probability", "probability_AO")) %>% drop_na(mean) 

Coef_stations_sub <- Coef_stations_sub %>% arrange(station_number) %>% convert_as_factor(station_number) %>% drop_na(station_number)

custom_labeller <- function(variable) {
  # Map your parameters to desired facet titles
  titles <- c("probability" = "<b>a)</b>",
              "probability_AO" = "<b>b)</b>")
  
  # Return the corresponding title for each parameter
  return(titles[variable])
}

for (facet_name in facet_names$parameter) {
  p_subset <- Coef_stations_sub %>% filter(parameter == facet_name)
  
  if(nrow(p_subset) > 2) {
  ggsave(
    filename = paste0(facet_name, ".png"),
    
    plot = ggplot(
      p_subset,
      aes(
        x = station_number,
        y = mean,
        fill = treat,
        label = significance
      )
    ) +
      scale_fill_colorblind() +
      geom_bar(stat = "identity", position = dodge) +
      geom_linerange(limits, position = dodge, linewidth = 0.25) +
      ylab(facet_name) +
      ggtitle("") +
      xlab("station") +
      facet_wrap(~ species, scales = "free_y",
      labeller = labeller(species = custom_labeller)) +
      theme(
        plot.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_markdown()
      ),
    height = 2.5,
    width = 4.5,
    unit = "in",
    path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/abiotic_parameters")
} }

pacman::p_load(patchwork)


custom_labeller <- function(variable) {
  # Map your parameters to desired facet titles
  titles <- c("NH4" = "NH<sub>4</sub><sup>+</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
              "NO3" = "NO<sub>3</sub><sup>-</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
              "sal" = "Salinity",
              "silicate" = "Silicate (<i>&mu;</i>mol L<sup>-1</sup>)",
              "temp" = "Temperature (°C)",
              "TN" = "Total nitrogen (<i>&mu;</i>mol L<sup>-1</sup>)",
              "PO4" = "PO<sub>4</sub><sup>3-</sup> (<i>&mu;</i>mol L<sup>-1</sup>)",
              "limiting_conditions" = "probability of nutrient<br>limiting conditions",
              "strat" = "probability of<br>stratification",
              "chl" = "Chl-a (<i>&mu;</i>mol L<sup>-1</sup>)",
              "wind_ms" = "wind speed (m s<sup>-1</sup>)",
              "C" = "Carbon (<i>&mu;</i>mol L<sup>-1</sup>)",
              "DIP" = "DIP (<i>&mu;</i>mol L<sup>-1</sup>)",
              "DIN" = "DIN (<i>&mu;</i>mol L<sup>-1</sup>)",
              "NO3_PO4" = "DIN:PO<sub>4</sub><sup>3-</sup>",
              "Si_N" = "Si:N",
              "DON" = "DON (<i>&mu;</i>mol L<sup>-1</sup>)",
              "DON_DIN" = "DON:DIN",
              "C_chl" = "C:Chl")
  
  # Return the corresponding title for each parameter
  return(titles[variable])
}

create_plot_abiotic <- function(facet_name) {
  
  # Use the custom labeller function to get the y-axis label
  y_label <- custom_labeller(facet_name)
  
  p_subset <- Coef_stations_sub %>% 
    filter(parameter == facet_name, species == "probability") %>% 
    drop_na(mean)
  ggplot(
    p_subset,
    aes(
      x = station_number,
      y = mean,
      fill = treat,
      label = significance
    )
  ) +
    scale_fill_colorblind() +
    geom_bar(stat = "identity", position = dodge) +
    geom_linerange(limits, position = dodge, linewidth = 0.5) +
    ylab(y_label) +
    # geom_text(position = dodge, vjust = 1, size = 2) +
    ggtitle("") +
    xlab("station") +
    theme(
      plot.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 90),
      axis.title.y = element_markdown(),
      strip.text = element_markdown(),
      strip.background = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      strip.text.x = element_markdown(hjust = 0, margin=margin(l=0))# Adjust overall margins
    )
}

# Save each facet plot
facet_names <- c("NO3", "PO4", "temp", "C_chl")

plots <- lapply(facet_names, create_plot_abiotic)
pacman::p_load(patchwork)

facet_plots_grid <-
  plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plot_layout(
    ncol = 2,
    nrow = 2,
    axis_titles = "collect",
    guides = "collect"
  ) &
  theme(legend.position = "bottom")
# Save the grid of plots
ggsave(
  filename = "facet_plots_grid.png",
  plot = facet_plots_grid,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/abiotic_parameters",
  units = "in",
  width = 6,
  height = 6
)

# Create method figure as part of figure 2 in the manuscript: 
# Function to calculate Gaussian distribution
gaussian_distribution <- function(x, mean, sd) {
  exp(-(x - mean)^2 / (2 * sd^2)) / (sd * sqrt(2 * pi))
}

# Generate data
x_values <- seq(-5, 5, length.out = 1000)
y_values <- gaussian_distribution(x_values, mean = 0, sd = 1)

# Scale y-values to ensure the maximum is 1
y_values <- y_values / max(y_values)

# Find t1 and t2
t1 <- min(x_values[y_values > 0.1])
t2 <- max(x_values[y_values > 0.1])

# Find tmax and pmax
tmax <- x_values[which.max(y_values)]
pmax <- max(y_values)

# Create data frame for plotting
df <- data.frame(x = x_values, y = y_values)

plot_list <- list()

manuscript_subset <- data.frame()

probparm_station_doyly <- left_join(probparm_station_doyly, station_and_waterbody)
probparm_station_doyly_fit <- left_join(probparm_station_doyly_fit, station_and_waterbody)

# plot doyly data of each water body and export 
for (each_water_body in unique(probparm_station_doyly$water_body)) {
  p.subset <-
    probparm_station_doyly %>% filter(water_body == each_water_body & station != "STRETUDDEN") %>%
    mutate(doy = as.numeric(doy)) %>% drop_na(water_body)
  p.subset2 <-
    probparm_station_doyly_fit %>% filter(water_body == each_water_body) %>% drop_na(water_body)
  plot <- create_plot(
    p.subset,
    p.subset2,
    "station",
    "doy",
    paste0(each_water_body, " doyly"),
    "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/Prob/Seasonal_variation/water_body_doyly",
    paste0(each_water_body, ".png")  # Specify the complete filename here
  )
  # Merge the current list of plots with the overall list
  plot_list <- c(plot_list, plot)
}

label = "D = t<sub>2</sub> - t<sub>1</sub>"
label2 = "p<sub>max</sub>"

# Plot
doy_method <- ggplot(df, aes(x, y)) +
  geom_line() +
  geom_segment(
    y = 0.1,
    yend = 0.1,
    linetype = "dashed",
    x = t1,
    xend = t2
  ) +
  geom_segment(
    y = 0,
    yend = 0.1,
    linetype = "dashed",
    x = t1,
    xend = t1
  ) +
  geom_segment(
    y = 0,
    yend = 0.1,
    linetype = "dashed",
    x = t2,
    xend = t2
  ) +
  geom_richtext(
    x = tmax,
    y = 0.13,
    label = label,
    fill = NA,
    label.color = NA,
    size = 4
  ) +
  geom_richtext(
    x = tmax,
    y = 1.04,
    label = label2,
    fill = NA,
    label.color = NA,
    size = 4
  ) +
  scale_x_continuous(breaks = c(tmax, t1, t2),
                     labels = c("t<sub>max</sub>", "t<sub>1</sub>", "t<sub>2</sub>")) +
  scale_y_continuous(breaks = c(0.1, seq(0, 1, by = 0.5)),
                     limits = c(0, 1.05)) +  # Add breaks for y-axis
  labs(x = "time (doy)", y = "probability") +
  ggtitle("a)") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = 'white'),
    legend.title = element_blank(),
    axis.text.x = element_markdown(),
    strip.background = element_blank(),
    axis.title = element_text(),
    axis.text = element_text(),
    plot.title = element_text()
  ) 


fig2_manuscript <-
  ggarrange(
    doy_method + ggtitle("a)") + theme(plot.title = element_text(
      hjust = 0, size = 12, face = "bold"
    )),
    plot_list$`estuary doyly` + ggtitle("b)") + ylim(c(0, 0.65)) + xlab("time (doy)") + theme(plot.title = element_text(
      hjust = 0, size = 12, face = "bold"
    )),
    plot_list$`open doyly` + ggtitle("c)") + ylim(c(0, 0.65)) + xlab("time (doy)") + theme(plot.title = element_text(
      hjust = 0, size = 12, face = "bold"
    )),
    plot_list$`coastal doyly` + ggtitle("d)") + ylim(c(0, 0.65)) + xlab("time (doy)") + theme(plot.title = element_text(
      hjust = 0, size = 12, face = "bold"
    )),
    ncol = 2,
    nrow = 2
  )

ggsave(
  "fig2_manuscript.png",
  fig2_manuscript,
  path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/",
  dpi = 300,
  width = 6,
  height = 6,
  units = "in"
)

manuscript_stations <-
  probparm_station_yearly %>% filter(
    station %in% c(
      "ARH170006",
      "VIB3708",
      "NOR409",
      "ANHOLT E",
      "Heiligendamm",
      "KBH431",
      "Pricken",
      "VEJ0006870"
    )
  )

manuscript_stations <- calculate_upr_lwr(manuscript_stations) %>% left_join(unique_stations)

# Create a list to store the individual plots
plots <- list()  
plots2 <- list()

# Loop to create and store the individual plots of each station in manuscript_stations subset 
for (each_station in unique(manuscript_stations$station)) {
  manuscript_subset <-
    manuscript_stations %>% filter(station == each_station)
  p.subset2 <-
    predicted_probs_station_yearly %>% filter(station == each_station)
  
  station_number <- manuscript_subset$station_number[1]  # Get the station_number for the current station
  
  gg <- create_plot3(
    manuscript_subset,
    p.subset2,
    "station",
    "year",
    station_number,  # Pass the station_number value
    each_station
  )
  
  # Store the individual plot in the list
  plots[[length(plots) + 1]] <- ggplotGrob(gg)
  plots2[[length(plots2) + 1]] <- gg
  
 }
pacman::p_load(patchwork)
Fig_manuscript <- plots2[[8]] + plots2[[5]] + plots2[[2]] + plots2[[7]] + plot_layout(ncol = 1, axes = "collect")
Fig_manuscript2 <- plots2[[6]] + plots2[[1]] + plots2[[4]] + plots2[[3]] + plot_layout(ncol = 1, axes = "collect")

ggsave("Fig_prediction_1b.png",
       Fig_manuscript,
       path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/manuscript",
       dpi = 300,
       width = 4,
       height = 8,
       units = "in"
)
ggsave("Fig_prediction_2b.png",
       Fig_manuscript2,
       path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/manuscript",
       dpi = 300,
       width = 4,
       height = 8,
       units = "in"
)

  ## BUBBLE PLOTS ############
  
  ## find years of presence alternatively to reduce standard deviations
  years_of_presence_station <- data.frame(station = " ", year = " ")
  for(object in unique(all_data$station)){
    years_of_presence_df <- all_data[all_data$station == object,] %>% dplyr::select(year, station, probability) %>% 
      filter(probability == "present") %>% as.data.frame(station = station, year = year) %>% dplyr::select(station, year) %>% distinct()
    years_of_presence_station <- rbind(years_of_presence_station, years_of_presence_df)
  }
  
  years_of_presence_station$year <- as.numeric(years_of_presence_station$year)
  
  register_google(key = "AIzaSyCFXkNzarqALqL3qD3Jt5sv-LvPWjLc1mA")
  colnames(all_data)[colnames(all_data) %in% c("lon", "lat")] <-
    c("longitude", "latitude")
  colnames(filtered_data)[colnames(filtered_data) %in% c("lon", "lat")] <-
    c("longitude", "latitude")
  
  sq_map3 <-
    get_map(
      location = c(9, 54, 16, 60),
      maptype = "satellite",
      source = "google",
      color = "color",
      darken = 0.6,
      zoom = 6
    )
  
  # Define your custom color palette
  colors <- c(colorblind_pal()(4))[c(4, 2, 3, 1)]
  custom_colors <-
    c(
      "Alexandrium" = colors[1],
      "Alexandrium pseudogonyaulax" = colors[3],
      "Alexandrium ostenfeldii" = colors[2]
    )
  custom_colors2 <-
    c(
      "Alexandrium" = colors[1],
      "Alexandrium pseudogonyaulax" = colors[3],
      "Alexandrium ostenfeldii" = colors[2],
      "Prorocentrum" = colors[4]
    )
  
  # Define custom color labels with italics
  custom_color_labels <- c(
    "Alexandrium" = expression(italic("Alexandrium")),
    "Alexandrium pseudogonyaulax" = expression(italic("Alexandrium pseudogonyaulax")),
    "Alexandrium ostenfeldii" = expression(italic("Alexandrium ostenfeldii"))
  )
  
  
# with first year of appearance 
first_year_of_appearance <- all_data %>% group_by(station) %>% filter(probability == "present") %>% aggregate(year ~ station, FUN = "min")

first_year <- all_data %>%  group_by(probability, station, longitude, latitude, region) %>%
  filter(probability == "present" & month > 5 & month < 10) %>% dplyr::summarise(data = mean(cells_L))
first_year <- full_join(first_year, first_year_of_appearance)

# Create a loop to generate maps for each year
generate_map <- function(yr) {
  current_stations  <- first_year[first_year$year < yr, ]
  current_stations2  <- first_year[first_year$year == yr, ]
  
  p <- ggmap(sq_map3) +
    geom_point(data=current_stations , aes(x = longitude, y = latitude, col = data), size = 3, show.legend = T) +
    geom_point(data=current_stations2 , aes(x = longitude, y = latitude, col = data), size = 5, shape = 18, show.legend = T) +
    labs(title = "First appearance including mean cell counts June-September",
         x = "Longitude", y = "Latitude",
         col = "years") +
    theme_minimal() +
    scale_color_viridis_c(breaks = seq(0, 5000, 1000), limits = c(0, 5000)) +
    theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'), legend.title = element_blank(), strip.background = element_blank()) +
    annotate("text", x = 8, y = 53.5, label = yr, col= "red", size = 10)
  
  # Save the plot as an image (e.g., PNG)
  filename <- paste0("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/map_", yr, ".png")
  ggsave(filename, p)
}

# Generate and save maps for each year
years <- seq(1997, 2023, 1)
lapply(years, generate_map)


png_files <- list.files(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/", full.names = TRUE)

# Read the PNG images using magick and save as gif
animated_gif <- image_read(png_files) %>%
  image_animate(fps = 0.5) %>% save_animation(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/station_animation.gif")


# with all years of appearance and increasing point size relative to the amount of observations per year

all_years <- data.frame()
for(each_station in unique(all_data$station)){
  for(each_year in unique(all_data %>% filter(station == each_station & probability == "present" ) %>% pull(year))){
    all_years2 <- all_data %>% group_by(station, latitude, longitude, year, probability) %>% filter(station == each_station & probability == "present" & year == each_year & month < 10 & month > 5)  %>%
      dplyr::summarise(present_appearances = n(), cells_L = mean(cells_L)) %>% mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))
    all_years <- rbind(all_years, all_years2)
    
  }
}

all_years2 <- data.frame()
for(each_station in unique(all_data$station)){
  for(each_year in unique(all_data %>% filter(station == each_station & probability == "absent" & year > 1996) %>% pull(year))){
    all_years3 <- all_data %>% group_by(station, latitude, longitude, year, probability) %>% filter(station == each_station & probability == "absent" & year == each_year & month < 10 & month > 5)  %>%
      dplyr::summarise(data = n()) %>% mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))
    all_years2 <- rbind(all_years2, all_years3)
    
  }
}

colnames(all_years2) <- c("station", "latitude", "longitude", "year", "probability2", "data2")

test <- inner_join(all_years, all_years2)
test$ratio <- test$data / test$data2

# Create a loop to generate maps for each year
generate_map <- function(yr) {
  current_stations  <- all_years[all_years$year < yr, ]
  current_stations2  <- all_years[all_years$year == yr, ]
  
  p <- ggmap(sq_map3) +
    geom_point(data=current_stations , aes(x = longitude, y = latitude, col = data), size = 3, show.legend = T) +
    geom_point(data=current_stations2 , aes(x = longitude, y = latitude, col = data), size = 5, shape = 18, show.legend = T) +
    labs(title = "Year wise occurences with amount of observations - only 06-09",
         x = "Longitude", y = "Latitude",
         col = "years") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'), legend.title = element_blank()) +
    scale_color_viridis_c(breaks = seq(0, 10, 2), limits = c(0, 10)) +
    annotate("text", x = 8, y = 53.5, label = yr, col= "red", size = 10)
  
  # Save the plot as an image (e.g., PNG)
  filename <- paste0("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/map_", yr, ".png")
  ggsave(filename, p)
}

# Generate and save maps for each year
years <- seq(1997, 2023, 1)
lapply(years, generate_map)


png_files <- list.files(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/", full.names = TRUE)

# Read the PNG images using magick and save as gif 
animated_gif <- image_read(png_files) %>%
  image_animate(fps = 0.5) %>% save_animation(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif2/station_animation2.gif")


## with mean cell counts library(ggmap)

generate_map <- function(yr) {
  current_stations <- all_years[all_years$year == yr, ]
  
  p <- ggmap(sq_map3) +
    geom_point(data = current_stations, aes(x = longitude, y = latitude, col = cells_L), size = 3, show.legend = TRUE) +
    labs(title = "Year wise occurrences with mean cells/L - only 06-09",
         x = "Longitude", y = "Latitude",
         col = "years") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          legend.title = element_blank()) +
    scale_color_viridis_c(breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
    scale_fill_viridis_c(breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
    annotate("text", x = 8, y = 53.5, label = yr, col = "red", size = 10) +
    guides(fill = "none")  # Remove the legend for the gradient background
  
  # Save the plot as an image (e.g., PNG)
  filename <- paste0("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/map_", yr, ".png")
  ggsave(filename, p)
}

# Generate and save maps for each year
years <- seq(1997, 2021, 1)
lapply(years, generate_map)

png_files <- list.files(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/", full.names = TRUE)

# Read the PNG images using magick and save as gif
animated_gif <- image_read(png_files) %>%
  image_animate(fps = 2) %>%
  save_animation(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif2/station_animation_cellsL.gif")


## with mean cell counts library(ggmap)

generate_map <- function(yr) {
  current_stations <- all_years[all_years$year == yr, ]
  
  p <- ggmap(sq_map3) +
    geom_point(data = current_stations, aes(x = longitude, y = latitude, col = cells_L), position = position_jitter(width = 0.1, height = 0.1, seed = 100), size = 3, show.legend = TRUE) +
    labs(title = expression(paste("Expansion of", italic(" A. pseudogonyaulax"))),
         x = "Longitude", y = "Latitude",
         col = "years") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          legend.title = element_blank()) +
    scale_color_viridis_c(breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
    scale_fill_viridis_c(breaks = seq(0, 35000, 5000), limits = c(0, 35000)) +
    annotate("text", x = 8, y = 53.5, label = yr, col = "red", size = 10) +
    guides(fill = "none")  # Remove the legend for the gradient background
  
  # Save the plot as an image (e.g., PNG)
  filename <- paste0("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/map_", yr, ".png")
  ggsave(filename, p)
}

# Generate and save maps for each year
years <- seq(1997, 2021, 1)
lapply(years, generate_map)

png_files <- list.files(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/", full.names = TRUE)

# Read the PNG images using magick and save as gif
animated_gif <- image_read(png_files) %>%
  image_animate(fps = 2) %>%
  save_animation(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif2/station_animation_cellsL2.gif")


## With ratio present/absent and over the years

# Create a loop to generate maps for each year
generate_map <- function(yr) {
  current_stations  <- test[test$year < yr, ]
  current_stations2  <- test[test$year == yr, ]
  
  p <- ggmap(sq_map3) +
    geom_point(data=current_stations , aes(x = longitude, y = latitude, col = ratio), size = 3, show.legend = T) +
    geom_point(data=current_stations2 , aes(x = longitude, y = latitude, col = ratio), size = 5, shape = 18, show.legend = T) +
    labs(title = "Year wise occurences with ratio of observations present/absent - only 06-09",
         x = "Longitude", y = "Latitude",
         col = "years") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'), legend.title = element_blank()) +
    scale_color_viridis_c(breaks = seq(0, 5, 1), limits = c(0, 5)) +
    annotate("text", x = 8, y = 53.5, label = yr, col= "red", size = 10)
  
  # Save the plot as an image (e.g., PNG)
  filename <- paste0("C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/map_", yr, ".png")
  ggsave(filename, p)
}

# Generate and save maps for each year
years <- seq(1997, 2023, 1)
lapply(years, generate_map)


png_files <- list.files(path = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif/", full.names = TRUE)

# Read the PNG images using magick and save as gif 
animated_gif <- image_read(png_files) %>%
  image_animate(fps = 0.5) %>% save_animation(file = "C:/Users/krist/OneDrive/Dokumente/AWI/Promotion/Data Analysis/gif2/station_animation3.gif")


