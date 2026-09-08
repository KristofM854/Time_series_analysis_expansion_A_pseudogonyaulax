# Time series analysis of *Alexandrium pseudogonyaulax* in Northern European waters

R code accompanying Möller et al. (2026), *Harmful Algae*. Archived, not maintained.

Analysis of four national phytoplankton monitoring programmes (Norway, Denmark,
Sweden, Germany) to test whether *Alexandrium pseudogonyaulax* has expanded its
range and seasonal window across Northern European waters.

## Citation

> Möller, K., Carstensen, J., Jakobsen, H., Engesmo, A., Karlson, B. (2026).
> Time series analysis of the toxic dinoflagellate *Alexandrium pseudogonyaulax*
> across Northern European waters. *Harmful Algae*.
> https://doi.org/10.1016/j.hal.2026.103060

## Data

The analysis runs on national monitoring data that are not redistributed here.
Obtain them from the source programmes:

| Country | Source | Files expected by the scripts |
|---|---|---|
| Norway | Norwegian phytoplankton and hydrography monitoring | `norway/phyto_aug24.txt`, `norway/nutrients.txt`, `norway/norway_ctd_full.txt` |
| Denmark | ODA / Danish national monitoring (NOVANA) | `DK_counts.txt`, `DK_CTD.txt`, `DK_water_quality.txt`, `DK_secci_kd.txt`, `wind_data.txt`, `stationen_wind.txt` |
| Sweden | SHARKweb (SMHI) | `sharkweb_phyto_new.txt`, `sharkweb_phyto_2024.txt`, `sharkweb_phys2.txt`, `sharkweb_phys_2024.txt` |
| Germany | IOW ODIN2 | `odin2_2024-01-31_095801.txt`, `odin2_2025-09-12_070935_red.txt` |

Place them in the repository root, with the Norwegian files in a `norway/`
subdirectory. `Station_details.txt` holds station metadata.

## Contents

| File | Purpose |
|---|---|
| `Time_series_analysis.R` | Main pipeline. Loads and harmonises the four monitoring datasets, joins abiotic and phytoplankton records within a one-day window, merges programmes, then fits the models and builds the manuscript figures. |
| `Time_series_analysis_custom_functions.R` | Package loading and ~35 helper functions (seawater density, stratification index, station clustering, model fitting and plotting wrappers). Sourced by the other two scripts. |
| `web_export.R` | Optional. Writes `stations.json`, `yearly_probability.json` and `meta.json` for a web front end. Run only after the main script has produced `all_data.txt` and `filtered_data.txt`. |

The main script is organised into commented sections. The intermediate products
`norway_combined.txt`, `denmark_combined.txt`, `sweden_combined.txt`,
`germany_combined.txt`, `all_data.txt` and `filtered_data.txt` are written to the
repository root as the pipeline runs. Section `read in all_data to start here`
lets you resume from `all_data.txt` without repeating the harmonisation step.

## Methods implemented

- Station harmonisation across programmes, including DBSCAN clustering of
  nearby stations and a fuzzy join of abiotic and biotic samples on date.
- Water column density and stratification index from temperature and salinity.
- Binomial GLMs of presence/absence against environmental predictors, with
  lagged predictors and harmonic (sine/cosine) seasonal terms.
- Cyclic-spline GAMs (`mgcv`) of occurrence probability against day of year,
  temperature and salinity.
- Seasonal probability surfaces, occurrence heatmaps for
  *A. pseudogonyaulax* and *A. ostenfeldii*, and station maps (`ggOceanMaps`).

## Running it

```r
# R >= 4.2, inside RStudio
source("Time_series_analysis_custom_functions.R")
install_packages()   # pacman::p_load, installs what is missing
source("Time_series_analysis.R")
```

Requires `data.table`, `dbscan`, `dplyr`, `flextable`, `fishmethods`,
`fuzzyjoin`, `furrr`, `geosphere`, `ggOceanMaps`, `ggplot2`, `ggpubr`, `ggtext`,
`ggthemes`, `ggspatial`, `mgcv`, `officer`, `patchwork`, `purrr`, `readr`,
`readxl`, `rstatix`, `rstudioapi`, `scales`, `stringr`, `tidyverse`,
`viridisLite`. `install_packages()` handles all of them.

## Known limitations

- The working directory is set via `rstudioapi::getActiveDocumentContext()`, so
  the scripts only run inside RStudio. Replace with `here::here()` or a manual
  `setwd()` to run headless.
- The main script is a linear pipeline rather than a set of functions. Run it
  top to bottom or resume from a written intermediate.
- Station code assignments in `web_export.R` are duplicated from the main
  script and must be kept in sync manually.
- No random seed is set. DBSCAN clustering is deterministic here, but
  reproduce with the same package versions to be safe.

## Status

Archived on publication. Kept for transparency and reuse of the harmonisation
and modelling code. Issues are not monitored. For questions about the study,
contact the corresponding author.

## License

Code: MIT. Monitoring data belong to the respective national programmes and are
subject to their own terms.
