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

The analysis runs on national monitoring data. None of it is included in this
repository (nor ever was), and the Norwegian and Danish data may not be
redistributed.

| Country | Source | Access |
|---|---|---|
| Norway | Norwegian national phytoplankton and hydrography monitoring | Not redistributable; from the data holder on request (see the paper's data availability statement) |
| Denmark | Danish national monitoring (NOVANA) | Not redistributable; from the data holder on request (see the paper's data availability statement) |
| Sweden | SHARK, Swedish Meteorological and Hydrological Institute (SMHI) | Public: [SHARKweb](https://sharkweb.smhi.se) |
| Germany | Leibniz Institute for Baltic Sea Research (IOW) | Public: [ODIN2](https://odin2.io-warnemuende.de); download dates are in the file names |

The scripts expect one folder per country in the repository root:

```
norway/    phyto_aug24.txt, nutrients.txt, norway_ctd_full.txt
denmark/   DK_counts.txt, DK_CTD.txt, DK_water_quality.txt, DK_secci_kd.txt,
           wind_data.txt, stationen_wind.txt
sweden/    sharkweb_phyto_new.txt, sharkweb_phyto_2024.txt,
           sharkweb_phys2.txt, sharkweb_phys_2024.txt
germany/   odin2_2024-01-31_095801.txt, odin2_2025-09-12_070935_red.txt,
           Station_details.txt (station metadata)
```

`.gitignore` excludes these folders and every intermediate and output file,
so that no monitoring data can be committed by accident.

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

Open the `time-series-expansion.Rproj` file in RStudio. This sets the
working directory to the repository root, which is what the scripts'
relative paths assume. Then run any script top to bottom.

```r
# R >= 4.2, inside RStudio
source("Time_series_analysis_custom_functions.R")
install_packages()   # pacman::p_load, installs what is missing
source("Time_series_analysis.R")
```

Requires `data.table`, `dbscan`, `dplyr`, `flextable`, `fishmethods`,
`fuzzyjoin`, `furrr`, `geosphere`, `ggOceanMaps`, `ggplot2`, `ggpubr`, `ggtext`,
`ggthemes`, `ggspatial`, `here`, `mgcv`, `officer`, `patchwork`, `purrr`,
`readr`, `readxl`, `rstatix`, `scales`, `stringr`, `tidyverse`,
`viridisLite`. `install_packages()` handles all of them.

## Known limitations

- The working directory is resolved via `here::here()`, so it no longer
  requires RStudio to locate the project root. `Time_series_analysis.R` still
  calls `.rs.restartR()` twice between pipeline stages to free memory — an
  RStudio-only function — so running the full script end-to-end still
  requires RStudio; to run a stage headless, replace that call with
  `rm(list = ls()); gc()` or run each section as a separate `Rscript`
  invocation.
- The main script is a linear pipeline rather than a set of functions. Run it
  top to bottom or resume from a written intermediate.
- Station code assignments in `web_export.R` are duplicated from the main
  script and must be kept in sync manually.
- Random seeds are fixed for the station-level parameter estimates
  (`set.seed(1234)`) and for the bootstrap (10,000 resamples) of the
  season's start, end and peak probability (`set.seed(123)`). DBSCAN
  clustering is deterministic. Package versions are not recorded, so results
  can differ slightly with other versions.

## Status

Archived on publication. Kept for transparency and reuse of the harmonisation
and modelling code. Issues are not monitored. For questions about the study,
contact the corresponding author.

## License

Code: MIT. Monitoring data belong to the respective national programmes and are
subject to their own terms.
