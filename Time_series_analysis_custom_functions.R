# Time series analysis regarding the expansion of Alexandrium pseudogonyaulax in northern European waters
# by Kristof Möller (Alfred Wegener Institut), Jacob Carstensen and Hans Jakobsen (each Aarhus University)
# as part of the Phd-thesis of Kristof Möller
# file to store all custom functions of All_monitoring_stations R-file


# Create function to install needed packages #####
install_packages <- function() {
  if (!require("pacman"))
    install.packages("pacman")
  
  # Use pacman to load add-on packages as desired
  pacman::p_load(
    psych,
    tidyverse,
    rstatix,
    data.table,
    geosphere,
    stats,
    MuMIn,
    data.table,
    gganimate,
    animation,
    magick,
    ggthemes,
    readr,
    ggplot2,
    dplyr,
    lubridate,
    stringi,
    readr,
    magrittr,
    car,
    fishmethods,
    rnaturalearthdata,
    sf,
    ggspatial,
    ggmap,
    gridExtra,
    gridpattern,
    ggpubr,
    RColorBrewer,
    flextable,
    officer,
    elliptic,
    gam,
    plyr,
    MASS,
    ggcorrplot,
    formatR,
    microbenchmark,
    roxygen2,
    mgcv,
    splines,
    ggtext,
    scales
  )
}

# Function to check and average rows with one day before or after
average_adjacent_days <- function(df, station_col) {
  df <- df %>%
    dplyr::mutate(date = as.Date(date),
                  year = as.numeric(year),
                  month = as.numeric(month),
                  day = as.numeric(day)) %>%
    drop_na(date) %>%
    dplyr::arrange(dplyr::across(all_of(station_col)), date) %>%
    group_by(across(all_of(station_col)), date) %>%
    mutate(grp = cumsum(c(1, diff(date) > 1)))
  
  unique_species <- na.omit(unique(df$species))
  result_list <- list()
  
  for (sp in unique_species) {
    sp_df <- df %>% arrange(desc(species == sp))
    
    sp_df2 <- sp_df %>%
      ungroup() %>%
      group_by(across(all_of(station_col)), grp, .add = TRUE) %>%
      dplyr::mutate(row_count = n()) %>%
      ungroup()
    
    sp_df3 <- sp_df2 %>%
      filter(row_count > 1) %>%
      group_by(across(all_of(station_col)), grp, .add = TRUE) %>%
      dplyr::reframe(
        across(where(is.numeric), ~ ifelse(all(is.na(.)), NA, mean(., na.rm = TRUE))),
        date = coalesce(first(date[!is.na(probability)]), first(date)),
        probability = ifelse(all(is.na(probability)), NA_character_, first(na.omit(as.character(probability)))),
        probability_AO = ifelse(all(is.na(probability_AO)), NA_character_, first(na.omit(as.character(probability_AO)))),
        probability_AT = ifelse(all(is.na(probability_AT)), NA_character_, first(na.omit(as.character(probability_AT)))),
        probability_AM = ifelse(all(is.na(probability_AM)), NA_character_, first(na.omit(as.character(probability_AM)))),
        probability_Aspp = ifelse(all(is.na(probability_Aspp)), NA_character_, first(na.omit(as.character(probability_Aspp)))),
        strat = ifelse(all(is.na(strat)), NA_character_, first(na.omit(as.character(strat)))),
        species = first(species)
      ) %>%
      dplyr::select(-grp) %>%
      ungroup()
    
    sp_result <- full_join(sp_df3, sp_df2 %>% filter(row_count == 1)) %>%
      dplyr::select(-grp, -row_count)
    
    result_list[[sp]] <- sp_result
  }
  
  result <- dplyr::bind_rows(result_list)
  
  return(result)
}




combine_close_stations <- function(df, threshold = 0.025) {
  library(dplyr)
  
  # Create a unique dataframe with station, lat, and lon columns
  unique_df <- df %>%
    dplyr::ungroup() %>%
    dplyr::select(station, lat, lon) %>%
    unique() %>%
    drop_na(lat)
  
  # Create a list to store station groups
  station_groups <- list()

  # Loop through each station
  for (i in 1:nrow(unique_df)) {
    # Check if the station has already been assigned to a group
    if (!unique_df$station[i] %in% unlist(station_groups)) {
      # Find stations within the threshold
      close_stations <- unique_df$station[abs(unique_df$lat - unique_df$lat[i]) <= threshold &
                                            abs(unique_df$lon - unique_df$lon[i]) <= threshold]
      
      # Add the group of close stations to the list
      station_groups[[length(station_groups) + 1]] <- close_stations
    }
  }
  
  # Combine station names within each group
  combined_stations <- sapply(station_groups, function(x) paste(x, collapse = ", "))
  
  # Create a new dataframe with combined station names, mean lat/lon
  new_df <- data.frame(
    station = unlist(combined_stations),
    lat = sapply(station_groups, function(x) mean(unique_df$lat[unique_df$station %in% x], na.rm = TRUE)),
    lon = sapply(station_groups, function(x) mean(unique_df$lon[unique_df$station %in% x], na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  new_df$old_station <- sapply(strsplit(new_df$station, ", "), function(x) x[1])

  # Update unique_df with combined station names and mean lat/lon
  unique_df$old_station <- NA
  for (i in 1:nrow(unique_df)) {
    station_name <- unique_df$station[i]
    match_indices <- sapply(new_df$station, function(x) str_detect(x, fixed(station_name)))
    if (!is.na(match_indices) && any(match_indices)) {
      match_row <- new_df[match_indices, ]
      unique_df$station[i] <- match_row$station[1]
      unique_df$lat[i] <- match_row$lat[1]
      unique_df$lon[i] <- match_row$lon[1]
      unique_df$old_station[i] <- station_name
    } else {
      unique_df$station[i] <- station_name
      unique_df$lat[i] <- unique_df$lat[i]
      unique_df$lon[i] <- unique_df$lon[i]
      unique_df$old_station[i] <- station_name
    }
  }
  
  return(unique_df)
}


# Function to find the closest stations of wind stations and phytoplankton stations (used for Danish data set)
find_closest_stations <- function(df1, df2) {
  closest_stations <-
    data.frame(station1 = character(),
               station2 = character(),
               distance = numeric())
  
  for (j in 1:nrow(df2)) {
    station2 <- df2[j, ]
    closest_station1 <- NULL
    min_distance <- Inf
    
    for (i in 1:nrow(df1)) {
      station1 <- df1[i, ]
      
      # Calculate the distance between station1 and station2
      distance <-
        distHaversine(station1[c("Lon", "Lat")], station2[c("lon", "lat")])
      
      # Update closest_station1 and min_distance if a smaller distance is found
      if (distance < min_distance) {
        closest_station1 <- station1$station
        min_distance <- distance
      }
    }
    
    # Store the closest station1 for the current station2
    closest_stations <-
      rbind(
        closest_stations,
        data.frame(
          station1 = closest_station1,
          station2 = station2$station,
          distance = min_distance
        )
      )
  }
  
  return(closest_stations)
}

# Define a function to combine 'strat' values
combine_strat <- function(strat_values) {
  if ("stratified" %in% strat_values & "not stratified" %in% strat_values) {
    return("stratified/not stratified")
  } else if ("stratified" %in% strat_values) {
    return("stratified")
  } else if ("not stratified" %in% strat_values) {
    return("not stratified")
  } else {
    return(NA)
  }
}

# Define a function to combine 'limiting_conditions' values
combine_limiting_conditions <- function(limiting_conditions_values) {
  if ("yes" %in% limiting_conditions_values & "no" %in% limiting_conditions_values) {
    return("yes/no")
  } else if ("yes" %in% limiting_conditions_values) {
    return("yes")
  } else if ("no" %in% limiting_conditions_values) {
    return("no")
  } else {
    return(NA)
  }
}

# Function to change DMM to DM format
dmm_to_dd <- function(dmm_coordinates) {
  dmm_parts <- strsplit(dmm_coordinates, " ") # Split degrees and minutes
  degrees <- as.numeric(dmm_parts[[1]][1])    # Convert degrees to numeric
  minutes <- as.numeric(dmm_parts[[1]][2])    # Convert minutes to numeric
  dd_coordinates <- degrees + (minutes / 60)  # Calculate decimal degrees
  return(dd_coordinates)
}

# add lag columns for parameters of interest. Function checks whether there is lagged data available in the respective lagged timeframes
# and otherwise prints NA
create_lagged_columns <- function(data, column_name, num_weeks) {
  lagged_values <- vector("double", length = nrow(data))
  
  for (i in 2:nrow(data)) {
    if (!is.na(data$date[i]) && !is.na(data$date[i - 1]) &&
        data$date[i] - data$date[i - 1] <= num_weeks * 7) {
      lagged_values[i] <- data[[column_name]][i - 1]
    } else {
      lagged_values[i] <- NA
    }
  }
  
  col_name <- paste0(column_name, "_lag_", num_weeks)
  data[[col_name]] <- lagged_values
  
  return(data)
}

# assign regions to each station depending on their geographical position
assign_regions_to_stations <- function(data, region_boxes) {
  regions <- character(length = nrow(data))
  
  for (i in 1:nrow(data)) {
    matching_regions <- character()
    
    for (j in 1:nrow(region_boxes)) {
      if (!is.na(data$lat[i]) &&
          !is.na(region_boxes$min_lat[j]) &&
          !is.na(region_boxes$max_lat[j]) &&
          !is.na(data$lon[i]) &&
          !is.na(region_boxes$min_lon[j]) &&
          !is.na(region_boxes$max_lon[j]) &&
          between(data$lat[i],
                  region_boxes$min_lat[j],
                  region_boxes$max_lat[j]) &&
          between(data$lon[i],
                  region_boxes$min_lon[j],
                  region_boxes$max_lon[j])) {
        matching_regions <- c(matching_regions, region_boxes$region[j])
      }
    }
    
    if (length(matching_regions) > 0) {
      regions[i] <-
        matching_regions[1]  # Assign the first matching region
    } else {
      regions[i] <- "unassigned"
    }
  }
  
  # Add the 'regions' column to 'data'
  data$regions <- regions
  
  return(data)
}

# Define a helper function for log transformation
log_transform <- function(x) {
  result <- ifelse(x != 0 & !is.infinite(x), log(x), NA)
  return(result)
}

# check_and_fit <- function(data,
#                           each_group,
#                           station_or_region,
#                           result_month_name,
#                           result_year_name,
#                           predicted_probs_name,
#                           results_list_name,
#                           result_doy_name,
#                           result_doy_gam_name,
#                           probability_column) {
#   if (!exists(results_list_name, envir = .GlobalEnv)) {
#     assign(
#       results_list_name,
#       list(
#         result_month = data.frame(),
#         result_year = data.frame(),
#         predicted_probs = data.frame()
#       ),
#       envir = .GlobalEnv
#     )
#   }
#   
#   # Initialize local data frames
#   result_month <- data.frame()
#   result_year <- data.frame()
#   result_doy <- data.frame()
#   result_doy_gam <- data.frame()
#   predicted_probs_df <- data.frame()
#   
#   results_list <- get(results_list_name, envir = .GlobalEnv)
#   # Determine the column name to filter and group by based on the grouping_var
#   if (station_or_region == "station") {
#     grouping_var <- "station"
#   } else if (station_or_region == "region") {
#     grouping_var <- "regions"
#   }
#   data <- data %>%
#     filter(!is.na(!!rlang::sym(probability_column)) &
#              station %in% stations_matching_condition)
#   
#   present_years_only <- data %>%
#     filter(!!rlang::sym(grouping_var) == each_group)
#   present_years_only <-
#     present_years_only %>% filter(year >= (
#       present_years_only %>% filter(!!rlang::sym(probability_column) == 1) %>% pull(year) %>% min()
#     ))
#   
#   if (sum(present_years_only[[probability_column]] == 0) >= 2 &
#       sum(present_years_only[[probability_column]] == 1) >= 2 &
#       length(unique(present_years_only$year)) > 3 &
#       length(unique(present_years_only$month)) > 3) {
#     mod <- as.formula(sprintf("factor(%s) ~  factor(month)  + 0", probability_column))
#     M1 <-
#       glm(formula = mod,
#           data = present_years_only,
#           family = "binomial") # monthly probability pattern
#     
#     mod_doy <- as.formula(sprintf("factor(%s) ~ factor(doy) + 0", probability_column))
#     M1a <-
#       glm(formula = mod_doy,
#           data = present_years_only,
#           family = "binomial") # doyly probability pattern
#     
#     if (length(levels(as.factor(
#       as.character(present_years_only$limiting_conditions)
#     ))) >= 2) {
#       M1_lim <-
#         glm(
#           factor(limiting_conditions) ~ factor(month) + 0,
#           data = present_years_only,
#           family = "binomial"
#         ) # monthly probability pattern
#     }
#     
#     if (length(levels(as.factor(
#       as.character(present_years_only$strat)
#     ))) >= 2) {
#       M1_strat <-
#         glm(factor(strat) ~ factor(month) + 0,
#             data = present_years_only,
#             family = "binomial") # monthly probability pattern
#     }
#     mod_gam <- as.formula(sprintf("factor(%s) ~ s(doy, bs = 'cp')" , probability_column))
#     
#     M1a_gam <- gam(
#       data = present_years_only,
#       formula = mod_gam,
#       family = binomial,
#       knots = list(x = c(0, 365))
#     )
#    
#     M1a_gam_predict <- predict(M1a_gam, se.fit = TRUE)
#     
#     mod_year <- as.formula(sprintf("factor(%s) ~ factor(year) + 0", probability_column))
#     M2 <-
#       glm(formula = mod_year,
#           data = present_years_only,
#           family = "binomial") # yearly probability pattern
#     
#     # Convert month column to character
#     present_years_only$month <-
#       as.character(present_years_only$month)
#     
#     CI_M1 = confint(M1, level = 0.95)[, 1]
#     # Bind coefficients and meta-data in result_month/result_year
#     result_month <- data.frame(
#       data = exp(coef(M1)) / (1 + exp(coef(M1))),
#       month = unlist(str_extract_all(
#         names(coef(M1)), "\\d+\\.?\\d*|\\.\\d+"
#       )),
#       CI = exp(CI_M1) / (1 + exp(CI_M1))
#     )
#     
#     if (exists("M1_lim")) {
#       limiting_conditions_df <- data.frame(
#         month = unlist(str_extract_all(
#           names(coef(M1_lim)), "\\d+\\.?\\d*|\\.\\d+"
#         )),
#         limiting_conditions = exp(coef(M1_lim)) / (1 + exp(coef(M1_lim)))
#       )
#     }
#     if (exists("M1_strat")) {
#       strat_df <- data.frame(month = unlist(str_extract_all(
#         names(coef(M1_strat)), "\\d+\\.?\\d*|\\.\\d+"
#       )),
#       strat = exp(coef(M1_strat)) / (1 + exp(coef(M1_strat))))
#     }
#     # Merge the data frames for result_month
#     if (exists("limiting_conditions_df")) {
#       result_month <-
#         merge(result_month,
#               limiting_conditions_df,
#               by = "month",
#               all = TRUE)
#     }
#     if (exists("strat_df")) {
#       result_month <-
#         merge(result_month, strat_df, by = "month", all = TRUE)
#     }
#     
#     result_month[[grouping_var]] = as.factor(each_group)
#     
#     CI_M2 = confint(M2, level = 0.95)[, 1]
#     
#     # Create data frames for result_year
#     result_year <- data.frame(
#       data = exp(coef(M2)) / (1 + exp(coef(M2))),
#       year = as.numeric(unlist(
#         str_extract_all(names(coef(M2)), "\\d+\\.?\\d*|\\.\\d+")
#       )),
#       CI = exp(CI_M2) / (1 + exp(CI_M2))
#     )
#     
#     result_year[[grouping_var]] = as.factor(each_group)
#     
#     # Bind coefficients and meta-data in result_month/result_year
#     result_doy <- bind_rows(result_doy,
#                             data.frame(
#                               data = exp(coef(M1a)) / (1 + exp(coef(M1a))),
#                               doy = unlist(str_extract_all(
#                                 names(coef(M1a)), "\\d+\\.?\\d*|\\.\\d+"
#                               ))
#                             ))
#     result_doy[[grouping_var]] = as.factor(each_group)
#     
#     # Bind coefficients and meta-data in result_month/result_year
#     result_doy_gam <- data.frame(data = NA,
#                                  doy = present_years_only$doy,
#                                  CI = NA,
#                                  p_val = NA)
#     if (exists("M1a_gam")) {
#       result_doy_gam$p_val = summary(M1a_gam)$p.pv[[1]]
#     }
#     
#     if (summary(M1a_gam)$p.pv < 0.1) {
#       result_doy_gam$data <-
#         exp(M1a_gam_predict$fit) / (1 + exp(M1a_gam_predict$fit))
#       result_doy_gam$CI <-
#         exp(M1a_gam_predict$se.fit) / (1 + exp(M1a_gam_predict$se.fit))
#     }
#     
#     result_doy_gam[[grouping_var]] = as.factor(each_group)
#     
#     # mod_pred <-
#     #   as.formula(
#     #     sprintf(
#     #       "cbind(as.numeric(as.character((%s))), 1- as.numeric(as.character((%s)))) ~ year",
#     #       probability_column,
#     #       probability_column
#     #     )
#     #   )
#     
#     model <-
#       glm(formula = cbind(logistic, 1 - logistic) ~ year,
#           data = present_years_only,
#           family = binomial)
#     
#     # Create a sequence of years for prediction
#     min_year <- min(present_years_only$year)
#     max_year <- max(present_years_only$year)
#     years <- seq(min_year, max_year, by = 1)
#     
#     # Predict probabilities using the logistic regression model
#     predicted_probs <-
#       predict(model, newdata = data.frame(year = years), type = "response")
#     
#     # Create a dataframe with predicted probabilities and station information
#     predicted_probs_df <-  bind_rows(
#       predicted_probs_df,
#       data.frame(year = years,
#                  predicted_probability = predicted_probs)
#     )
#     predicted_probs_df[[grouping_var]] = as.factor(each_group)
#     
#     # Bind the results to the specified output data frames within results_list
#     results_list[[result_month_name]] <-
#       bind_rows(results_list[[result_month_name]], result_month)
#     results_list[[result_year_name]] <-
#       bind_rows(results_list[[result_year_name]], result_year)
#     results_list[[predicted_probs_name]] <-
#       bind_rows(results_list[[predicted_probs_name]], predicted_probs_df)
#     results_list[[result_doy_name]] <-
#       bind_rows(results_list[[result_doy_name]], result_doy)
#     results_list[[result_doy_gam_name]] <-
#       bind_rows(results_list[[result_doy_gam_name]], result_doy_gam)
#     
#     assign(results_list_name, results_list, envir = .GlobalEnv)
#     
#     return(results_list_name)
#   }
# }

check_and_fit <- function(data,
                          each_group,
                          station_or_region,
                          result_month_name,
                          result_year_name,
                          predicted_probs_name,
                          results_list_name,
                          result_doy_name,
                          result_doy_gam_name) {
  if (!exists(results_list_name, envir = .GlobalEnv)) {
    assign(
      results_list_name,
      list(
        result_month = data.frame(),
        result_year = data.frame(),
        predicted_probs = data.frame()
      ),
      envir = .GlobalEnv
    )
  }
  
  # Initialize local data frames
  result_month <- data.frame()
  result_year <- data.frame()
  result_doy <- data.frame()
  result_doy_gam <- data.frame()
  predicted_probs_df <- data.frame()
  
  results_list <- get(results_list_name, envir = .GlobalEnv)
  # Determine the column name to filter and group by based on the grouping_var
  if (station_or_region == "station") {
    grouping_var <- "station"
  } else if (station_or_region == "region") {
    grouping_var <- "regions"
  }
  data <- data %>%
    filter(!is.na(probability) &
             station %in% stations_matching_condition)
  
  present_years_only <- data %>%
    filter(!!rlang::sym(grouping_var) == each_group) %>% drop_na(cells_L)
  present_years_only2 <-
    present_years_only %>% filter(year >= (
      present_years_only %>% filter(probability == 1) %>% pull(year) %>% min()
    ))
  
  # # Define a function to perform the operations
  # process_data <- function(data) {
  #   data %>%
  #     dplyr::select(station, date, year, month, doy, day, probability, limiting_conditions, strat, temp) %>%
  #     distinct() %>%
  #     group_by(date) %>%
  #     filter(!(any(probability == 1) & probability == 0)) %>%
  #     slice(ifelse(any(probability == 1), which(probability == 1), 1)) %>%
  #     ungroup()
  # }
  # 
  # # Apply the function to both datasets
  # present_years_only <- process_data(present_years_only)
  # 
  print(each_station) 
  if (sum(present_years_only$probability == 0) >= 2 &
      # at least two absent entries
      sum(present_years_only$probability == 1) >= 2 &
      # at least three absent entries
      length(unique(present_years_only$year)) > 3 &
      # at least three distinct years
      length(unique(present_years_only$month)) > 3) {
    # at least two distinct months)
    # if (length(levels(as.factor(
    #   as.character(present_years_only$probability)
    # ))) >= 2) {
    M1 <-
      glm(factor(probability) ~ factor(month) + 0,
          data = present_years_only %>% drop_na(probability),
          family = "binomial") # monthly probability pattern
    
    M1a <-
      glm(factor(probability) ~ factor(doy) + 0,
          data = present_years_only %>% drop_na(probability),
          family = "binomial") # doyly probability pattern

    if (length(levels(as.factor(
      as.character(present_years_only$limiting_conditions)
    ))) >= 2) {
      M1_lim <-
        glm(
          factor(limiting_conditions) ~ factor(month) + 0,
          data = present_years_only %>% drop_na(limiting_conditions),
          family = "binomial"
        ) # monthly probability pattern
    }
    
    if (length(levels(as.factor(
      as.character(present_years_only$strat)
    ))) >= 2) {
      M1_strat <-
        glm(factor(strat) ~ factor(month) + 0,
            data = present_years_only %>% drop_na(strat),
            family = "binomial") # monthly probability pattern
    }
    
    M1a_gam <- gam(
      data = present_years_only %>% drop_na(probability),
      probability ~ s(doy, bs = "cp"),
      family = binomial,
      knots = list(x = c(0, 365))
    )
    M1a_gam_predict <- predict(M1a_gam, se.fit = TRUE)
    
    M2 <-
      glm(factor(probability) ~ factor(year) + 0,
          data = present_years_only  %>% drop_na(probability),
          family = "binomial") # yearly probability pattern
    
    # Convert month column to character
    present_years_only$month <-
      as.character(present_years_only$month)
    
    CI_M1 = confint(M1, level = 0.95)[, 1]
    # Bind coefficients and meta-data in result_month/result_year
    result_month <- data.frame(
      data = exp(coef(M1)) / (1 + exp(coef(M1))),
      month = unlist(str_extract_all(
        names(coef(M1)), "\\d+\\.?\\d*|\\.\\d+"
      )),
      CI = exp(CI_M1) / (1 + exp(CI_M1))
    )
    
    if (exists("M1_lim")) {
      limiting_conditions_df <- data.frame(
        month = unlist(str_extract_all(
          names(coef(M1_lim)), "\\d+\\.?\\d*|\\.\\d+"
        )),
        limiting_conditions = exp(coef(M1_lim)) / (1 + exp(coef(M1_lim)))
      )
    }
    if (exists("M1_strat")) {
      strat_df <- data.frame(month = unlist(str_extract_all(
        names(coef(M1_strat)), "\\d+\\.?\\d*|\\.\\d+"
      )),
      strat = exp(coef(M1_strat)) / (1 + exp(coef(M1_strat))))
    }
    # Merge the data frames for result_month
    if (exists("limiting_conditions_df")) {
      result_month <-
        merge(result_month,
              limiting_conditions_df,
              by = "month",
              all = TRUE)
    }
    if (exists("strat_df")) {
      result_month <-
        merge(result_month, strat_df, by = "month", all = TRUE)
    }
    
    result_month[[grouping_var]] = as.factor(each_group)
    
    CI_M2 = confint(M2, level = 0.95)[, 1]
    
    # Create data frames for result_year
    result_year <- data.frame(
      data = exp(coef(M2)) / (1 + exp(coef(M2))),
      year = as.numeric(unlist(
        str_extract_all(names(coef(M2)), "\\d+\\.?\\d*|\\.\\d+")
      )),
      CI = exp(CI_M2) / (1 + exp(CI_M2))
    )
    
    result_year[[grouping_var]] = as.factor(each_group)
    
    # Bind coefficients and meta-data in result_month/result_year
    result_doy <- bind_rows(result_doy,
                            data.frame(
                              data = exp(coef(M1a)) / (1 + exp(coef(M1a))),
                              doy = unlist(str_extract_all(
                                names(coef(M1a)), "\\d+\\.?\\d*|\\.\\d+"
                              ))
                            ))
    result_doy[[grouping_var]] = as.factor(each_group)
    
    # Bind coefficients and meta-data in result_month/result_year
    result_doy_gam <- data.frame(data = NA,
                                 doy = present_years_only$doy,
                                 CI = NA,
                                 p_val = NA)
    if (exists("M1a_gam")) {
      result_doy_gam$p_val = summary(M1a_gam)$p.pv[[1]]
    }
    
    if (summary(M1a_gam)$p.pv < 0.1) {
      result_doy_gam$data <-
        exp(M1a_gam_predict$fit) / (1 + exp(M1a_gam_predict$fit))
      result_doy_gam$CI <-
        exp(M1a_gam_predict$se.fit) / (1 + exp(M1a_gam_predict$se.fit))
    }
    
    result_doy_gam[[grouping_var]] = as.factor(each_group)
    
    model <-
      glm(cbind(as.numeric(as.character(probability)), 1 - as.numeric(as.character(probability))) ~ year,
          data = present_years_only2 %>% filter(year >= min((present_years_only2 %>% filter(year >= 2000))$year)) %>% drop_na(probability),
          family = binomial)
    
    # Create a sequence of years for prediction
    min_year <- min((present_years_only2 %>% filter(probability == 1 & year >= 2000))$year)
    max_year <- max(present_years_only2$year)
    max_year <- ifelse(max_year >= 2022, 2022, max_year)
    years <- seq(min_year, max_year, by = 1)
    
    # Predict probabilities using the logistic regression model
    predicted_probs <-
      predict(model, newdata = data.frame(year = years), type = "response")
    
    # Create a dataframe with predicted probabilities and station information
    predicted_probs_df <-  bind_rows(
      predicted_probs_df,
      data.frame(year = years,
                 predicted_probability = predicted_probs)
    )
    predicted_probs_df[[grouping_var]] = as.factor(each_group)
    
    # Bind the results to the specified output data frames within results_list
    results_list[[result_month_name]] <-
      bind_rows(results_list[[result_month_name]], result_month)
    results_list[[result_year_name]] <-
      bind_rows(results_list[[result_year_name]], result_year)
    results_list[[predicted_probs_name]] <-
      bind_rows(results_list[[predicted_probs_name]], predicted_probs_df)
    results_list[[result_doy_name]] <-
      bind_rows(results_list[[result_doy_name]], result_doy)
    results_list[[result_doy_gam_name]] <-
      bind_rows(results_list[[result_doy_gam_name]], result_doy_gam)
    
    assign(results_list_name, results_list, envir = .GlobalEnv)
    
    return(results_list_name)
  }
}

check_and_fit2 <- function(data,
                           grouping_var,
                           each_group,
                           result_month_name,
                           each_year,
                           results_list_name,
                           year_column,
                           probability_columns) {
  # Check if the results list exists, and if not, initialize it
  if (!exists(results_list_name, envir = .GlobalEnv)) {
    assign(results_list_name, list(), envir = .GlobalEnv)
  }
  
  results_list <- get(results_list_name, envir = .GlobalEnv)
  
  for (prob_col in probability_columns) {
    # Initialize local data frame
    result_month <- data.frame()
    
    data_subset <- data %>%
      filter(!is.na(!!rlang::sym(prob_col)) & !!rlang::sym(year_column) == each_year)
    
    if (grouping_var == "each_region") {
      data_subset <- data_subset %>% filter(region == each_group)
    }
    
    # Exclude stations without both "absent" and "present" observations in the current year
    stations_with_both <- data_subset %>%
      group_by(station) %>%
      filter(all(c(0, 1) %in% !!rlang::sym(prob_col))) %>%
      pull(station) %>% unique()
    
    year_to_fit <- data_subset %>%
      filter(station %in% stations_with_both)
    
    if (sum(year_to_fit[[prob_col]] == 0) >= 2 &
        sum(year_to_fit[[prob_col]] == 1) >= 2 &
        length(unique(year_to_fit$month)) > 3) {
      mod <- as.formula(sprintf("factor(%s) ~ factor(month) + 0", prob_col))
      M3 <-
        glm(formula = mod,
            data = year_to_fit,
            family = "binomial") # monthly probability pattern
      
      CI = confint(M3, level = 0.95)[, 1]
      
      result_month <- bind_rows(
        result_month,
        data.frame(
          data = exp(coef(M3)) / (1 + exp(coef(M3))),
          month = unlist(str_extract_all(
            names(coef(M3)), "\\d+\\.?\\d*|\\.\\d+"
          )),
          year = each_year,
          region = ifelse(
            as.character(grouping_var) %in% c("each_region", "region"),
            each_group,
            NA
          ),
          CI = exp(CI) / (1 + exp(CI))
        )
      )
      
      # Add the result to the results list
      results_list[[paste0(result_month_name, "_", prob_col)]] <-
        bind_rows(results_list[[paste0(result_month_name, "_", prob_col)]], result_month)
    }
  }
  
  assign(results_list_name, results_list, envir = .GlobalEnv)
  return(results_list_name)
}

calculate_seasonal_mean <- function(data,
                                    probability_columns) {
  
  seasonal_means <- list()
  result_coef_list <- list() # Initialize an empty list
  
  # Prepare parameters
  parameters <-
    c("NH4",
      "NO3",
      "sal",
      "temp",
      "TN",
      "PO4",
      "chl",
      # "wind_ms",
      # "silicate",
      # "C",
      "DIP",
      "DIN",
      "NO3_PO4",
      # "Si_N",
      # "DON",
      # "DON_DIN",
      "C_chl"
      )
  
  for (each_station in unique(data$station)) {
    # Filter data for the current station
    station_data <- data %>%
      filter(station == each_station)
    
    for (prob_col in probability_columns) {
      min_year_of_presence <- station_data %>% filter(!!rlang::sym(prob_col) == 1) %>% pull(year) %>% min()
      ifelse(min_year_of_presence < 2000, 2000, min_year_of_presence)
      present_years_only <-
        station_data %>%
        filter(year >= min_year_of_presence & month >= 5 & month <= 10)
      
      # Perform the first glm for seasonal probability
      if (length(levels(as.factor(as.character(present_years_only[[prob_col]])))) >= 2) {
        mod <- as.formula(sprintf("factor(%s) ~ factor(month) + 0", prob_col))
        seasonal_mean_glm <- glm(formula = mod,
                                 data = present_years_only,
                                 family = "binomial")
        coef_exp <- exp(coef(seasonal_mean_glm))
        result_coef_list[[prob_col]] <- data.frame(data = coef_exp / (1 + coef_exp)) %>%
          rownames_to_column("time") %>%
          filter(time %like% "month") %>%
          mutate(time = parse_number(time), station = each_station, parameter = prob_col)
      }
    
    if (length(levels(as.factor(as.character(present_years_only$limiting_conditions)))) >= 2) {
      M1_lim <- glm(factor(limiting_conditions) ~ factor(month) + 0,
                    data = present_years_only,
                    family = "binomial")
      coef_exp2 <- exp(coef(M1_lim))
      result_coef_list$limiting_conditions <- data.frame(data = coef_exp2 / (1 + coef_exp2)) %>%
        rownames_to_column("time") %>%
        filter(time %like% "month") %>%
        mutate(time = parse_number(time), station = each_station, parameter = "limiting_conditions")
    }
    
    if (length(levels(as.factor(as.character(present_years_only$strat)))) >= 2) {
      M1_strat <- glm(factor(strat) ~ factor(month) + 0,
                      data = present_years_only,
                      family = "binomial")
      coef_exp3 <- exp(coef(M1_strat))
      result_coef_list$strat <- data.frame(data = coef_exp3 / (1 + coef_exp3)) %>%
        rownames_to_column("time") %>%
        filter(time %like% "month") %>%
        mutate(time = parse_number(time), station = each_station, parameter = "strat")
    }
    }
    
    # Join the results from the list
    result_coef_all <- Reduce(full_join, result_coef_list)
    
    results_coef <- lapply(parameters, function(x) NULL)  
    
    for (each_parameter in parameters) {
      print(each_station)
      print(each_parameter)
      
      present_years_only <-
        station_data %>%
        filter(year >= 2000 & month >= 5 & month <= 10)
      
        param_data <- present_years_only %>% drop_na(!!sym(each_parameter)) %>%
          filter(!is.infinite(!!sym(each_parameter)))
        # Identify and handle extreme outliers
        param_data <- param_data %>%
          filter(!!sym(each_parameter) < quantile(!!sym(each_parameter), 0.99, na.rm = TRUE))

        if (nlevels(factor(param_data[[prob_col]])) == 2 & length(unique(param_data[[each_parameter]])) > 2) {
          # mod <- as.formula(sprintf("%s ~ factor(month) + factor(%s) + 0", each_parameter, prob_col))
          mod <- as.formula(sprintf("%s ~ factor(month)  + 0", each_parameter))
          
          if (!(each_parameter %in% c("temp", "wind_ms", "sal", "DON_DIN", "NO3_PO4", "Si_N", "DON", "chl", "C_chl"))) {
            # Check if there are negative values in the parameter column
            if (any(param_data[[each_parameter]] < 0)) {
              # Skip fitting the model if there are negative values
              next
            }
            # add small constant to handle 0s
            param_data[[each_parameter]] <-
              param_data[[each_parameter]] + min(param_data[[each_parameter]][param_data[[each_parameter]] > 0], na.rm = TRUE) / 10
            seasonal_mean_glm2 <- glm(data = param_data,
                                      formula = mod,
                                      family = gaussian(link = "log"))
          } else {
            seasonal_mean_glm2 <- glm(data = param_data,
                                      formula = mod,
                                      family = gaussian)
          }
          
          # Extract model coefficients and residuals
          model_coefficients <- coef(seasonal_mean_glm2)
          # model_residuals <- residuals(seasonal_mean_glm2)
          
          # Calculate the smearing estimate
          # smearing_estimate <- exp(mean(model_residuals))
          
          # Apply the smearing estimate to the model coefficients, except for "temp" and "wind_ms"
          if (!(each_parameter %in% c("temp", "wind_ms", "sal", "DON_DIN", "NO3_PO4", "Si_N", "DON", "chl", "C_chl"))) {
            result_coef_month2 <- data.frame(data = exp(model_coefficients)) %>%
              rownames_to_column("time") %>%
              filter(time %like% "month") %>%
              mutate(time = parse_number(time), station = each_station, parameter = each_parameter)
          } else {
            result_coef_month2 <- data.frame(data = model_coefficients) %>%
              rownames_to_column("time") %>%
              filter(time %like% "month") %>%
              mutate(time = parse_number(time), station = each_station, parameter = each_parameter)
          }
          
          # Add the result to the seasonal_means dataframe
          results_coef[[each_parameter]] <- bind_rows(results_coef[[each_parameter]], result_coef_month2)
        }
      }
    
    # Remove NULL elements from results_coef
    results_coef <- results_coef[!sapply(results_coef, is.null)]
    
    results_coef <- Reduce(full_join, results_coef)
    if (length(results_coef) > 0)  {
      results_coef_all2 <- full_join(results_coef, result_coef_all)
    } else {
      results_coef_all2 <- result_coef_all
    }
    
    seasonal_means[[each_station]] <- results_coef_all2
  }
  
  # Return the updated results dataframe
  return(seasonal_means)
}


# function to find all years with minimum one present observation of A. pseudogonyaulax
# restricting the dataset to these years reduces standard deviations
find_years_of_presence <- function(data, probability_columns) {
  years_of_presence <- list()
  
  for (prob_col in probability_columns) {
    years_of_presence_df <- data.frame(station = " ", year = " ")
    
    for (each_station in unique(data$station)) {
      station_years <- data %>%
        filter(station == each_station) %>%
        dplyr::select(year, station, !!rlang::sym(prob_col)) %>%
        filter(!!rlang::sym(prob_col) == "present") %>%
        as.data.frame(station = station, year = year) %>%
        dplyr::select(station, year) %>%
        distinct()
      years_of_presence_df <- rbind(years_of_presence_df, station_years)
    }
    
    years_of_presence_df <- years_of_presence_df[-1, ]
    years_of_presence_df$year <- as.numeric(years_of_presence_df$year)
    years_of_presence[[prob_col]] <- years_of_presence_df
  }
  
  return(years_of_presence)
}
# function to find all years with minimum one present observation of A. pseudogonyaulax
# restricting the dataset to these years reduces standard deviations
find_years_since_first_observation <- function(data, probability_columns) {
  years_since_first_observation <- list()
  
  for (prob_col in probability_columns) {
    years_df <- data %>%
      filter(!!rlang::sym(prob_col) == 1) %>%
      group_by(station) %>%
      mutate(min_year_present = min(year)) %>%
      filter(year >= min_year_present) %>%
      dplyr::select(station, year, min_year_present) %>%
      distinct() %>%
      mutate(year_since_first_observation = year - min_year_present + 1,
             species = prob_col)
    
    years_since_first_observation[[prob_col]] <- years_df
  }
  
  years_since_first_observation_df <- do.call(rbind, years_since_first_observation)
  years_since_first_observation_df <- years_since_first_observation_df %>%
    dplyr::select(station, year, year_since_first_observation, min_year_present, species)
  
  return(years_since_first_observation_df)
}

  
# # Create sinusoidal and cosinusoidal functions
sin_function <- function(doy, amplitude, phase) {
   amplitude * sin(2 * pi * doy / 365 + phase)
}

cos_function <- function(doy, amplitude, phase) {
  amplitude * cos(2 * pi * doy / 365 + phase)
}

# Create sinusoidal and cosinusoidal functions
# sin_function <- function(doy) {
# sin(2 * pi * doy / 365 )
# }
# 
# cos_function <- function(doy) {
# cos(2 * pi * doy / 365 )
# }

# Define a function to fit the model
fit_model <- function(data, parameter) {
  # Use summary statistics to estimate starting parameters
  amplitude_init <- 2 
  phase_init <- 0  #

  nls_formula <- as.formula(
    sprintf(
      "%s ~ sin_function(doy, amplitude, phase) + cos_function(doy, amplitude, phase)",
      parameter
    )
  )
  
  nls_fit <- nls(
    nls_formula,
    data = data,
    start = c(amplitude = amplitude_init, phase = phase_init),
    control = nls.control(maxiter = 50, warnOnly = T, minFactor = 10^-7),
    lower = c(0.2, 0),
    upper = c(4, 5)
  )
  
  return(coef(nls_fit))
}
process_data <- function(data,
                         years_of_presence2,
                         parameters_of_interest,
                         station_or_region_col,
                         each_station,
                         each_region,
                         probability_columns) {
  parameter_results <- list()
  
  for (prob_col in probability_columns) {
    parameter_results[[prob_col]] <- lapply(parameters_of_interest, function(parameter) {
      print(parameter)
      param_data <- data %>%
        drop_na(parameter)
      
      years_of_presence2 <- years_of_presence2 %>% filter(species == prob_col) %>% pull(year)
      
      # Conditionally apply the year filter if years_of_presence2 is not NULL
      if (!is.null(years_of_presence2)) {
        param_data <- param_data %>%
          filter(year %in% years_of_presence2 & month <= 10 & month >= 5)
      }
      
      if (!(parameter %in% c(parameters_of_interest[grep("temp", parameters_of_interest)], "strat", "limiting_conditions"))) {
        param_data[[parameter]] <- log_transform(param_data[[parameter]])
       }
      
      if (nlevels(factor(param_data[[prob_col]])) >= 2 &&
          length(levels(factor(
            param_data %>%
            filter(!!rlang::sym(prob_col) == 0) %>%
            pull(parameter)
          ))) >= 5 &&
          length(levels(factor(
            param_data %>%
            filter(!!rlang::sym(prob_col) == 1) %>%
            pull(parameter)
          ))) >= 5) {
        
        result <- fit_model(param_data, parameter)
        
        mod <- as.formula(
          sprintf(
            "%s ~ sin_function(doy, %f, %f) + cos_function(doy, %f, %f) + factor(%s) + 0",
            parameter,
            result[["amplitude"]],
            result[["phase"]],
            result[["amplitude"]],
            result[["phase"]],
            prob_col
          )
        )
        
        M <- glm(data = param_data,
                 formula = mod,
                 family = gaussian)
        
        if (parameter %in% c("strat", "limiting_conditions")){
          M <- glm(data = param_data,
                   formula = parameter ~ year + month + factor(probability) + 0,
                   family = "binomial")
        }
        
        print(summary(M))
        
        if (length(summary(M)$coefficients[, 1]) >= 4 &&
            !is.na(summary(M)$coefficients[3, 4]) &&
            !is.na(summary(M)$coefficients[4, 4])) {
          coefficients_data <- data.frame(
            "parameter" = rep(parameter, each = 2),
            mean = as.numeric(summary(M)$coefficients[3:4, 1]),
            sd = as.numeric(summary(M)$coefficients[3:4, 2]),
            p_value = as.numeric(summary(M)$coefficients[3:4, 4]),
            p_sin = summary(M)$coefficients[1, 4],
            p_cos = summary(M)$coefficients[2, 4],
            treat = rep(c("absent", "present"), each = 1),
            amplitude = rep(result[["amplitude"]], each = 2),
            phase = rep(result[["phase"]], each = 2),
            species = rep(prob_col, each = 2)
          )
          if (station_or_region_col == "regions") {
            coefficients_data$regions = rep(each_region, each = 2)
          } else if (station_or_region_col == "station") {
            coefficients_data$station = rep(each_station, each = 2)
          }
          
          sample_size_data <- data.frame(
            "parameter" = rep(parameter, each = 2),
            present = sum(param_data[[prob_col]] == 1),
            absent = sum(param_data[[prob_col]] == 0),
            species = rep(prob_col, each = 2)
          )
          if (station_or_region_col == "regions") {
            sample_size_data$regions = rep(each_region, each = 2)
          } else if (station_or_region_col == "station") {
            sample_size_data$station = rep(each_station, each = 2)
          }
          
          return(
            list(
              coefficients_data = coefficients_data,
              sample_size_data = sample_size_data
            )
          )
        }
      } else {
        return(NULL)
      }
    })
  }
  
  # Return the captured output
  return(parameter_results)
}


backtransformation_of_log_transformed_coefficients <- function(Coef_stations_back, parameters_of_interest, grouping_var) {
  # Initialize an empty data frame to store backtransformed values
  Coef_stations_back3 <- data.frame()
  
  # Determine the column name to filter and group by based on the grouping_var
  if (grouping_var == "each_station") {
    filter_col <- "station"
  } else if (grouping_var == "each_region") {
    filter_col <- "regions"
  }
  
  for (grouping_value in unique(Coef_stations_back[[filter_col]])) {
    for (each_parameter in unique(Coef_stations_back$parameter)) {
      parameter_subset <- Coef_stations_back %>%
        filter(!!sym(filter_col) == grouping_value, parameter == each_parameter)
      
      unique_species <- unique(parameter_subset$species)
      unique_treats <- unique(parameter_subset$treat)
      
      for (each_species in unique_species) {
        for (each_treat in unique_treats) {
          Coef_stations_back2_sub <- parameter_subset %>%
            filter(species == each_species, treat == each_treat)
          
          if (nrow(Coef_stations_back2_sub) > 0 && !is.na(Coef_stations_back2_sub$sample_size) && Coef_stations_back2_sub$sample_size > 1) {
            # Extract mean and SD values
            mean_val <- Coef_stations_back2_sub$mean
            sd_val <- Coef_stations_back2_sub$sd
            
            if (each_parameter %in% c(unique(Coef_stations_back$parameter)[grep("temp", unique(Coef_stations_back$parameter))],
                                      "strat",
                                      "limiting_conditions")) {
              new_row <- data.frame(
                parameter = rep(each_parameter, each = 3),
                data = c(mean_val, sd_val, Coef_stations_back2_sub$p_value),
                type = c("mean", "sd", "p_value"),
                treat = rep(each_treat, each = 3),
                n = Coef_stations_back2_sub$sample_size[1],
                species = rep(each_species, each = 3)
              )
              new_row[[filter_col]] <- rep(grouping_value, each = 3)
              
              # Append the new row to the new data frame
              Coef_stations_back3 <- bind_rows(Coef_stations_back3, new_row)
            } else if (each_parameter != c("temp", "strat", "limiting_conditions") && !is.na(mean_val)) {
              # Apply bt.log() function to the mean and SD
              backtransformed_val <- bt.log(mean_val, sd_val, Coef_stations_back2_sub$sample_size)
              
              # Create a new row with the backtransformed value and replicate the original row information
              new_row <- data.frame(
                parameter = rep(each_parameter, each = 3),
                data = c(backtransformed_val[[1]], backtransformed_val[[4]], Coef_stations_back2_sub$p_value),
                type = c("mean", "sd", "p_value"),
                treat = rep(each_treat, each = 3),
                n = Coef_stations_back2_sub$sample_size[1],
                species = rep(each_species, each = 3)
              )
              new_row[[filter_col]] <- rep(grouping_value, each = 3)
              
              # Append the new row to the new data frame
              Coef_stations_back3 <- bind_rows(Coef_stations_back3, new_row)
            }
          }
        }
      }
    }
  }
  
  return(Coef_stations_back3)
}



# Define a function to calculate upr and lwr limit of the confidence interval of monthly and yearly probabilities of A. pseudogonyaulax presence
calculate_upr_lwr <- function(df) {
  if (nrow(df) > 0 && !is.null(df$CI)) {
    df$upr <- pmin(df$data + df$CI, 1)
    df$lwr <- pmax(df$data - df$CI, 0)
  } else {
    df$upr <- NA
    df$lwr <- NA
  }
  return(df)
}

# Function to create ggplot of monthly/yearly station/region probability data
create_plot <- function(data, model, group_var, x_var, title_var, save_path, filename) {
  
  # List to store plots
  plots <- list()
  
  gg <- ggplot(data, aes(x = !!sym(x_var), y = data, group = !!sym(group_var))) +
    ylab("probability") +
    xlab(x_var) +
    ggtitle(title_var) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 90),
      strip.background = element_blank(),
      axis.text = element_text()
    )

  if (x_var == "month") {
    gg <-
      gg + scale_x_discrete(breaks = 1:12,
                            labels = month.abb[1:12],
                            limits = factor(1:12)) + geom_line(aes(y = upr), col = "red") +
      geom_line(aes(y = lwr), col = "red") +
      geom_ribbon(aes(ymin = lwr, ymax = upr),
                  fill = "grey",
                  alpha = 0.25) +
      geom_line() +
      geom_point()
      
  }
  if (x_var == "year" && any(grepl("month", filename))) {
    gg <-
      gg + scale_x_discrete(breaks = 1:12,
                            labels = month.abb[1:12],
                            limits = factor(1:12)) + geom_line(aes(y = upr), col = "red") +
      geom_line(aes(y = lwr), col = "red") +
      geom_ribbon(aes(ymin = lwr, ymax = upr),
                  fill = "grey",
                  alpha = 0.25) +
      geom_line() +
      geom_point()
    
  }
  if (x_var == "year" &&
      any(grepl("year", save_path)) && !is.null(model)) {
    gg <-
      gg + geom_line(data = model,
                     aes(x = year, y = predicted_probability),
                     col = "blue") + geom_line(aes(y = upr), col = "red") +
      geom_line(aes(y = lwr), col = "red") +
      geom_ribbon(aes(ymin = lwr, ymax = upr),
                  fill = "grey",
                  alpha = 0.25) +
      geom_line() +
      geom_point()
    
  }
  if (x_var == "doy" & !(grepl("estuary|open|coastal", title_var)))
 {    
    gg <-
      gg + geom_line(data = model, aes(x = doy, y = data), col = "blue") +
      geom_point() +
      xlab(paste0("time (", x_var, ")"))
    # +
    #   geom_ribbon(data = model,
    #               aes(ymin = data - CI, ymax = data + CI),
    #               col = "grey")
  }
  if (grepl("estuary|open|coastal", title_var)) {
    gg <-
      gg + geom_line(data = model, aes(x = doy, y = data), col = "blue") 
    plots[[title_var]] <- gg
  }
  if (exists("manuscript_subset")) {
  if (identical(data, manuscript_subset)) {  # Check if the data is manuscript_subset    
    gg <- gg +
      scale_x_continuous(breaks = seq(2000, 2020, by = 2), limits = c(2000, 2020)) +
      scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
      theme(plot.title = element_blank()) +
      geom_text(aes(x = 2005, y = 0.6, label = title_var), size = 5) +
      geom_point()
    
  } }
  
  ggsave(filename = filename, plot = gg, path = save_path)
  
  # Return the list of plots
  return(plots)
}

calculate_seasonal_probability <- function(df, prob_col) {
  result_all <- data.frame()
  
  for (water_body_level in unique(df$water_body)) {
    present_years_only <- df %>%
      filter(water_body == water_body_level) %>%
      filter(station %in% unique_stations$station) %>%
      filter(year >= (df %>%
                        filter(!!rlang::sym(prob_col) == 1, water_body == water_body_level) %>%
                        pull(year) %>%
                        min())) %>%
      drop_na(prob_col)
    if (sum(present_years_only[[prob_col]] == 0) >= 2 &
        sum(present_years_only[[prob_col]] == 1) >= 2 &
        length(unique(present_years_only$year)) > 3 &
        length(unique(present_years_only$month)) > 3) {
      M1 <- glm(formula = as.formula(sprintf("%s ~ factor(month) + 0", prob_col)),
                data = present_years_only,
                family = "binomial")
      
      CI_M1 = confint(M1, level = 0.95)[, 1]
      
      result_month <- data.frame(
        data = exp(coef(M1)) / (1 + exp(coef(M1))),
        month = unlist(str_extract_all(names(coef(M1)), "\\d+\\.?\\d*|\\.\\d+")),
        CI = exp(CI_M1) / (1 + exp(CI_M1))
      )
      
      result_month <- calculate_upr_lwr(result_month) %>%
        mutate(station = "all stations", species = prob_col, water_body = water_body_level)
      
      result_all <- rbind(result_all, result_month)
    }}
  
  return(result_all)
}

create_plot2 <- function(data, model, group_var, x_var, title_var, each_region) {
  
  station_label <- ifelse(each_region == "Heiligendamm (Seebruecke)", "Heiligendamm", ifelse(each_region == "Heiligendamm_(Seebruecke)", "Heiligendamm", each_region))
  
  gg <-
    ggplot(data, aes(
      x = !!sym(x_var),
      y = data,
      group = !!sym(group_var)
    )) +
    geom_point()  +
    ylab("probability") +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = 'white'),
      legend.title = element_blank(),
      strip.background = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    ) +
    labs(x = NULL) +
    scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1)) +
    theme(plot.title = element_blank()) +
    geom_text(aes(x = 2006, y = 0.7, label = station_label), size = 3.5) +
    geom_line(data = model,
              aes(x = year, y = predicted_probability),
              col = "blue") + geom_line(aes(y = upr), col = "red") +
    geom_line(aes(y = lwr), col = "red") +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                fill = "grey",
                alpha = 0.25) +
    geom_line() +
    scale_x_continuous(breaks = seq(2000, 2020, by = 5), limits = c(2000, 2020)) +
    xlab(if (each_region %in% c("Heiligendamm_(Seebruecke)", "VEJ0006870")) x_var else NULL)
  
  if (!(each_region %in% c("Heiligendamm_(Seebruecke)", "VEJ0006870"))) {
    gg <- gg +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  if ((each_region %in% c("NOR409", "SLÄGGÖ", "KBH431", "Heiligendamm_(Seebruecke)"))) {
    gg <- gg +
      scale_y_continuous(
        limits = c(0, 0.7),
        breaks = seq(0, 0.7, by = 0.1),
        position = "right"
      )
    
  }
  
  return(gg)  # Return the ggplot object
}


create_plot3 <- function(data, model, group_var, x_var, station_number, each_region) {  
  gg <-
    ggplot(data, aes(
      x = !!sym(x_var),
      y = data,
      group = !!sym(group_var)
    )) +
    geom_point(size = 2)  +
    theme(
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_markdown(size = 12),
      axis.text.x = element_markdown(size = 10),
      axis.text.y.left = element_markdown(size = 10, margin = margin(r = 10)),
      panel.background = element_rect(fill = 'white'),
      legend.position = "none",
      plot.title.position = "plot",
      plot.margin = unit(c(0.2, 1, 0.2, 1), "cm")) +
    labs(x = NULL) +
    scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.2)) +
    geom_line(data = model,
              aes(x = year, y = predicted_probability),
              col = "blue", linewidth = 1) + 
    geom_line(aes(y = upr), col = "red") +
    geom_line(aes(y = lwr), col = "red", linewidth = 1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                fill = "grey",
                alpha = 0.25) +
    geom_line() +
    scale_x_continuous(breaks = seq(2000, 2020, by = 5), limits = c(2000, 2020)) +
    xlab(if (each_region %in% c("Heiligendamm", "VEJ0006870")) x_var else NULL)
  
  if (!(each_region %in% c("Heiligendamm", "VEJ0006870"))) {
    gg <- gg +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  if ((each_region %in% c("Heiligendamm", "VEJ0006870"))) {
    gg <- gg +
      theme(
        axis.text.x = element_markdown(size = 12))
  }
  
    if ((each_region %in% c("ANHOLT E", "Pricken", "KBH431", "Heiligendamm"))) {
    gg <- gg +
      scale_y_continuous(
        limits = c(0, 0.7),
        breaks = seq(0, 0.7, by = 0.2),
        position = "right"
      ) +
      annotate("text", x = 2005, y = 0.5, label = station_number, fontface = "bold", size = 6) +
      theme(axis.text.y.right = element_markdown(size = 10, margin = margin(l = 10)))
      
    
  }
  if (!(each_region %in% c("ANHOLT E", "Pricken", "KBH431", "Heiligendamm"))) {
    gg <- gg +
      annotate("text", x = 2017, y = 0.5, label = station_number, fontface = "bold", size = 6) 
  }
  
  return(gg)  # Return the ggplot object
}

