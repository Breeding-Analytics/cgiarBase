
summary_stats <- function(data, genotype, trial, traits, rep = NULL, block = NULL,
                          row = NULL, col = NULL, missing_data_threshold) {
  # Validate mandatory parameters
  if (is.null(data)) stop("Error: 'data' cannot be NULL.")
  if (is.null(genotype) || !genotype %in% names(data)) stop("Error: 'genotype' column not found.")
  if (is.null(traits) || !all(traits %in% names(data))) stop("Error: One or more 'traits' columns not found.")

  # # Automatically handle missing optional parameters
  # optional_params <- list(rep = rep, block = block, row = row, col = col)
  # for (param in names(optional_params)) {
  #   if (is.null(optional_params[[param]])) {
  #     message(paste0("No '", param, "' name column provided. Creating 'NA' column."))
  #     data[[param]] <- NA
  #   } else if (!optional_params[[param]] %in% names(data)) {
  #     stop(paste("Error: '", optional_params[[param]], "' column not found."))
  #   }
  # }

  for (i in traits) {
    if (!i %in% names(data)) {
      stop(paste("No '", i, "' column found"))
    }
    class_trait <- data[[i]] %>% class()
    if (!class_trait %in% c("numeric", "integer")) {
      stop(
        paste0("The class of the trait '", i, "' should be numeric or integer.")
      )
    }
  }

  # Identify traits with names that don't conform to R's naming conventions
  invalidTraitNames <- traits[make.names(traits) != traits]

  # Update 'traits' to include only those with valid R names
  validTraitNames <- traits[make.names(traits) == traits]

  # If there are invalid trait names, inform the user which ones are being removed
  if (length(invalidTraitNames) != 0) {
    warningMessage <- paste("Invalid trait names detected and removed:",
                            paste(invalidTraitNames, collapse=", "),
                            "Please check your trait names.")
    message(warningMessage)

    # Check if the cleaning process results in an empty list of traits
    if (length(validTraitNames) == 0) {
      stop("All provided trait names were invalid. Please check the variable names.")
    }
  }

  # Return the processed dataset
  # return(data)

  stats_descri <- data %>%
    select(!!sym(trial), all_of(traits)) %>%
    pivot_longer(cols = all_of(traits), names_to = "traits", values_to = "value") %>%
    group_by(!!sym(trial), traits) %>%
    summarise(
      Min = min(value, na.rm = TRUE),
      Mean = mean(value, na.rm = TRUE),
      Median = median(value, na.rm = TRUE),
      Max = max(value, na.rm = TRUE),
      SD = sd(value, na.rm = TRUE),
      CV = ifelse(Mean == 0, NA_real_, SD / Mean), # Handle division by zero
      n = sum(!is.na(value)),
      n_miss = sum(is.na(value)),
      miss_perc = (n_miss / n),
      .groups = "drop"
    )


  # Create heatmap with annotations
  p1 <- ggplot(stats_descri, aes(x = traits , y = trial_name, label = round(miss_perc,2),  fill = miss_perc ))+
    geom_tile(color = "gray")+
    geom_text(color = "white")+
    theme_minimal(base_size = 13)+
    labs(title = "Missing data", x = "", y = "") +
    theme(axis.text.x = element_text(hjust = 1 , angle = 75, size = 16),
          axis.text.y = element_text(size = 16))

  p1 <- plotly::ggplotly(p1)

  # Identify trials with high perc of missing data for specified traits
  missing_data_trials <- stats_descri %>%
    filter(miss_perc >= missing_data_threshold) %>% # Combine filter conditions for efficiency
    dplyr::distinct(trial_name, traits, miss_perc)

  if(nrow(missing_data_trials) == 0){
    message(paste0("No traits with more than ", (missing_data_threshold)*100, "%",  " of missing data"))
    data_filtered <- data
  }else{
    message(paste0("The following traits have been removed from the corresponding trial"))
    data_filtered <- data %>%
      mutate(across(all_of(unique(missing_data_trials$traits)), ~ifelse(trial_name %in% unique(missing_data_trials$trial_name), NA, .)))

  }


  objt <- list(
    stats_descri = stats_descri,
    missing_data_plot = p1,
    missing_data_trials_removed = missing_data_trials,
    data_cleaned = data_filtered,
    inputs = list(
      genotype = genotype,
      trial = trial,
      traits = traits,
      rep = rep,
      block = block,
      row = row,
      col = col
    )
  )

  return(objt)

}

