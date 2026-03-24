screen_sta_random_effects <- function(random_effects, dat, min_levels = NULL) {
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data.frame", call. = FALSE)
  }

  if (length(random_effects) == 0) {
    return(list(
      kept = character(0),
      dropped = character(0),
      summary = data.frame(
        effect = character(0),
        exists = logical(0),
        n_non_missing = integer(0),
        n_levels = integer(0),
        min_required = integer(0),
        passed_min_levels = logical(0),
        equivalent_to = character(0),
        kept = logical(0),
        reason = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  random_effects <- unique(as.character(random_effects))

  if (is.null(min_levels)) {
    min_levels <- setNames(rep(2L, length(random_effects)), random_effects)
  } else {
    if (is.list(min_levels)) {
      min_levels <- unlist(min_levels)
    }
    min_names <- names(min_levels)
    min_levels <- as.integer(min_levels)
    names(min_levels) <- min_names
  }

  get_min_required <- function(effect, min_levels) {
    if (!is.null(names(min_levels)) && effect %in% names(min_levels)) {
      return(as.integer(min_levels[[effect]]))
    }
    2L
  }

  normalize_effect <- function(x) {
    x <- as.character(x)
    x[trimws(x) == ""] <- NA_character_
    x
  }

  is_equivalent_partition <- function(x, y) {
    cc <- stats::complete.cases(x, y)
    if (sum(cc) == 0) {
      return(FALSE)
    }

    x <- normalize_effect(x[cc])
    y <- normalize_effect(y[cc])

    cc2 <- stats::complete.cases(x, y)
    if (sum(cc2) == 0) {
      return(FALSE)
    }

    x <- x[cc2]
    y <- y[cc2]

    # Same partition up to relabeling:
    # each level of x maps to exactly one level of y, and vice versa.
    x_to_y <- tapply(y, x, function(z) length(unique(z)))
    y_to_x <- tapply(x, y, function(z) length(unique(z)))

    all(x_to_y == 1L) && all(y_to_x == 1L)
  }

  is_identical_coding <- function(x, y) {
    cc <- stats::complete.cases(x, y)
    if (sum(cc) == 0) {
      return(FALSE)
    }

    x <- normalize_effect(x[cc])
    y <- normalize_effect(y[cc])

    cc2 <- stats::complete.cases(x, y)
    if (sum(cc2) == 0) {
      return(FALSE)
    }

    x <- x[cc2]
    y <- y[cc2]

    identical(x, y)
  }

  summary_list <- vector("list", length(random_effects))
  names(summary_list) <- random_effects

  kept_effects <- character(0)

  for (effect in random_effects) {
    exists <- effect %in% colnames(dat)

    if (!exists) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = FALSE,
        n_non_missing = 0L,
        n_levels = 0L,
        min_required = get_min_required(effect, min_levels),
        passed_min_levels = FALSE,
        equivalent_to = NA_character_,
        kept = FALSE,
        reason = "column not found in dataset",
        stringsAsFactors = FALSE
      )
      next
    }

    x <- normalize_effect(dat[[effect]])
    n_non_missing <- sum(!is.na(x))
    n_levels <- length(unique(stats::na.omit(x)))
    min_required <- get_min_required(effect, min_levels)
    passed_min <- n_levels >= min_required

    if (n_non_missing == 0L) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = TRUE,
        n_non_missing = 0L,
        n_levels = 0L,
        min_required = min_required,
        passed_min_levels = FALSE,
        equivalent_to = NA_character_,
        kept = FALSE,
        reason = "all values are missing",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!passed_min) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = TRUE,
        n_non_missing = n_non_missing,
        n_levels = n_levels,
        min_required = min_required,
        passed_min_levels = FALSE,
        equivalent_to = NA_character_,
        kept = FALSE,
        reason = paste0("fewer than ", min_required, " observed levels"),
        stringsAsFactors = FALSE
      )
      next
    }

    equivalent_to <- NA_character_
    reason <- "kept"

    if (length(kept_effects) > 0) {
      for (prev_effect in kept_effects) {
        y <- normalize_effect(dat[[prev_effect]])

        if (is_identical_coding(x, y)) {
          equivalent_to <- prev_effect
          reason <- paste0("identical coding to kept effect '", prev_effect, "'")
          break
        }

        if (is_equivalent_partition(x, y)) {
          equivalent_to <- prev_effect
          reason <- paste0("confounded with kept effect '", prev_effect, "' (same partition)")
          break
        }
      }
    }

    keep_this <- is.na(equivalent_to)

    if (keep_this) {
      kept_effects <- c(kept_effects, effect)
    }

    summary_list[[effect]] <- data.frame(
      effect = effect,
      exists = TRUE,
      n_non_missing = n_non_missing,
      n_levels = n_levels,
      min_required = min_required,
      passed_min_levels = TRUE,
      equivalent_to = equivalent_to,
      kept = keep_this,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  list(
    kept = kept_effects,
    dropped = setdiff(random_effects, kept_effects),
    summary = summary_df
  )
}


screen_sta_residual_terms <- function(residual_terms, dat, min_levels = NULL) {
  if (!is.data.frame(dat)) {
    stop("'dat' must be a data.frame", call. = FALSE)
  }

  residual_terms <- unique(as.character(residual_terms))

  if (length(residual_terms) == 0) {
    return(list(
      kept = character(0),
      dropped = character(0),
      use_residual = FALSE,
      summary = data.frame(
        effect = character(0),
        exists = logical(0),
        n_non_missing = integer(0),
        n_levels = integer(0),
        min_required = integer(0),
        kept = logical(0),
        reason = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  if (is.null(min_levels)) {
    min_levels <- setNames(rep(2L, length(residual_terms)), residual_terms)
  } else {
    if (is.list(min_levels)) {
      min_levels <- unlist(min_levels)
    }
    if (is.null(names(min_levels))) {
      stop("'min_levels' must be a named vector or named list.", call. = FALSE)
    }
    min_names <- names(min_levels)
    min_levels <- as.integer(min_levels)
    names(min_levels) <- min_names
  }

  get_min_required <- function(effect, min_levels) {
    if (effect %in% names(min_levels)) {
      return(as.integer(min_levels[[effect]]))
    }
    2L
  }

  normalize_effect <- function(x) {
    x <- as.character(x)
    x[trimws(x) == ""] <- NA_character_
    x
  }

  summary_list <- vector("list", length(residual_terms))
  names(summary_list) <- residual_terms

  kept_terms <- character(0)

  for (effect in residual_terms) {
    exists <- effect %in% colnames(dat)
    min_required <- get_min_required(effect, min_levels)

    if (!exists) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = FALSE,
        n_non_missing = 0L,
        n_levels = 0L,
        min_required = min_required,
        kept = FALSE,
        reason = "column not found in dataset",
        stringsAsFactors = FALSE
      )
      next
    }

    x <- normalize_effect(dat[[effect]])
    n_non_missing <- sum(!is.na(x))
    n_levels <- length(unique(stats::na.omit(x)))

    if (n_non_missing == 0L) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = TRUE,
        n_non_missing = 0L,
        n_levels = 0L,
        min_required = min_required,
        kept = FALSE,
        reason = "all values are missing",
        stringsAsFactors = FALSE
      )
      next
    }

    if (n_levels < min_required) {
      summary_list[[effect]] <- data.frame(
        effect = effect,
        exists = TRUE,
        n_non_missing = n_non_missing,
        n_levels = n_levels,
        min_required = min_required,
        kept = FALSE,
        reason = paste0("fewer than ", min_required, " observed levels"),
        stringsAsFactors = FALSE
      )
      next
    }

    kept_terms <- c(kept_terms, effect)

    summary_list[[effect]] <- data.frame(
      effect = effect,
      exists = TRUE,
      n_non_missing = n_non_missing,
      n_levels = n_levels,
      min_required = min_required,
      kept = TRUE,
      reason = "kept",
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  list(
    kept = kept_terms,
    dropped = setdiff(residual_terms, kept_terms),
    use_residual = length(kept_terms) == length(residual_terms),
    summary = summary_df
  )
}
