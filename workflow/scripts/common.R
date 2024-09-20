source("scripts/plotting/theme.R")

library(MatchIt)

logbin <- function(data, scale = 1, zeros = FALSE) {
  # Generate breakpoints of exponentially increasing width
  # for use in binning logarithmic data.
  if (scale < 1) {
    stop("Function requires scale >= 1.")
  }

  count <- table(data)
  smax <- max(data)

  jmax <- ceiling(log(smax) / log(scale))
  if (zeros) {
    binedges <- scale^(0:jmax)
    binedges[1] <- 0
  } else {
    binedges <- scale^(1:jmax)
  }
  binedges <- unique(as.integer(binedges))
  x <- sqrt(binedges[-length(binedges)] * (binedges[-1] - 1))
  return(x)
}

load_aggregate_df <- function(letter_path, nonletter_path) {
  library(dplyr)
  library(readr)

  letters <- read_csv(letter_path, col_types = cols()) %>%
    select(
      original_id, original_year, venue,
      impact_2year, impact_3year, impact_4year, lag
    ) %>%
    rename(id = original_id, year = original_year) %>%
    mutate(type = "letter")

  nonletters <- read_csv(nonletter_path, col_types = cols()) %>%
    select(id, year, venue, impact_2year, impact_3year, impact_4year) %>%
    mutate(type = "article")

  df <- data.table::rbindlist(
    list(letters, nonletters),
    use.names = TRUE,
    fill = TRUE
  ) %>%
    mutate(
      venue = factor(venue, levels = venue_levels()),
      type = factor(type)
    )

  return(df)
}

perform_matching <- function(df, cite_tolerance, year_tolerance) {
  # This function encapsulates the logic of performing the matching using
  # R's `MatchIt` package.

  # allow some flexibility on year, but focus on impact norm...
  formula <- treat ~ year + impact_norm

  # Must have exact match on field and venue...
  exact_formula <- ~ field + quarter + multifield

  # Perform the matching
  match <- matchit(
    formula,
    data = df,
    method = "nearest",
    exact = exact_formula,
    caliper = c(
      impact_norm = cite_tolerance, # allow 5% variance in log-scaled citations
      year = year_tolerance # allow 3 years flexibility in matching year...
    ),
    ratio = 1, # 1 nearest match for each record
    discard = "both" # remove instances where a match could not be found
  )

  # Extract the matched dataset
  return(match.data(match))
}