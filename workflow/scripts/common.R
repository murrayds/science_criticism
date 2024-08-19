source("scripts/plotting/theme.R")

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
    select(original_id, original_year, venue, impact_2year, impact_3year, impact_4year, lag) %>%
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