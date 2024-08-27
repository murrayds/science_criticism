suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(MatchIt)


source("scripts/common.R")

letters <- read_csv(snakemake@input[[1]], col_types = cols())
articles <- read_csv(snakemake@input[[2]], col_types = cols())
histories <- read_csv(snakemake@input[[3]], col_types = cols())

#
# Get script parameters from Snakemake
#

# Either "f"(irst), "m"(iddle), or "l"(ast) to mark which authorship position
# we are examining...
authorship <- substring(snakemake@wildcards[[1]], 1, 1)
cite_tolerance <- 0.15
prod_tolerance <- 1

#
# Define other constants
#
kStartYear <- 2000
kEndYear <- 2018

kWindowAroundLetter <- 4
kPublicationWindow <- 5

kMinProductivity <- 5
kMinCareerAge <- 5

venues <- c("Nature", "Science", "PNAS", "PRL")

# First, iterate over each venue
df_matched_parts <- lapply(venues, function(event_venue) {
  # Now iterate over each year
  all_matches <- lapply(c(kStartYear:kEndYear), function(event_year) {
    print(paste0("Starting for venue: ", event_venue, ", year: ", event_year))
    treatment_authors <- histories %>%
      inner_join(letters, by = c("PaperId" = "id")) %>%
      filter(
        venue == event_venue, # focus on a specific venue...
        letter_year == event_year # go 1 year at a time to perform matching...
      ) %>%
      # we are only interested in first/last authors
      filter(author_position == authorship) %>%
      select(AuthorId) %>%
      mutate(
        treat = 1,
        center_year = event_year
      )

    # If there are no authors to match this year, skip and return NULL
    if (count(treatment_authors) == 0) {
      return(NULL)
    }

    # Identify articles published in the same venue around the same time as
    # the critical letter
    articles_near_letter <- articles %>%
      filter(
        venue == event_venue, # filter to same venue...
        year >= event_year - kWindowAroundLetter, # filter to relevant years
        year <= event_year + kWindowAroundLetter
      ) %>%
      # Set the "center" to the publication of this paper...
      # Add 1 to be in line with a critical letter...
      mutate(center_year = year + 1)

    # Now use these related papers as a pool from which to extract the
    # candidate authors for matching
    control_candidates <- histories %>%
      inner_join(articles_near_letter, by = c("PaperId" = "id")) %>%
      filter(author_position == authorship) %>%
      select(AuthorId, center_year) %>%
      mutate(treat = 0) %>%
      # make sure that the treatment authors aren't included in the candidates
      filter(!AuthorId %in% treatment_authors$AuthorId)

    df_all <- data.table::rbindlist(
      list(treatment_authors, control_candidates),
      use.names = TRUE
    ) %>%
      mutate(
        treat = factor(treat)
      )

    # Now lets calculate features...
    # One thing that would help to filter...focus only on ~3-5 years before &
    # after...get average based on "recent"  #truncate career age perhaps???
    prepared_df <- histories %>%
      inner_join(df_all, by = "AuthorId", relationship = "many-to-many") %>%
      group_by(AuthorId) %>%
      mutate(
        career_age = center_year - min(Year)
      ) %>% # focus on the timespan around the comment...
      filter(
        Year >= (center_year - kPublicationWindow),
        Year <= (center_year + kPublicationWindow),
        # need people with a minimum career age...
        career_age > kMinCareerAge
      ) %>%
      mutate(
        partition = ifelse(Year <= center_year, "before", "after")
      ) %>%
      ungroup() %>%
      group_by(AuthorId, partition, treat) %>%
      summarize(
        total_prod = n(),
        num_years = length(unique(Year)),
        career_age = first(career_age),
        lead_prod = sum(author_position %in% c("f", "l")) / num_years,
        # for productivity, divide by total years in range to get an average,
        # the number of years may be unbalanced so we want to smooth it out
        frac_prod = sum(1 / num_authors) / num_years,
        impact = mean(log(impact_3year_norm), na.rm = TRUE),
        impact_raw = mean(impact_3year_norm, na.rm = TRUE),
        field = first(DescTools::Mode(field, na.rm = TRUE)),
      ) %>%
      select(-num_years) %>%
      ungroup() %>%
      filter(total_prod > kMinProductivity) %>%
      mutate(
        impact_norm = impact / mean(impact, na.rm = TRUE),
        field = factor(field),
      ) %>%
      filter(
        !is.na(impact_norm),
        !is.na(field)
      ) %>%
      group_by(AuthorId) %>%
      # make sure we still have a before/after for every author
      filter("before" %in% partition & "after" %in% partition) %>%
      ungroup()


    # allow some flexibility on year, but focus on impact norm...
    formula <- treat ~ frac_prod + lead_prod + impact_norm + field + career_age

    # Perform the matching
    match_result <- tryCatch({
      matchit(
        formula,
        data = prepared_df %>% filter(partition == "before"),
        method = "nearest",
        caliper = c(
          impact_norm = cite_tolerance,
          frac_prod = prod_tolerance,
          lead_prod = 1.0
        ),
        ratio = 1, # 1 nearest match for each record
        discard = "both" # remove instances where a match could not be found
      )
    },
    error = function(cond) {
      # Throwing errors is common since we may not find matches.
      # In this case, return NULL
      match_result <- NULL
      return(NULL)
    })

    if (is.null(match_result)) {
      print("No matches found...")
      return(NULL)
    }

    
    # Extract the matched dataset
    matched_part <- MatchIt::match.data(match_result) %>%
      group_by(subclass) %>%
      filter(n() > 1) %>% # keep only those where a match was made...
      # ignore unmatched duplicate records (usually for other fields)
      filter(!is.na(partition)) %>%
      ungroup() %>%
      mutate(
        match.group = paste0(event_venue, ".", event_year, ".", subclass),
        event_year = event_year,
        venue = event_venue
      ) %>%
      select(AuthorId, venue, match.group, field) %>%
      left_join(prepared_df, by = c("AuthorId", "field"))

    return(matched_part)
  }) # end event_year loop

  df_matched_part <- data.table::rbindlist(all_matches)
  return(df_matched_part)
}) # end event_venue loop

# merge the dataframes into one
df_matched <- data.table::rbindlist(df_matched_parts)

write_csv(df_matched, snakemake@output[[1]])
