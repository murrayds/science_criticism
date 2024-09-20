suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(MatchIt)

source("scripts/common.R")

#
# Get script parameters from Snakemake
#

# the target venue to match
target_venue <- snakemake@wildcards[[1]]

# Either "f"(irst), "m"(iddle), or "l"(ast) to mark which authorship position
# we are examining...
authorship <- substring(snakemake@wildcards[[2]], 1, 1)
cite_tolerance <- as.numeric(snakemake@wildcards[[3]])
prod_tolerance <- as.numeric(snakemake@wildcards[[4]])


# Load the datasets and filter to the relevant venue
letters <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  filter(venue == target_venue)
articles <- read_csv(snakemake@input[[2]], col_types = cols()) %>%
  filter(venue == target_venue)
histories <- read_csv(snakemake@input[[3]], col_types = cols())

#
# Define other constants
#
kStartYear <- 2000
kEndYear <- 2018

kWindowAroundLetter <- 4
kPublicationWindow <- 5

kMinProductivity <- 5
kMinCareerAge <- 5

# Iterate over each year
df_matched_part <- lapply(c(kStartYear:kEndYear), function(event_year) {
  print(paste0("Starting for venue: ", target_venue, ", year: ", event_year))
  treatment_authors <- histories %>%
    inner_join(letters, by = c("PaperId" = "id")) %>%
    filter(
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
      career_age = first(career_age),
      lead_prod = sum(author_position %in% c("f", "l")) / (max(Year) - min(Year)),
      # for productivity, divide by total years in range to get an average,
      # the number of years may be unbalanced so we want to smooth it out
      frac_prod = sum(1 / num_authors) / (max(Year) - min(Year)),
      impact = mean(impact_3year_norm, na.rm = TRUE),
      impact_raw = mean(impact_3year_norm, na.rm = TRUE),
      field = first(DescTools::Mode(field, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    filter(total_prod > kMinProductivity) %>%
    mutate(
      # in case the impact is zero, replace with a very small number
      impact = ifelse(impact == 0, 0.00001, impact),
      impact_norm = log(impact) / log(mean(impact, na.rm = TRUE)),
      prod_norm = log(frac_prod) / log(mean(frac_prod, na.rm = TRUE)),
      lead_prod_norm = log(lead_prod) / log(mean(lead_prod, na.rm = TRUE)),
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
  formula <- treat ~ prod_norm + lead_prod_norm + impact_norm + field + career_age

  # Perform the matching
  match_result <- tryCatch({
    matchit(
      formula,
      data = prepared_df %>% filter(partition == "before"),
      method = "nearest",
      caliper = c(
        impact_norm = cite_tolerance,
        frac_prod = prod_tolerance,
        lead_prod_norm = 0.15
      ),
      ratio = 1, # 1 nearest match for each record
      discard = "both" # remove instances where a match could not be found
    )
  },
  error = function(cond) {
    # Throwing errors is common since we may not find matches.
    # In this case, return NULL
    message(conditionMessage(cond))
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
      match.group = paste0(target_venue, ".", event_year, ".", subclass),
      venue = target_venue,
      event_year = event_year,
    ) %>%
    select(AuthorId, venue, match.group, field) %>%
    left_join(prepared_df, by = c("AuthorId", "field"))

  return(matched_part)
}) # end event_year loop

df_matched <- data.table::rbindlist(df_matched_part)

write_csv(df_matched, snakemake@output[[1]])
