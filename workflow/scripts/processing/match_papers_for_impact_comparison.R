suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(MatchIt)

source("scripts/common.R")

letters <- read_csv(snakemake@input[[1]], col_types = cols())
articles <- read_csv(snakemake@input[[2]], col_types = cols())
fields <- read_csv(snakemake@input[[3]], col_types = cols())


# This is the amount of time POST CRITICAL LETTER at which we should measure
# the impact.
impact_delay <- as.numeric(snakemake@wildcards[[1]])
impact_before_lower_limit <- as.numeric(snakemake@wildcards[[2]])
cite_tolerance <- as.numeric(snakemake@wildcards[[3]])
year_tolerance <- as.numeric(snakemake@wildcards[[4]])

# Apply basic filtering and pre-processing first...
letters <- letters %>%
  filter(
    original_year <= (2021 - impact_delay),
    !is.na(month)
  ) %>%
  select(-id, -letter_id, -letter_year) %>%
  rename(id = original_id, year = original_year) %>%
  mutate(type = "Treatment")

articles <- articles %>%
  filter(
    year <= (2021 - impact_delay),
    !is.na(month)
  ) %>%
  mutate(type = "Control")

prepare_df_for_matching <- function(x, letters, articles, fields) {
  # This function will encapsulate the logic of preparing a single
  # dataframe for matching.
  #
  # We will perform the matching for each possible value of "lag",
  # which designates the time lag between the publication of the
  # critical letter and the original targeted publication

  letters_filt <- letters %>%
    filter(lag == x) %>%
    select(-lag) %>%
    rename(
      impact_before = paste0("impact_", ifelse(x == 0, 1, x), "year"),
      impact_after = paste0(
        "impact_",
        ifelse(x == 0, 1, x) + impact_delay,
        "year"
      )
    ) %>%
    mutate(
      impact_before = impact_before - num_critical_letters,
      impact_after = impact_after - num_critical_letters
    ) %>%
    select(-num_critical_letters)

  articles_filt <- articles %>%
    rename(
      impact_before = paste0(
        "impact_", ifelse(x == 0, 1, x), "year"
      ),
      impact_after = paste0(
        "impact_", ifelse(x == 0, 1, x) + impact_delay, "year"
      )
    )

  df <- data.table::rbindlist(
    list(letters_filt, articles_filt),
    use.names = TRUE
  ) %>%
    mutate(
      year = as.integer(year),
      treat = ifelse(type == "Treatment", 1, 0),
      quarter = factor(paste0("Q", floor((month - 1) / 3) + 1))
    ) %>%
    # remove instances of no citations, since they don't give us any reliable
    # info for matching
    filter(impact_before >= impact_before_lower_limit)

  # Now we will incorperate the field-level information...
  df_prepared <- fields %>% 
    inner_join(df, by = "id") %>%
    group_by(id) %>%
    filter(!duplicated(field_level0)) %>%
    mutate(
      multifield = n() > 1
    ) %>%
    sample_n(1) %>% # select only one level0 field
    ungroup() %>%
    select(-field) %>%
    rename(field = field_level0) %>%
    group_by(venue) %>%
    mutate(
      impact_norm = log(impact_before) / mean(log(impact_before))
    ) %>%
    # return only a subset of necessary columns
    select(
      id, venue, field, year, quarter, multifield, # basic metadata
      type, treat, # the type they are in (treat is just binary (1/0))
      impact_before, impact_after, impact_norm # impact metrics
    ) %>%
    mutate(lag = x)

  return(df_prepared)
}


perform_matching <- function(df) {
  # This function encapsulates the logic of performing the matching using
  # R's `MatchIt` package.

  # allow some flexibility on year, but focus on impact norm...
  formula <- treat ~ year + impact_norm

  # Must have exact match on field and venue...
  exact_formula <- ~ field + quarter + venue + multifield

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

# For values of `lag` between 1 and 10, prepare the data frame
# and perform the matching`
df_matched_list <- lapply(c(0:6), function(x) {
  prepared_df_part <- prepare_df_for_matching(x, letters, articles, fields)
  matched_df_part <- perform_matching(prepared_df_part)
  return(matched_df_part)
})

# aggregate the results as a single dataframe
df_matched <- data.table::rbindlist(
  df_matched_list,
  use.names = TRUE,
  fill = TRUE
)

# perform final processing
df_matched <- df_matched %>%
  rowwise() %>%
  mutate(
    growth = impact_after / impact_before,  # calculate the growth
    match.group = paste0(lag, ".", subclass) # set unique group ID
  ) %>%
  ungroup() %>%
  arrange(venue, match.group, type) %>%
  #filter(!duplicated(id)) %>% # remove duplicated records
  group_by(match.group) %>%
  filter(n() > 1) # keep only groups for which a match was identified

write_csv(df_matched, snakemake@output[[1]])
