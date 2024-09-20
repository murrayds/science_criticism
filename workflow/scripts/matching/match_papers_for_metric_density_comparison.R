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
cite_tolerance <- as.numeric(snakemake@wildcards[[2]])
year_tolerance <- as.numeric(snakemake@wildcards[[3]])

# Apply basic filtering and pre-processing first...
letters <- letters %>%
  filter(
    !is.na(month)
  ) %>%
  select(-id, -letter_id, -letter_year) %>%
  rename(id = original_id, year = original_year) %>%
  mutate(
    impact = impact_3year - num_critical_letters
  ) %>%
  select(-num_critical_letters) %>%
  mutate(type = "Treatment") %>%
  select(id, venue, type, year, month, impact)

articles <- articles %>%
  filter(
    !is.na(month)
  ) %>%
  rename(
    impact = impact_3year
  ) %>%
  mutate(type = "Control") %>%
  select(id, venue, type, year, month, impact)


df <- data.table::rbindlist(
  list(letters, articles),
  use.names = TRUE
) %>%
  mutate(
    year = as.integer(year),
    treat = ifelse(type == "Treatment", 1, 0),
    quarter = factor(paste0("Q", floor((month - 1) / 3) + 1))
  ) %>%
  # remove instances of no citations, since they don't give us any reliable
  # info for matching
  filter(impact >= 5)

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
    impact_norm = log(impact) / mean(log(impact))
  ) %>%
  # return only a subset of necessary columns
  select(
    id, venue, field, year, quarter, multifield, # basic metadata
    type, treat, # the type they are in (treat is just binary (1/0))
    impact, impact_norm # impact metrics
  )

# Perform the matching...
df_matched <- perform_matching(
  df_prepared,
  cite_tolerance = cite_tolerance,
  year_tolerance = year_tolerance
)

# perform final processing
df_matched <- df_matched %>%
  rowwise() %>%
  mutate(
    match.group = paste0(venue, ".", subclass) # set unique group ID
  ) %>%
  ungroup() %>%
  arrange(venue, match.group, type) %>%
  #filter(!duplicated(id)) %>% # remove duplicated records
  group_by(match.group) %>%
  filter(n() > 1) # keep only groups for which a match was identified

write_csv(df_matched, snakemake@output[[1]])
