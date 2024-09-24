library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Load the datasets and filter to the relevant venue
letters <- read_csv(snakemake@input[[1]], col_types = cols())

match_files <- snakemake@input[c(2:length(snakemake@input))]

# Either `frac_prod`` or `impact_raw``
metric_to_show <- snakemake@wildcards[[1]]
# Hold the count of letters for each venue which is needed to calculate the
# percentage of total matches
letter_counts <- letters %>%
  group_by(venue) %>%
  summarize(
    count_candidates = length(unique(letter_id)),
    .groups = "drop"
  )

# Now, we iterate through each file, calculating information about
# the parameters, matches made, and effects observed
matched <- data.table::rbindlist(lapply(match_files, function(f) {
  df <- read_csv(f, col_types = cols())

  authorship <- sub(".*_(first|last)_authors.*", "\\1", f)
  impact <- as.numeric(sub(".*_(\\d+\\.\\d+)impact_.*", "\\1", f))
  productivity <- as.numeric(sub(".*_(\\d+\\.\\d+)prod.*", "\\1", f))

  # extract the file parameters
  t_matched <- df %>%
    rowwise() %>%
    rename(`metric` = metric_to_show) %>%
    select(AuthorId, venue, treat, match.group, partition, metric) %>%
    pivot_wider(names_from = partition, values_from = metric) %>%
    mutate(change = after / before) %>%
    filter(!is.na(change)) %>%
    group_by(venue, match.group) %>%
    filter(n() > 1) %>%
    group_by(venue, match.group, treat) %>%
    arrange(venue, match.group, treat)

  print(head(t_matched))
  # Conduct statistical tests and save results
  tests <- t_matched %>%
    mutate(
      venue = factor(venue),
      treat = factor(treat)
    ) %>%
    group_by(venue) %>%
    arrange(match.group, treat) %>%
    do(
      w = wilcox.test(
        (change) ~ treat,
        data = .,
        paired = TRUE
      ),
      t = t.test(
        (change) ~ treat,
        data = .,
        paired = TRUE
      )
    ) %>%
    summarize(
      venue = venue,
      t.statistic = t$statistic,
      t.estimate = t$estimate,
      t.p.value = round(t$p.value, 4),
      wilcox.statistic = w$statistic,
      wilcox.p.value = round(w$p.value, 4),
    )

  t_matched %>%
    group_by(venue) %>%
    summarize(
      num_matched = length(unique(match.group)),
      authorship = authorship,
      impact = impact,
      productivity = productivity,
      .groups = "drop"
    ) %>%
    left_join(tests, by = "venue")
}))

matched <- matched %>%
  left_join(letter_counts, by = "venue")

write_csv(matched, snakemake@output[[1]])