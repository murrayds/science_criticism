library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


letters <- read_csv(snakemake@input[[1]], col_types = cols())
match_files <- snakemake@input[c(2:length(snakemake@input))]

# First, build a dataframe calculating the "impact before" 
# values for every letter
letters <- data.table::rbindlist(lapply(c(0:12), function(x) {
  letters %>%
    filter(lag == x) %>%
    rename(
      impact_before = paste0("impact_", ifelse(x == 0, 1, x), "year")
    ) %>%
    select(id, venue, original_year, month, lag, impact_before)
}))

# Now, we iterate through each file, calculating information about
# the parameters, matches made, and effects observed
matched <- data.table::rbindlist(lapply(match_files, function(f) {
  df <- read_csv(f, col_types = cols())

  delay <- as.numeric(sub(".*_(\\d+)delay_.*", "\\1", f))
  impact <- as.numeric(sub(".*_(\\d+\\.\\d+)impact_.*", "\\1", f))
  year <- as.numeric(sub(".*_(\\d+)year\\.csv$", "\\1", f))

  # extract the file parameters
  t_matched <- df %>%
    group_by(venue) %>%
    summarize(
      num_matched = sum(type == "Treatment")
    ) %>%
    mutate(
      delay = delay,
      impact = impact,
      year = year
    )

  # Conduct statistical tests and save results
  tests <- df %>%
    mutate(venue = factor(venue)) %>%
    group_by(venue) %>%
    arrange(match.group, type) %>%
    do(
      w = wilcox.test(
        (growth) ~ type,
        data = .,
        paired = TRUE
      ),
      t = t.test(
        (growth) ~ type,
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

  # Now examine which records were excluded due to 
  # filters, and which were not matched
  letters %>%
    mutate(
      year_oob = !(original_year <= (2021 - delay)),
      month_isna = is.na(month),
      lag_oob = lag > 6,
      low_impact = impact_before < 5
    ) %>%
    group_by(venue) %>%
    summarize(
      count_total = n(),
      count_candidates = count_total -
        sum(year_oob, month_isna, lag_oob, low_impact),
    ) %>%
    left_join(t_matched, by = "venue") %>%
    left_join(tests, by = "venue")
}))

write_csv(matched, snakemake@output[[1]])