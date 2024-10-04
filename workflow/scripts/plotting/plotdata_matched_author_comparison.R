library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

source("scripts/common.R")
source("scripts/plotting/theme.R")

df <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  collapse_aps()

# Either `frac_prod`` or `impact_raw``
wc_authorship <- snakemake@wildcards[[1]]
variable_name <- snakemake@wildcards[[2]]
if (variable_name == "citations") {
  variable_name <- "impact_raw"
} else if (variable_name == "prod") {
  variable_name <- "frac_prod"
}

df_matched <- df %>%
  rename(var = variable_name) %>%
  rowwise() %>%
  mutate(venue = factor(venue, levels = rev(venue_levels()))) %>%
  select(AuthorId, venue, treat, match.group, partition, var) %>%
  pivot_wider(names_from = partition, values_from = var) %>%
  mutate(change = after / before) %>%
  filter(!is.na(change)) %>%
  group_by(venue, match.group) %>%
  filter(n() > 1) %>%
  ungroup()


# Construct the data to be plotted
plotdata <- df_matched %>%
  group_by(venue, match.group) %>%
  arrange(venue, match.group, treat) %>%
  summarize(
    ratio = last(change) - first(change)
  ) %>%
  group_by(venue) %>%
  summarize(
    mean_diff = mean(ratio),
    interval = 1.96 * sd(ratio) / sqrt(n()),
    upper = mean_diff + interval,
    lower = mean_diff - interval,
  )

# Perform statistical test...
tests <- df_matched %>%
  mutate(venue = factor(venue)) %>%
  group_by(venue) %>%
  arrange(match.group, treat) %>%
  do(
    t = t.test(
      (change) ~ treat,
      data = .,
      paired = TRUE
    )
  ) %>%
  summarize(
    venue = venue,
    t.p.value = round(t$p.value, 4),
  )

plotdata <- plotdata %>%
  left_join(tests, by = "venue") %>%
  mutate(
    authorship = wc_authorship,
    metric = variable_name
  )

write_csv(plotdata, snakemake@output[[1]])