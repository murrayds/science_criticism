library(readr)
suppressPackageStartupMessages(library(dplyr))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df_matched <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  collapse_aps()


# We will assemble data for each lag parameter...
lag0 <- df_matched %>%
  filter(lag == 0) %>%
  mutate(lagtype = "lag0")

lag1plus <- df_matched %>%
  filter(lag > 0) %>%
  mutate(lagtype = "lag1plus")

df_matched$lagtype <- "lagall"

df_matched <- data.table::rbindlist(list(df_matched, lag1plus, lag0))

# First, perform statistical test...
tests <- df_matched %>%
  mutate(venue = factor(venue)) %>%
  group_by(venue, lagtype) %>%
  arrange(match.group, type) %>%
  do(
    t = t.test(
      (growth) ~ type,
      data = .,
      paired = TRUE
    )
  ) %>%
  summarize(
    venue = venue,
    lagtype = lagtype,
    t.p.value = round(t$p.value, 4),
  )

plotdata <- df_matched %>%
  mutate(
    venue = factor(venue, levels = rev(venue_levels())),
    type = factor(type, levels = c("Treatment", "Control"))
  ) %>%
  group_by(venue, type, lagtype) %>%
  summarize(
    mu_growth = mean(growth),
    interval = sd(growth) / sqrt(n()),
    upper = mu_growth + interval,
    lower = mu_growth - interval
  )

plotdata <- df_matched %>%
  arrange(venue, match.group, type) %>%
  group_by(venue, lagtype, match.group) %>%
  summarize(
    diff = last(growth) - first(growth),
    .groups = "drop"
  ) %>%
  group_by(venue, lagtype) %>%
  summarize(
    mean_diff = mean(diff),
    interval = 1.96 * sd(diff) / sqrt(n()),
    upper = mean_diff + interval,
    lower = mean_diff - interval,
    .groups = "drop"
  ) %>%
  mutate(
    venue = factor(venue, levels = rev(venue_levels())),
  ) %>%
  left_join(tests, by = c("venue", "lagtype"))

write_csv(plotdata, snakemake@output[[1]])