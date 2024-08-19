library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]])
fields <- read_csv(snakemake@input[[3]], col_types = cols())

# Break into 6 discrete time periods. Merge in level 1 (granular)
# fields. Calculate percentile rank of 2-year impact based on the
# distribution within each venue, time period, and field. When a paper
# is assigned multiple level 1 fields, take the average across all fields
df_ranked <- fields %>%
  filter(level == 1) %>%
  inner_join(df, by = "id") %>%
  mutate(year = as.integer(year)) %>%
  mutate(period = cut(year, breaks = 4)) %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(impact_4year)
  ) %>%
  group_by(venue, type, lag, id) %>%
  summarize(
    mu_rank = mean(rank)
  )

# Construct the data for the plot. Examine only the letters and see which 
# percentile bin they disproporinately fall within.
plotdata <- df_ranked %>%
  ungroup() %>%
  filter(type == "letter") %>%
  mutate(
    category = cut(
      mu_rank,
      breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
      include.lowest = TRUE
    )
  ) %>%
  group_by(venue) %>%
  mutate(total = n()) %>%
  group_by(venue, category) %>%
  summarize(
    prop = n() / first(total)
  ) %>%
  ungroup() %>%
  # We lose a few records when we filter NAs, its because 
  # they could not be matched to a level1 field in the time period
  filter(!is.na(category))


p <- plotdata %>%
  ggplot(aes(x = category, y = prop, fill = venue)) +
  geom_col(position = position_dodge(0.5), width = 0.5) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  scale_fill_manual(
    values = venue_colors()
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0.01, 0)) +
  theme_criticism() +
  theme(
    legend.position.inside = c(0.20, 0.75),
  ) +
  xlab("Percentile rank of targeted papers' impact") +
  ylab("%")

ggsave(p, filename = snakemake@output[[1]], width = 4, height = 4, bg = "white")
