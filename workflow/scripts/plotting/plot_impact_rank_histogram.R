library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")

letters <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  select(original_id, original_year, venue, impact_2year) %>%
  rename(id = original_id, year=original_year) %>%
  mutate(type = "letter")

articles <- read_csv(snakemake@input[[2]], col_types = cols()) %>%
  select(id, year, venue, impact_2year) %>%
  mutate(type = "article")

fields <- read_csv(snakemake@input[[3]], col_types = cols())

df <- data.table::rbindlist(list(letters, articles)) %>%
  mutate(
    venue = factor(venue, levels = venue_levels())
  )

# Break into 6 discrete time periods. Merge in level 1 (granular)
# fields. Calculate percentile rank of 2-year impact based on the
# distribution within each venue, time period, and field. When a paper
# is assigned multiple level 1 fields, take the average across all fields
df_ranked <- df %>%
  mutate(year = as.integer(year)) %>%
  mutate(period = cut(year, breaks = 4)) %>%
  left_join(fields %>% filter(level == 1), by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(impact_2year)
  ) %>%
  group_by(id) %>%
  summarize(
    mu_rank = mean(rank),
    venue = first(venue),
    type = first(type)
  )

# Construct the data for the plot. Examine only the letters and see which 
# percentile bin they disproporinately fall within.
plotdata <- df_ranked %>%
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
  geom_hline(aes(yintercept = 0.2, linetype = "Evenly distributed")) +
  scale_linetype_manual(
    name = "",
    values = 2,
  ) +
  scale_fill_manual(
    values = venue_colors()
  ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_criticism() +
  theme(
    legend.position = c(0.30, 0.75),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black")
  ) +
  xlab("Percentile rank of targeted papers' impact") +
  ylab("%")

ggsave(p, filename = snakemake@output[[1]], width = 4, height = 4, bg = "white")
