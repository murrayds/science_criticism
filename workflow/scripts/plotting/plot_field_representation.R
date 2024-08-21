library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]])
fields <- read_csv(snakemake@input[[3]], col_types = cols())

# First, format the fields dataframe. Collapse to
# the level0 field, which is what we would like
# to focus on here.
df_fields <- fields %>%
  inner_join(df, by = "id") %>%
  group_by(id) %>%
  filter(!duplicated(field_level0)) %>%
  mutate(multifield = n() > 1) %>%
  ungroup() %>%
  select(-field) %>%
  rename(field = field_level0) %>%
  mutate(
    field = factor(
      field,
      levels = names(field_values()),
      labels = field_values()
    )
  )

# Work with non-letters and letters separately. Calculate the proportion of
# documents within each type that are associated with each field.
# Note, we are using full-counting, so papers appear in multiple fields...
baseline <- df_fields %>%
  filter(type == "article") %>%
  group_by(venue) %>%
  mutate(total = n()) %>%
  group_by(venue, field) %>%
  summarize(
    prop = n() / first(total),
    type = "Not targeted",
    count = n()
  )
# Lets just do an over/under-representation

targeted <- df_fields %>%
  filter(type == "letter") %>%
  group_by(venue) %>%
  mutate(total = n()) %>%
  group_by(venue, field) %>%
  summarize(
    prop = n() / first(total),
    type = "targeted",
    count = n()
  ) %>%
  filter(count > 10) # filter fields with a min number of documents


plotdata <- data.table::rbindlist(list(baseline, targeted)) %>%
  group_by(venue, field) %>%
  arrange(venue, field, type) %>%
  summarize(
    todrop = ifelse("targeted" %in% type, FALSE, TRUE),
    count = last(count),
    ratio = last(prop) / first(prop),
    ratio_sign = ifelse(ratio > 1, TRUE, FALSE)
  ) %>%
  filter(!todrop) # remove fields where letters are absent...

# Construct the plot...
p <- plotdata %>%
  ggplot(aes(x = ratio, y = field, fill = venue, alpha = ratio_sign)) +
  geom_col() +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_text(
    aes(
      label = paste0("n=", count),
      x = ifelse(ratio < 1, 1, ratio)
    ),
    size = 2.5, nudge_x = 0.2, hjust = 0
  ) +
  scale_x_continuous(limits = c(0, 7)) +
  scale_fill_manual(values = venue_colors(), guide = "none") +
  scale_alpha_manual(values = c(0.75, 1.0), guide = "none") +
  facet_wrap(~venue, nrow = 1) +
  theme_criticism() +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("(% in targets) / (% in non-targets)")


ggsave(p, filename = snakemake@output[[1]], width = 7, height = 4, bg = "white")