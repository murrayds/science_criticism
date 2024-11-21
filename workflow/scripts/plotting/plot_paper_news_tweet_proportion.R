library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

matched <- read_csv(snakemake@input[[1]], col_types = cols())
features <- read_csv(snakemake@input[[2]], col_types = cols())


plotdata <- matched %>%
  collapse_aps() %>%
  left_join(features, by = "id") %>%
  group_by(venue) %>%
  mutate(
    has_tweet = Tweet_Count > 0,
    has_news = Newsfeed_Count > 0,
    type = factor(type, levels = c("Control", "Treatment"))
  ) %>%
  select(id, venue, type, has_tweet, has_news) %>%
  pivot_longer(cols = c(has_tweet, has_news), names_to = "metric") %>%
  group_by(venue, type, metric) %>%
  summarize(
    prop = sum(value, na.rm = TRUE) / n(),
    .groups = "drop"
  ) %>%
  mutate(
    venue = factor(venue, levels = rev(venue_levels())),
    metric = factor(
      metric,
      levels = c("has_tweet", "has_news"),
      labels = c("(A) At least one tweet", "(B) At least one news mention")
    )
  )

p <- plotdata %>%
  ggplot(aes(x = prop, y = venue, fill = venue, alpha = type)) +
  geom_col(position = position_dodge(0.8), color = "black") +
  geom_text(
    aes(
      label = formatC(round(prop, 3), digits = 3, format = "f"),
      x = ifelse(
        prop < 0.08,
        prop + 0.05, prop - 0.05
      )
    ),
    position = position_dodge(0.8),
    size = 3.5
  ) +
  facet_wrap(~metric) +
  scale_alpha_manual(values = c(0.6, 1.0), guide = "none") +
  scale_fill_manual(values = venue_colors(), guide = "none") +
  theme_criticism() +
  theme(
    axis.title.y = element_blank(),
    panel.spacing = unit(1.5, "lines")
  ) +
  xlab("Proportion of all papers")

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 6,
  height = 4,
  bg = "white"
)