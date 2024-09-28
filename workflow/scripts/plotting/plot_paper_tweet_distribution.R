library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]]) %>%
  collapse_aps()

features <- read_csv(snakemake@input[[3]], col_types = cols())

plotdata <- features %>%
  inner_join(df, by = "id") %>%
  collapse_aps() %>%
  filter(Tweet_Count > 0) %>%
  group_by(venue) %>%
  mutate(
    top_10pct = Tweet_Count >= quantile(Tweet_Count, prob = 0.90),
    top_05pct = Tweet_Count >= quantile(Tweet_Count, prob = 0.95),
    top_01pct = Tweet_Count >= quantile(Tweet_Count, prob = 0.99)
  ) %>%
  pivot_longer(
    cols = c(top_10pct, top_05pct, top_01pct),
    names_to = "metric"
) %>%
  group_by(venue, type, metric) %>%
  summarize(
    prop_hit = round(sum(value) / n(), 3)
  ) %>%
  filter(type == "letter") %>%
  mutate(
    metric = factor(
      metric,
      levels = c("top_10pct", "top_05pct", "top_01pct"),
      labels = c(
        "Among top 10%",
        "Among top 5%",
        "Among top 1%"
      )
    ),
    venue = factor(venue, levels = rev(venue_levels()))
  )

lines <- data.frame(
  "metric" = unique(plotdata$metric),
  "intercept" = c(0.01, 0.05, 0.1)
)


p <- plotdata %>%
  ggplot(aes(x = prop_hit, y = venue, fill = venue)) +
  geom_col(position = position_dodge(0.8)) +
  geom_text(
    aes(
      label = formatC(round(prop_hit, 3), digits = 3, format = "f"),
      x = ifelse(
        prop_hit < 0.08,
        prop_hit + 0.06, prop_hit - 0.05
      )
    ),
    position = position_dodge(0.8),
    size = 3.5
  ) +
  geom_vline(data = lines, aes(xintercept = intercept), linetype = "solid") +
  facet_wrap(~metric) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.31)) +
  scale_fill_manual(values = venue_colors()) +
  theme_criticism() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.spacing = unit(1.5, "lines")
  ) +
  xlab("Proportion among papers that received a critical letter")

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 6,
  height = 4,
  bg = "white"
)