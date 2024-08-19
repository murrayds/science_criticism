library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")

df_matched <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  mutate(
    treat = factor(type, levels = c("Treatment", "Control")),
    venue = factor(venue, levels = venue_levels())
  )

plotlabs <- df_matched %>%
  group_by(venue, treat) %>%
  summarize(
    mu = round(mean(impact_before), 1),
    n = n(),
    label = paste0("μ=", mu, "\nn=", n)
  )

p <- df_matched %>%
  ggplot(aes(x = impact_before, color = venue, linetype = treat)) +
  geom_density(alpha = 0.5) +
  geom_label(data = plotlabs, aes(x = 50, y = 1, label = label), hjust = 0) +
  facet_grid(treat ~ venue) +
  scale_color_manual(values = venue_colors()) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_log10() +
  theme_criticism() +
  theme(
    legend.position = "none"
  ) +
  ylab("Density") +
  xlab("Impact prior to critical letter")

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 7,
  height = 5,
  bg = "white"
)