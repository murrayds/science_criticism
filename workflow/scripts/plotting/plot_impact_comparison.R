library(readr)
library(xtable)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")

df_matched <- read_csv(snakemake@input[[1]], col_types = cols())

# First, perform statistical test...
tests <- df_matched %>%
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
    wilxoc.p.value = round(w$p.value, 4),
  )

plotdata <- df_matched %>%
  mutate(
    venue = factor(venue, levels = venue_levels()),
    type = factor(type, levels = c("Treatment", "Control"))
  ) %>%
  group_by(venue, type) %>%
  summarize(
    mu_growth = mean(growth),
    interval = sd(growth) / sqrt(n()),
    upper = mu_growth + interval,
    lower = mu_growth - interval
  ) 

# Format the test results that will be displayed on the plot...
plottests <- plotdata %>%
  group_by(venue) %>%
  summarize(
    type = "Treatment",
    mx = max(upper) 
  ) %>%
  left_join(tests, by = "venue") %>%
  mutate(
    signif = case_when(
          t.p.value > 0.05 ~ "",
          t.p.value > 0.01 ~ "*",
          t.p.value > 0.001 ~ "**",
          !is.na(t.p.value) ~ "***",
          TRUE ~ NA_character_
        )
  )

# Now construct the plot object
pooled_plot <- plotdata %>%
  ggplot(aes(x = venue, y = mu_growth, color = venue, shape = type)) +
  geom_point(position = position_dodge(width=0.5), size = 4) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper), 
    width=0, 
    position=position_dodge(width=0.5)
  ) +
  geom_text(
    data = plottests, 
    aes(label = signif, y = mx), 
    nudge_y = 0.10,
    size = 6
  ) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values = venue_colors(), guide=FALSE) +
  theme_criticism() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.9),
    axis.title.x = element_blank(),
  ) +
  ylab("Δ Impact")


pairwise_plot <- df_matched %>%
  arrange(venue, match.group, type) %>%
  group_by(venue, match.group) %>%
  summarize(
    diff = last(growth) - first(growth)
  ) %>%
  group_by(venue) %>%
  summarize(
    mean_diff = mean(diff),
    interval = 1.96 * sd(diff) / sqrt(n()),
    upper = mean_diff + interval,
    lower = mean_diff - interval
  ) %>%
  mutate(
    venue = factor(venue, levels = venue_levels()),
  ) %>%
  left_join(plottests, by = "venue") %>%
  ggplot(aes(x = venue, y = mean_diff, color = venue)) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper), 
    width=0
  ) +
  geom_text(
    aes(label = signif, y = upper), 
    nudge_y = 0.075,
    size = 6
  ) +
  scale_color_manual(values = venue_colors(), guide=FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_criticism() +
  theme(
    axis.title.x = element_blank()
  ) +
  ylab("ΔTreatment - ΔControl")


# Save the figure output
ggsave(
    pooled_plot,
    filename = snakemake@output[[1]],
    width = 4,
    height = 4,
    bg = "white"
)

# Save the figure output
ggsave(
    pairwise_plot,
    filename = snakemake@output[[2]],
    width = 4,
    height = 4,
    bg = "white"
)

# Construct the latex table for the test results...
latex_table <- xtable(
  tests,
  align = c("llrrrrr"),
  digits = 4
)

print(
  latex_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[3]]
)