library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/common.R")
source("scripts/plotting/theme.R")

df <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  collapse_aps()

# Either `frac_prod`` or `impact_raw``
variable_name <- snakemake@wildcards[[1]]
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

# Format the results of the statistical tests for addition to the plot
plottests <- plotdata %>%
  group_by(venue) %>%
  summarize(
    treat = 1,
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

p <- plotdata %>%
  left_join(plottests, by = "venue") %>%
  ggplot(aes(x = mean_diff, y = venue, color = venue)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(
    aes(label = signif, x = upper),
    nudge_x = 0.075,
    nudge_y = -0.0575,
    size = 6,
  ) +
  scale_color_manual(values = venue_colors(), guide = "none") +
  theme_criticism() +
  theme(
    legend.position = c(0.25, 0.75),
    axis.title.x = element_blank()
  ) +
  ylab("ΔTreatment - ΔControl")

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 4,
  height = 4,
  bg = "white"
)