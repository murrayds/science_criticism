library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df_matched <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  collapse_aps() %>%
  mutate(
    treat = factor(type, levels = c("Treatment", "Control")),
    venue = factor(venue, levels = venue_levels())
  ) %>%
  filter(partition == "before") %>%
  mutate(
    treat = factor(treat, levels = c(1, 0), labels = c("Treatment", "Control")),
    venue = factor(venue, levels = venue_levels())
  ) %>%
  arrange(match.group) %>%
  group_by(match.group) %>%
  filter(n() >= 2)


wc_authorship <- snakemake@wildcards[[1]]
variable_name <- snakemake@wildcards[[2]]
if (variable_name == "citations") {
  variable_name <- "impact_raw"
  lab_xpos <- 2
  lab_ypos <- 0.5
  xlab <- "Average impact prior to critical letter"
} else if (variable_name == "prod") {
  variable_name <- "frac_prod"
  lab_xpos <- 0.05
  lab_ypos <- 1.5
  xlab <- "Average fractional productivity prior to critical letter"
}

df_matched <- df_matched %>%
  rename(metric = variable_name)

plotlabs <- df_matched %>%
  group_by(venue, treat) %>%
  summarize(
    mu = round(mean(metric), 3),
    n = n(),
    label = paste0("μ=", mu, "\nn=", n)
  )

p <- df_matched %>%
  ggplot(aes(x = metric, color = venue, linetype = treat)) +
  geom_density(alpha = 0.5) +
  geom_label(
    data = plotlabs,
    aes(label = label),
    x = lab_xpos, y = lab_ypos,
    hjust = 0,
    alpha = 0.8
  ) +
  facet_grid(treat ~ venue) +
  scale_color_manual(values = venue_colors()) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_log10() +
  theme_criticism() +
  theme(
    legend.position = "none"
  ) +
  ylab("Density") +
  xlab(xlab)

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 8,
  height = 4,
  bg = "white"
)