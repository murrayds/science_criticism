library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df <- read_csv(snakemake@input[[1]], col_types = cols())


plotdata <- df %>%
  mutate(
    venue = factor(venue, levels = rev(venue_levels())),
    factor = factor(
      factor,
      levels = c("genderFemale", "senioritySenior", "eliteElite"),
      labels = c(
        "Inferred gender is woman",
        "Is a senior researcher",
        "Affiliated with top university"
      )
    )
  )

p <- plotdata %>%
  ggplot(aes(x = AME, y = venue, color = venue)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~factor) +
  theme_criticism() +
  scale_color_manual(values = venue_colors(), legend) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.spacing = unit(1.5, "lines")
  ) +
  xlab("Average marginal effect on probability of receipt of criticism")

ggsave(p, filename = snakemake@output[[1]], width = 7, height = 4, bg = "white")
