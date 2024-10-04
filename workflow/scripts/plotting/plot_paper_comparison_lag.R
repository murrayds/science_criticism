library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

papers <- read_csv(snakemake@input[[1]], col_types = cols())

plotdata <- papers %>%
  mutate(
    signif = case_when(
      t.p.value > 0.05 ~ "",
      t.p.value > 0.01 ~ "*",
      t.p.value > 0.001 ~ "**",
      !is.na(t.p.value) ~ "***",
      TRUE ~ NA_character_
    ),
    venue = factor(venue, levels = rev(venue_levels())),
    lagtype = factor(
      lagtype,
      levels = c("lagall", "lag0", "lag1plus"),
      labels = c(
        "All",
        "Lag = 0",
        "Lag > 0"
      )
    )
  )

p <- plotdata %>%
  ggplot(aes(x = mean_diff, y = venue, color = venue)) +
  geom_point(size = 4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(
    aes(label = signif, x = upper),
    nudge_x = 0.075,
    nudge_y = -0.06,
    size = 6,
  ) +
  facet_wrap(~lagtype, scale = "free_x", nrow = 1) +
  scale_color_manual(values = venue_colors(), guide = "none") +
  theme_criticism() +
  theme(
    legend.position = c(0.25, 0.75),
    panel.spacing.x = unit(0.30, "cm", data = NULL),
    strip.text = element_text(angle = 0, hjust = 0),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "bold")
  ) +
  xlab("ΔTreatment - ΔControl")

# Save the figure output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 6,
  height = 3,
  bg = "white"
)
