library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

papers <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  filter(lagtype == "lagall") %>%
  select(venue, mean_diff, upper, lower, t.p.value) %>%
  mutate(
    panel = "papers"
  )

authors <- rbind(
  read_csv(snakemake@input[[2]], col_types = cols()),
  read_csv(snakemake@input[[3]], col_types = cols()),
  read_csv(snakemake@input[[4]], col_types = cols()),
  read_csv(snakemake@input[[5]], col_types = cols())
) %>%
  mutate(
    panel = paste0(authorship, "-", metric)
  ) %>%
  select(venue, mean_diff, upper, lower, t.p.value, panel)

plotdata <- rbind(authors, papers) %>%
  mutate(
    signif = case_when(
      t.p.value > 0.05 ~ "",
      t.p.value > 0.01 ~ "*",
      t.p.value > 0.001 ~ "**",
      !is.na(t.p.value) ~ "***",
      TRUE ~ NA_character_
    ),
    venue = factor(venue, levels = rev(venue_levels())),
    panel = factor(
      panel,
      levels = c(
        "papers",
        "first-frac_prod",
        "last-frac_prod",
        "first-impact_raw",
        "last-impact_raw"
      ),
      labels = c(
        "Paper impact",
        "First author productivity",
        "Last author productivity",
        "First author impact",
        "Last author impact"
      )
    )
  )

p <- plotdata %>%
  ggplot(aes(x = mean_diff, y = venue, color = venue)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text(
    aes(label = signif, x = upper),
    nudge_x = 0.075,
    nudge_y = -0.1,
    size = 6,
  ) +
  facet_wrap(~panel, scale = "free_x", nrow = 1) +
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
  width = 10,
  height = 3.5,
  bg = "white"
)
