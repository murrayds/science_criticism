library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

letters <- read_csv(snakemake@input[[1]], col_types = cols())
traj <- read_csv(snakemake@input[[2]], col_types = cols())

# For each paper, calculate citations recieved each year...
paper_cites <- traj %>%
  collapse_aps() %>%
  group_by(id, type, citing_year, venue) %>%
  summarize(
    cites = length(unique(citing_id)),
    .groups = "drop"
  )

# Link each letter to its original
letter_cites <- paper_cites %>%
  filter(type == "letter") %>%
  left_join(
    letters %>% select(letter_id, lag, original_id), by = c("id"="letter_id")
  ) %>%
  select(-type) %>%
  rename(cites.letter = cites)

# Construct plotdata. Merge letters with original, and then copmute ratio of
# citations received by each.
plotdata <- paper_cites %>%
  filter(type == "original") %>%
  rename(cites.original = cites) %>%
  left_join(
    letter_cites, by = c("id" = "original_id", "citing_year", "venue")
  ) %>%
  group_by(id, venue) %>%
  tidyr::fill(lag) %>%
  mutate(cites.letter = ifelse(is.na(cites.letter), 0, cites.letter)) %>%
  group_by(id, venue) %>%
  mutate(
    citing_year_norm = citing_year - min(citing_year),
    cites.frac = cites.letter / cites.original,
  ) %>%
  group_by(citing_year_norm, venue) %>%
  summarize(
    mu = mean(cites.frac),
    n = n()
  ) %>%
  filter(citing_year_norm <= 8)

# construct the plot
p <- plotdata %>%
  ggplot(aes(x = citing_year_norm, y = mu, color = venue)) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = venue_colors()
  ) +
  scale_x_continuous(breaks = c(0:10)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_criticism() +
  theme(
    legend.position = c(0.85, 0.85)
  ) +
  xlab("Years since publication of critical letter") +
  ylab("impact(letter) / impact(original)")

ggsave(p, filename = snakemake@output[[1]], width = 7, height = 4, bg = "white")