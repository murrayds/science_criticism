#' Plot Impact Features of Targeted Papers and Associated Critical Letters
#'
#' This script reads in data about targeted papers and their associated
#' critical letters, processes the data to calculate citation trajectories
#' and ratios, and then generates a plot showing these features. The plot
#' provides a visual comparison of the impact between the targeted papers
#' and their associated critical letters. The script uses several R packages
#' including readr, dplyr, and ggplot2, and also sources in custom themes and
#'  common functions from other scripts.
#'
#' @input A CSV file containing data about the targeted papers and their
#' associated critical letters.
#' @input A CSV file containing citation trajectories for the letter and
#' original papers.
#' @input A CSV file containing titles data to identify those that mention
#' a reply.
#' @output A plot saved as a PNG file.
#'
#' @examples
#' Rscript scripts/plotting/plot_cite_ratio.R
#'

library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

letters <- read_csv(snakemake@input[[1]], col_types = cols())

# Get the titles data so we can identify those that mentiona reply
replies <- read_csv(snakemake@input[[3]], col_types = cols()) %>%
  filter(grepl("reply", citing_title, fixed = TRUE)) %>%
  select(CitingPaperId)

# Load the citation trajectories for the letter & original papers...
dual <- read_csv(snakemake@input[[2]], col_types = cols()) %>%
  collapse_aps() %>%
  # remove replies
  anti_join(replies, by = c("citing_id" = "CitingPaperId")) %>%
  # remove self citations
  filter(!is_self_cite)

# For each paper, calculate citations recieved each year...
paper_cites <- dual %>%
  group_by(id, type, citing_year, venue) %>%
  summarize(
    cites = length(unique(citing_id)),
    .groups = "drop"
  )

# Link each letter to its original
letter_cites <- paper_cites %>%
  filter(type == "letter") %>%
  left_join(
    letters %>%
      select(letter_id, lag, original_id, letter_year),
    by = c("id" = "letter_id")
  ) %>%
  select(-type) %>%
  rename(cites.letter = cites)

# Construct plotdata. Merge letters with original, and then copmute ratio of
# citations received by each.
plotdata_ratio <- paper_cites %>%
  filter(type == "original") %>%
  rename(cites.original = cites) %>%
  left_join(
    letter_cites, by = c("id" = "original_id", "citing_year", "venue")
  ) %>%
  group_by(id, venue) %>%
  tidyr::fill(lag) %>%
  tidyr::fill(letter_year) %>%
  mutate(cites.letter = ifelse(is.na(cites.letter), 0, cites.letter)) %>%
  group_by(id, venue) %>%
  filter(citing_year >= letter_year) %>%
  mutate(
    delta_year = citing_year - min(citing_year), # years since publication...
    cites.frac = cites.letter / cites.original
  ) %>%
  group_by(venue, delta_year) %>%
  summarize(
    mu = mean(cites.frac)
  ) %>%
  mutate(
    metric = "ratio"
  )


# This takes some time to compute...
# Identify those papers that co-cite the target (original) paper along with
# the associated critical letter.
#
# There is likely a better way to do this, the current way I have of
# associating the letter with the original is a little hacky.
cocite <- dual %>%
  left_join(
    letters %>%
      select(original_id, letter_id) %>%
      distinct(original_id, letter_id, .keep_all = TRUE),
    by = c("id" = "letter_id")
  ) %>%
  rowwise() %>%
  mutate(group_id = coalesce(original_id, id)) %>%
  group_by(venue, citing_id, citing_year, group_id) %>%
  summarize(
    cocite = n() > 1,
    .groups = "drop"
  ) %>%
  select(venue, citing_id, group_id, citing_year, cocite)

# Collect the annual % of citations to the target (origiainL) paper
# that also cite the critical letter
plotdata_cocite <- cocite %>%
  left_join(
    letters %>% select(original_id, letter_year),
    by = c("group_id" = "original_id")
  ) %>%
  mutate(delta_year = citing_year - letter_year) %>%
  group_by(venue, group_id, delta_year) %>%
  summarize(
    total = n(),
    prop = sum(cocite) / n()
  ) %>%
  group_by(venue, delta_year) %>%
  summarize(
    mu = mean(prop)
  ) %>%
  mutate(metric = "cocite")

# Collect the yearly raw annual citations accumulated by the
# targeted (original) paper
plotdata_rawcite <- dual %>%
  filter(type == "original") %>%
  left_join(
    letters %>% select(original_id, letter_year),
    by = c("id" = "original_id")
  )  %>%
  mutate(delta_year = citing_year - letter_year) %>%
  group_by(venue, id, delta_year) %>%
  summarize(
    citations = n()
  ) %>%
  group_by(venue, delta_year) %>%
  summarize(
    mu = (mean(citations))
  ) %>%
  mutate(metric = "raw")

# Combine all the plots together
plotdata <- data.table::rbindlist(
  list(
    plotdata_rawcite,
    plotdata_cocite
  ),
  use.names = TRUE
)

# Do some processing on the names and such...
plotdata_final <- plotdata %>%
  mutate(
    metric = factor(
      metric,
      levels = c("raw", "cocite"),
      labels = c(
        "(A)",
        "(B)"
      )
    )
  ) %>%
  filter(
    delta_year >= 1,
    delta_year <= 8
  )


plotdata_slopes <- plotdata_final %>%
  group_by(venue, metric) %>%
  arrange(delta_year) %>%
  mutate(slope = (mu - lag(mu)) / (delta_year - lag(delta_year))) %>%
  summarize(
    avg_slope = round(mean(slope, na.rm = TRUE), 3),
    slope_lab = ifelse(
      avg_slope < 0,
      paste0("m = ", avg_slope),
      paste0("m = +", avg_slope)
    ),
    last_val = last(mu),
  )

# Build the plot
p <- plotdata_final %>%
  ggplot(aes(x = delta_year, y = mu, color = venue)) +
  geom_point() +
  geom_line() +
  # geom_text(
  #   data = plotdata_slopes,
  #   aes(
  #     x = 0, y = last_val,
  #     label = slope_lab
  #   ),
  #   hjust = 1,
  #   vjust = 1,
  #   size = 3
  # ) +
  facet_grid(metric ~ venue, scale = "free_y", switch = "y") +
  scale_color_manual(values = venue_colors()) +
  scale_y_continuous(
    position = "right"
  ) +
  scale_x_continuous(
    breaks = c(0, 2, 4, 6, 8)
  ) +
  theme_criticism() +
  theme(
    legend.position = "none",
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 0, vjust = 1,
      size = 12,
      face = "bold"
    ),
    strip.text.x = element_text(hjust = 0, size = 12),
    axis.title.y = element_blank(),
    panel.spacing.x = unit(0.30, "cm", data = NULL)
  ) +
  xlab("Years since receipt of critical letter")

# Save the output
ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 7.5,
  height = 3.0,
  bg = "white"
)