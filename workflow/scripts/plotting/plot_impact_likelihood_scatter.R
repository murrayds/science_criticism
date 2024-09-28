library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/plotting/theme.R")
source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]]) %>%
  collapse_aps()

# Construct the plot data
plotdata <- df %>%
  # For each venue, split into logarithmically-sized bins using the
  # logbin function
  group_by(venue) %>%
  mutate(
    bin = cut(
      impact_3year,
      # scale is an adjustible parameter controling how fast the width
      # of the log-scaled bins grow. Tuning it has an effect on the
      # results so record the value carefully. The resulting plot
      # should therefore be used for guidance, only.
      breaks = logbin(df$impact_3year, scale = 1.2),
      ordered_result = TRUE,
      dig.lab = 1
    ),
    # While here, also determine the count of letters and articles
    num_letters = sum(type == "letter"),
    num_articles = sum(type == "article")
  ) %>%
  filter(!is.na(bin)) %>%
  # Now, identify the edge of each bin, which will be used for plotting
  rowwise() %>%
  mutate(
    binedge = as.numeric(
      unlist(
        strsplit(
          gsub(
            "(?![,.])[[:punct:]]", "",
            as.character(bin),
            perl = TRUE
          ),
          ","
        )
      )[2])
  ) %>%
  # For each bin/venue, estimate the probability of a letter
  group_by(binedge, venue) %>%
  summarize(
    p_with_letter = sum(type == "letter") / first(num_letters),
    p_no_letter = sum(type == "article") / first(num_articles),
    p_ratio = p_with_letter / p_no_letter,
    count = n()
  ) %>%
  filter(p_ratio != 0 & !is.infinite(p_ratio))

# calculate correlation between the binedges and the likelihood ratio.
# We use log-transformed values to stay true to the same representation
# as the plot.
stats <- (plotdata %>%
  group_by(venue) %>%
  summarize(
    rho = cor(log10(p_ratio), log10(binedge), method = "pearson")
  ) %>%
  mutate(
    label = sprintf("%s, ρ = %.2f", venue, round(rho, 2))
  ))$label

# construct the plot
p <- plotdata %>%
  ggplot(aes(x = binedge, y = p_ratio, color = venue)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(
    stat = "smooth",
    formula = y ~ x,
    method = "loess",
    se = FALSE,
    alpha = 0.75,
    linewidth = 0.65
  ) +
  scale_color_manual(labels = stats, values = venue_colors()) +
  facet_wrap(~venue) +
  scale_y_log10() +
  scale_x_log10() +
  theme_criticism() +
  theme(
    legend.position.inside = c(0.85, 0.22)
  ) +
  ylab("r(k)") +
  xlab("3-year citation impact")

ggsave(p, filename = snakemake@output[[1]], width = 6, height = 4, bg = "white")
