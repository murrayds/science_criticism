library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/common.R")
source("scripts/plotting/theme.R")


letters <- read_csv(snakemake@input[[1]], col_types = cols())

replies <- read_csv(snakemake@input[[2]], col_types = cols()) %>%
  filter(grepl("reply", citing_title, fixed = TRUE)) %>%
  select(CitingPaperId)

traj <- read_csv(snakemake@input[[3]], col_types = cols()) %>%
  collapse_aps() %>%
  # remove replies
  anti_join(replies, by = c("citing_id" = "CitingPaperId")) %>%
  # remove self citations
  filter(!is_self_cite)


cocite <- traj %>%
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

emb <- read_csv(snakemake@input[[4]], col_types = cols()) %>%
  left_join(
    cocite,
    by = c("cited_id" = "group_id", "citing_id" = "citing_id")
  ) %>%
  filter(!is.na(cocite))

# Build the dataframe for plotting...
plotdata <- emb %>%
  group_by(cited_id) %>%
  mutate(
    rank = percent_rank(cosine_similarity)
  ) %>%
  filter(cocite == TRUE)

# Compute the Kolmogorov–Smirnov test statistics for display on the plot
plottests <- plotdata %>%
  filter(!is.na(rank)) %>%
  # KS test complains if there are exact ties. Because ties are innevitable
  # in the data, we add a tiny jutter to all values so they are not exactly
  # the same.
  mutate(rank = jitter(rank, factor = 1e-8)) %>%
  group_by(venue) %>%
  do(
    # KS test, compared to uniform distribution between 0 and 1
    ks = ks.test(.$rank, "punif", min = 0, max = 1, alternative = "less"),
    mu = mean(.$rank) * 100,
    n = length(unique(.$citing_id))
  ) %>%
  summarize(
    ks.1s.statistic = ks$statistic,
    ks.1s.p.value = ks$p.value,
    mu = mu,
    n = n,
    venue = venue
  ) %>%
  rowwise() %>%
  mutate(
    ks.1s.statistic = round(ks.1s.statistic, 4),
    ks.1s.p.value = round(ks.1s.p.value, 4)
  ) %>%
  rowwise() %>%
  mutate(
    label_mean = paste0(
      "μ = ",
      formatC(round(last(mu), 1), digits = 1, format = "f"),
      "%"
    ),
    label_tests = paste0(
      "1s KS, p = ",
      formatC(round(ks.1s.p.value, 2), digits = 2, format = "f")
    )
  )

# Now we construct the plot...
p <- plotdata %>%
  ggplot(aes(x = rank, fill = venue)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 0.5, color = "lightgrey", linewidth = 0.25) +
  geom_hline(yintercept = 1.0, color = "lightgrey", linewidth = 0.25) +
  facet_grid(. ~ venue, switch = "y") +
  geom_text(
    data = plottests,
    aes(label = label_mean),
    x = 0.6, y = 0.25, hjust = 0,
    size = 3.5
  ) +
  geom_text(
    data = plottests,
    aes(label = label_tests),
    x = 0.04, y = 1.3, hjust = 0,
    size = 3.5
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  scale_y_continuous(
    limits = c(0, 2.0),
    position = "right",
    expand = c(0, 0),
    breaks = c(0, 1, 2)
  ) +
  scale_fill_manual(values = venue_colors()) +
  theme_criticism() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(hjust = 0),
    axis.title.y = element_blank(),
    panel.spacing.y = unit(0.50, "cm", data = NULL),
    legend.position = "none"
  ) +
  xlab("Percentile rank")

ggsave(
  p,
  filename = snakemake@output[[1]],
  width = 10, height = 2.5,
  bg = "white"
)


#
# Now, lets save the table to a separate file...
#
library(xtable)

stats_table_1s <- plottests %>%
  select(venue, n, mu, ks.1s.statistic, ks.1s.p.value) %>%
  mutate(
    ks.1s.statistic = formatC(
      round(ks.1s.statistic, 4),
      digits = 4,
      format = "f"
    ),
    ks.1s.p.value = formatC(round(ks.1s.p.value, 4), digits = 4, format = "f"),
    mu = formatC(round(mu, 1), digits = 1, format = "f")
  ) %>%
  arrange(venue) %>%
  rename(
    `Venue` = venue,
    `N` = n,
    `Mean Rank` = mu,
    `Test Statistic` = ks.1s.statistic,
    `P Value` = `ks.1s.p.value`
  )

latex_table_1s <- xtable(
  stats_table_1s,
  align = c("llrrrr"),
  digits = 3
)

print(
  latex_table_1s,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[2]]
)
