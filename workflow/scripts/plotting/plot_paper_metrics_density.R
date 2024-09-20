library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

source("scripts/common.R")
source("scripts/plotting/theme.R")

# Get the base file...
df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]]) %>%
  mutate(year = as.integer(year)) %>%
  mutate(period = cut(year, breaks = 4)) # add the year/period distinction...

# Load the field file...
fields <- read_csv(snakemake@input[[3]], col_types = cols())

# Load the field categories of referenced papers
refs <- read_csv(snakemake@input[[4]], col_types = cols()) %>%
  rename(id = PaperId)

# Load the field categories of citing paper
cites <- read_csv(snakemake@input[[5]], col_types = cols()) %>%
  rename(id = CitedPaperId)

# Load the novelty of papers...
novelty <- read_csv(snakemake@input[[6]], col_types = cols())

# Load the matched records...
matched <- read_csv(snakemake@input[[7]], col_types = cols())

#
# DATA PREPARATION
#

# Organize by level 0 field...
df_byfield <- fields %>%
  inner_join(df, by = "id") %>%
  group_by(id) %>%
  filter(!duplicated(field_level0)) %>%
  ungroup() %>%
  select(-field) %>%
  rename(field = field_level0)

# Calculate simpson index of referenced papers
ref_diversity <- refs %>%
  left_join(df, by = "id") %>%
  group_by(id) %>%
  # Consdier only papers with at least 10 references...
  filter(length(unique(PaperReferenceId)) >= 10) %>%
  group_by(period, venue, type, id, field) %>%
  summarize(count = n()) %>%
  group_by(id) %>%
  summarize(
    # use reverse simpsons for consistent interpretation
    simpson = -sum((count / sum(count)) ^ 2)
  )

# Calculate simpson index of citing papers
cite_diversity <- cites %>%
  left_join(df, by = "id") %>%
  group_by(id) %>%
  # Consdier only papers with at least 10 citations...
  # Citations are counted within 5 years of the cited paper's publication
  filter(length(unique(CitingPaperId)) > 10) %>% # require at least 10 citations
  group_by(period, venue, type, id, field) %>%
  summarize(count = n()) %>%
  group_by(id) %>%
  summarize(
    # use reverse simpsons for consistent interpretation
    simpson = -sum((count / sum(count)) ^ 2)
  )

#
# REFERENCE DIVERSITY DATA
# 

# Calculate percentile ranks by venue, period, and level0 field
# for reference diversity
plotdata_ref_diversity <- df_byfield %>%
  left_join(ref_diversity, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(simpson)
  ) %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(rank = mean(rank, na.rm = TRUE)) %>%
  mutate(metric = "Reference Diversity")

#
# CITATION DIVERSITY DATA
# 

# Calculate percentile ranks by venue, period, and level0 field
# for citation diversity
plotdata_cite_diversity <- df_byfield %>%
  left_join(cite_diversity, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(simpson),
  ) %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(rank = mean(rank, na.rm = TRUE)) %>%
  mutate(metric = "Citation Diversity")

#
# NOVELTY DATA
#

# Calculate percentile ranks by venue, period, and level0 field
# for each papers' novelty
plotdata_novelty <- df_byfield %>%
  left_join(novelty, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    # use reverse novelty for consistent interpretation
    rank = percent_rank(-Zscore_10th)
  ) %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(
    rank = mean(rank),
  ) %>%
  mutate(metric = "Novelty")

#
# IMPACT DATA
#

# Calculate percentile ranks by venue, period, and level0 field
# for each papers' 3-year impact...
plotdata_impact <- df_byfield %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(impact_3year)
  ) %>%
  group_by(venue, type, id) %>%
  summarize(
    rank = mean(rank, na.rm = TRUE)
  ) %>%
  mutate(metric = "Impact")

# Aggregate all individual plotdata_* objects into a single dataframe
df_all <- data.table::rbindlist(
  list(
    plotdata_impact,
    plotdata_ref_diversity,
    plotdata_cite_diversity,
    plotdata_novelty
  ),
  use.names = TRUE
) %>%
  filter(!is.infinite(rank)) %>% # this happens in a couple of instances
  mutate(
    metric = factor(
      metric,
      levels = c(
        "Impact",
        "Reference Diversity",
        "Citation Diversity",
        "Novelty"
      ),
      labels = c(
        "Impact",
        "Reference\nDiversity",
        "Citation\nDiversity",
        "Novelty"
      ),
    )
  )

#
# MATCHING DATA
# 

# Identify records correponding to matched papers, with a similar 
# field, year, and citation impact to the papers targeted by criticism
plotdata_matched <- df_all %>%
  filter(type == "article") %>%
  inner_join(matched %>% select(id), by = "id") %>%
  filter(metric != "Impact")

# Now begin building a table holding the results of one-sample KS tests
# comparing the distribution of observed ranks of criticism-targeted papers
# against a uniform distribution
one_sample_stats <- df_all %>%
  # we only care about letters here...
  filter(type == "letter") %>%
  filter(!is.na(rank)) %>%
  # KS test complains if there are exact ties. Because ties are innevitable
  # in the data, we add a tiny jutter to all values so they are not exactly
  # the same.
  mutate(rank = jitter(rank, factor = 1e-8)) %>%
  group_by(venue, metric) %>%
  do(
    # KS test, copmared to uniform distribution between 0 and 1
    ks = ks.test(.$rank, "punif", min = 0, max = 1, alternative = "less"),
    mu = mean(.$rank) * 100,
    n = length(unique(.$id))
  ) %>%
  summarize(
    ks.1s.statistic = ks$statistic,
    ks.1s.p.value = ks$p.value,
    mu = mu,
    n = n,
    venue = venue, metric = metric
  ) %>%
  rowwise() %>%
  mutate(
    ks.1s.statistic = round(ks.1s.statistic, 3),
    ks.1s.p.value = round(ks.1s.p.value, 3)
  )

# Now construct a second table, this one containing the results of two-sample
# KS tests comparing the distribution of percentile ranks of the Letters
# against that of a marched population of non-letters with similar impact...
two_sample_stats <- df_all %>%
  # As in the 1-sample test, we need a small jitter value
  mutate(rank = jitter(rank, factor = 1e-8)) %>%
  inner_join(matched %>% select(id, match.group), by = "id") %>%
  filter(!is.na(rank)) %>%
  # Limit to only those records for which we could find a match...
  group_by(venue, metric, match.group) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(venue, metric, match.group),
    names_from = "type",
    values_from = "rank"
  ) %>%
  arrange(venue, match.group) %>%
  group_by(venue, metric) %>%
  do(
    # conduct the one-sided two-sample KS test
    ks = ks.test(.$letter, .$article, alternative = "less"),
    mu_letter = mean(.$letter) * 100,
    mu_article = mean(.$article) * 100,
    n = length(.$letter)
  ) %>%
  summarize(
    ks.2s.statistic = ks$statistic,
    ks.2s.p.value = round(ks$p.value, 3),
    mu_letter = mu_letter, mu_article = mu_article,
    venue = venue, metric = metric
  )

# Now, combine the two table into one and form the labels...
plotlabs <- one_sample_stats %>%
  left_join(two_sample_stats, by = c("venue", "metric")) %>%
  rowwise() %>%
  mutate(
    label_mean = paste0(
      "μ = ",
      formatC(round(last(mu), 1), digits = 1, format = "f"),
      "%"
    ),
    label_tests = paste0(
      "\n",
      "1s KS, p = ",
      formatC(round(last(ks.1s.p.value), 2), digits = 2, format = "f"),
      "\n"
    ),
    label_tests = ifelse(
      metric == "Impact",
      label_tests,
      paste0(
        label_tests,
        "2s KS, p = ",
        formatC(round(last(ks.2s.p.value), 2), digits = 2, format = "f")
      )
    )
  )

#
# Construct the plot object...
#
p <- df_all %>%
  filter(type == "letter") %>%
  filter(!is.na(rank)) %>%
  ggplot(aes(x = rank,  fill = venue)) +
  # Add a vertical line to mark the median for easier viewing
  geom_vline(xintercept = 0.5, color = "lightgrey", linewidth = 0.25) +
  # The main density
  geom_density(alpha = 0.2) +
  # A density for the matched records...
  geom_density(
    data = plotdata_matched,
    fill = NA,
    linetype = "dashed",
    color = "darkslategrey"
  ) +
  geom_text(
    data = plotlabs,
    aes(label = label_mean),
    x = 0.6, y = 0.25, hjust = 0,
    size = 3.5
  ) +
  geom_text(
    data = plotlabs,
    aes(label = label_tests),
    x = 0.04, y = 1.95, hjust = 0,
    size = 3.5
  ) +
  facet_grid(metric ~ venue, switch = "y") +
  scale_fill_manual(values = venue_colors(), guide = "none") +
  scale_x_continuous(
    expand = c(0, 0),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  scale_y_continuous(
    limits = c(0, 2.2),
    position = "right",
    expand = c(0, 0)
  ) +
  theme_criticism() +
  theme(
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0, hjust = 0),
    axis.title.y = element_blank(),
    panel.spacing.x = unit(0.30, "cm", data = NULL)
  ) +
  xlab("Percentile rank")


ggsave(p, filename = snakemake@output[[1]], width = 9, height = 6, bg = "white")

#
# Now, lets save the table to a separate file...
#
library(xtable)

stats_table_1s <- plotlabs %>%
  select(metric, venue, n, mu, ks.1s.statistic, ks.1s.p.value) %>%
  mutate(
    ks.1s.statistic = formatC(
      round(ks.1s.statistic, 4),
      digits = 4,
      format = "f"
    ),
    ks.1s.p.value = formatC(round(ks.1s.p.value, 4), digits = 4, format = "f"),
    mu = formatC(round(mu, 1), digits = 1, format = "f")
  ) %>%
  arrange(metric, venue) %>%
  rename(
    `Metric` = metric,
    `Venue` = venue,
    `N` = n,
    `Mean Rank` = mu,
    `Test Statistic` = ks.1s.statistic,
    `P Value` = `ks.1s.p.value`
  )

latex_table_1s <- xtable(
  stats_table_1s,
  align = c("lllrrrr"),
  digits = 3
)

print(
  latex_table_1s,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[2]]
)

stats_table_2s <- plotlabs %>%
  select(
    metric, venue, n,
    mu_letter, mu_article,
    ks.2s.statistic, ks.2s.p.value
  ) %>%
  mutate(
    ks.2s.statistic = formatC(
      round(ks.2s.statistic, 4),
      digits = 4,
      format = "f"
    ),
    ks.2s.p.value = formatC(round(ks.2s.p.value, 4), digits = 4, format = "f"),
    mu_letter = formatC(round(mu_letter, 1), digits = 1, format = "f"),
    mu_article = formatC(round(mu_article, 1), digits = 1, format = "f")
  ) %>%
  arrange(metric, venue) %>%
  rename(
    `Metric` = metric,
    `Venue` = venue,
    `N` = n,
    `Mean Rank (Criticism)` = mu_letter,
    `Mean Rank (¬Criticism)` = mu_article,
    `Test Statistic` = ks.2s.statistic,
    `P Value` = `ks.2s.p.value`
  )


latex_table_2s <- xtable(
  stats_table_2s,
  align = c("lllrrrrr"),
  digits = 3
)

print(
  latex_table_2s,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[3]]
)