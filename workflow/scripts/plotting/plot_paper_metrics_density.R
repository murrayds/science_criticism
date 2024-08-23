library(readr)
suppressPackageStartupMessages(library(dplyr))
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

# Calculate percentile ranks by venue, period, and level0 field
# for reference diversity
plotdata_ref_diversity <- df_byfield %>%
  left_join(ref_diversity, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(simpson)
  ) %>%
  filter(type == "letter") %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(rank = mean(rank, na.rm = TRUE)) %>%
  mutate(metric = "Reference Diversity")

# Calculate percentile ranks by venue, period, and level0 field
# for citation diversity
plotdata_cite_diversity <- df_byfield %>%
  left_join(cite_diversity, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(simpson),
  ) %>%
  filter(type == "letter") %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(rank = mean(rank, na.rm = TRUE)) %>%
  mutate(metric = "Citation Diversity")

# Calculate percentile ranks by venue, period, and level0 field
# for each papers' novelty
plotdata_novelty <- df_byfield %>%
  left_join(novelty, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    # use reverse novelty for consistent interpretation
    rank = percent_rank(-Atyp_10pct_Z)
  ) %>%
  filter(type == "letter") %>%
  # When multiple fields are present, use mean rank...
  group_by(venue, type, id) %>%
  summarize(
    rank = mean(rank),
  ) %>%
  mutate(metric = "Novelty")

# Calculate percentile ranks by venue, period, and level0 field
# for each papers' 3-year impact...
plotdata_impact <- df_byfield %>%
  left_join(novelty, by = "id") %>%
  group_by(venue, period, field) %>%
  mutate(
    rank = percent_rank(impact_3year)
  ) %>%
  filter(type == "letter") %>%
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

# Construct a table containing labels of the mean rank for each venue/metric,
# as well as the p-value result of a KS test of whether the percentile_ranks of
# the letters are drawn uniformly
plotlabs <- df_all %>%
  filter(!is.na(rank)) %>%
  # KS test complains if there are exact ties. Because ties are innevitable
  # in the data, we add a tiny jutter to all values so they are not exactly
  # the same.
  mutate(rank = jitter(rank, factor = 1e-8)) %>%
  group_by(venue, metric) %>%
  do(
    # KS test, copmared to uniform distribution between 0 and 1
    ks = ks.test(.$rank, "punif", min = 0, max = 1),
    mu = mean(.$rank) * 100,
    n = length(unique(.$id))
  ) %>%
  summarize(
    ks.statistic = ks$statistic,
    ks.p.value = ks$p.value,
    mu = mu,
    n = n,
    venue = venue, metric = metric
  ) %>%
  rowwise() %>%
  mutate(
    ks.statistic = round(ks.statistic, 1),
    p.value = round(ks.p.value, 3),
    # Add a measure of significance with asterisks...
    signif = case_when(
      ks.p.value > 0.05 ~ "",
      ks.p.value > 0.01 ~ "*",
      ks.p.value > 0.001 ~ "**",
      !is.na(ks.p.value) ~ "***",
      TRUE ~ NA_character_
    ),
    # Construct the label that will appear on the plot...
    label = paste0(
      "μ = ",
      formatC(round(last(mu), 1), digits = 1, format = "f"),
      "%",
      "\n",
      signif
    )
  )

# Construct the plot object...
p <- df_all %>%
  ggplot(aes(x = rank,  fill = venue)) +
  geom_density(alpha = 0.2) +
  geom_text(
    data = plotlabs,
    aes(label = label),
    x = 0.05, y = 1.55, hjust = 0,
  ) +
  facet_grid(metric ~ venue, switch = "y") +
  scale_fill_manual(values = venue_colors(), guide = FALSE) +
  scale_x_continuous(
    expand = c(0, 0),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  scale_y_continuous(limits = c(0, 2.0), position = "right") +
  theme_criticism() +
  theme(
    panel.grid = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    axis.title.y = element_blank(),
    panel.spacing.x = unit(0.30, "cm", data = NULL)
  ) +
  xlab("Percentile rank")


ggsave(p, filename = snakemake@output[[1]], width = 8, height = 5, bg = "white")

#
# Now, lets save the table to a separate file...
#
library(xtable)

stats_table <- plotlabs %>%
  select(metric, venue, n, mu, ks.statistic, ks.p.value) %>%
  mutate(
    ks.statistic = formatC(round(ks.statistic, 4), digits = 4, format = "f"),
    ks.p.value = formatC(round(ks.p.value, 4), digits = 4, format = "f"),
    mu = formatC(round(mu, 1), digits = 1, format = "f")
  ) %>%
  arrange(metric, venue) %>%
  rename(
    `Metric` = metric,
    `Venue` = venue,
    `N` = n,
    `Mean Rank` = mu,
    `Test Statistic` = ks.statistic,
    `P Value` = `ks.p.value`
  )

latex_table <- xtable(
  stats_table,
  align = c("lllrrrr"),
  digits = 3
)

print(
  latex_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[2]]
)