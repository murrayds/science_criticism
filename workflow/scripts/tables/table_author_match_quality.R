suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(xtable)

source("scripts/common.R")
source("scripts/plotting/theme.R")

m <- read_csv(snakemake@input[[1]], col_types = cols()) %>%
  collapse_aps() %>%
  filter(partition == "before") %>%
  mutate(
    treat = factor(treat, levels = c(1, 0), labels = c("Treatment", "Control")),
    venue = factor(venue, levels = venue_levels())
  ) %>%
  arrange(match.group) %>%
  group_by(match.group) %>%
  filter(n() >= 2)

tab <- m %>%
  select(venue, treat, career_age, total_prod, impact_raw, prod_norm) %>%
  pivot_longer(
    names_to = "metric",
    values_to = "value",
    c("career_age", "impact_raw", "prod_norm")
  ) %>%
  group_by(metric, venue, treat) %>%
  summarize(
    mu = round(mean(value), 3)
  ) %>%
  pivot_wider(names_from = "metric", values_from = "mu") %>%
  rename(
    `Career Age` = `career_age`,
    `Impact` = `impact_raw`,
    `Productivity` = `prod_norm`,
    `Journal` = venue,
    `Type` = `treat`
  )

# Construct the table
latex_table <- xtable(
  tab,
  align = c("lllrrr"),
  digits = 3,
)

# Output to a file
print(
  latex_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.colnames.function = function(x) { x },
  file = snakemake@output[[1]]
)