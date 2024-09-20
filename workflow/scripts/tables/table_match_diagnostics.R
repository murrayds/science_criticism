suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(xtable)

source("scripts/common.R")
source("scripts/plotting/theme.R")

agg <- read_csv(snakemake@input[[1]], col_types = cols())


mutate_counts_table <- function(df) {
  df %>%
    select(venue, count_candidates, num_matched, delay, impact, year) %>%
    rowwise() %>%
    mutate(
      value = paste0(
        count_candidates,
        " (",
        round(num_matched / count_candidates * 100, 1),
        "%)"
      )
    ) %>%
    select(-count_candidates, -num_matched)
}

mutate_tstats_table <- function(df) {
  print(names(df))
  df %>%
    select(venue, t.estimate, t.p.value, delay, impact, year) %>%
    rowwise() %>%
    mutate(
      value = paste0(
        round(t.estimate, 2),
        " (p = ",
        formatC(round(t.p.value, 3), format = "f", digits = 3),
        ")"
      )
    ) %>%
    select(-t.estimate, -t.p.value)
}

mutate_wilcox_table <- function(df) {
  print(names(df))
  df %>%
    select(venue, wilcox.p.value, delay, impact, year) %>%
    rowwise() %>%
    mutate(
      value = paste0(
        "p = ",
        formatC(round(wilcox.p.value, 3), format = "f", digits = 3)
      )
    ) %>%
    select(-wilcox.p.value)
}

# First, lets perform some formatting on the table
agg_formatted <- agg %>%
  mutate(
    venue = factor(venue, levels = venue_levels())
  )

param <- snakemake@wildcards[[1]]
if (param == "counts") {
  agg_formatted <- agg_formatted %>% mutate_counts_table()
} else if (param == "tstats") {
  agg_formatted <- agg_formatted %>% mutate_tstats_table()
} else if (param == "wilcox") {
  agg_formatted <- agg_formatted %>% mutate_wilcox_table()
}

agg_formatted <- agg_formatted %>%
  pivot_wider(names_from = venue, values_from = value) %>%
  select(delay, impact, year, venue_levels())



# We will narrow the table into three separate pieces, and
# for each we will vary only a single parameter. 
delay_table <- agg_formatted %>%
  filter(impact == 0.1, year == 1)
delay_table[nrow(delay_table) + 1, ] <- NA

cite_tolerance_table <- agg_formatted %>%
  filter(delay == 3, year == 1)
cite_tolerance_table[nrow(cite_tolerance_table) + 1, ] <- NA

year_tolerance_table <- agg_formatted %>%
  filter(delay == 3, impact == 0.1)

# Aggregate mini tables, perform final polish
tab <- data.table::rbindlist(
  list(delay_table, cite_tolerance_table, year_tolerance_table)
) %>%
  mutate(impact = ifelse(is.na(impact), NA, paste0(impact * 100, "%"))) %>%
  rename(
    `Delay` = delay,
    `Impact $\\pm$ $\\epsilon$` = impact,
    `Year $\\pm$ $\\epsilon$` = year
  )


# Construct the table
latex_table <- xtable(
  tab,
  align = c("lccccccccc"),
  digits = 0,
)

# Output to a file
print(
  latex_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.colnames.function = function(x) { x },
  file = snakemake@output[[1]]
)
