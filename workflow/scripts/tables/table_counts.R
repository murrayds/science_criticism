suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(readr)
library(xtable)

source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]])

tab <- df %>%
  filter(!is.na(venue)) %>%
  group_by(venue, type) %>%
  summarize(
    count = n()
  ) %>%
  pivot_wider(names_from = type, values_from = count) %>%
  rename(
    `Venue` = venue,
    `Articles` = `article`,
    `Critical letters` = `letter`
  )

# Construct the table
latex_table <- xtable(
  tab,
  align = "llrr",
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