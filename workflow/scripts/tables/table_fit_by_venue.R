suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
library(xtable)

source("scripts/common.R")

df <- load_aggregate_df(snakemake@input[[1]], snakemake@input[[2]]) %>%
  collapse_aps()
  
fields <- read_csv(snakemake@input[[3]], col_types = cols())

# Generate the raw and field-normalized impact for these data
df_ranked <- fields %>%
  filter(level == 1) %>%
  inner_join(df, by = "id") %>%
  mutate(year = as.integer(year)) %>%
  mutate(period = cut(year, breaks = 6)) %>%
  group_by(venue, period, field) %>%
  mutate(
    impact_denominator = mean(impact_2year)
  ) %>%
  group_by(id, venue, period, year, type) %>%
  dplyr::summarize(
    raw_impact = mean(impact_2year),
    norm_impact = mean(impact_2year / impact_denominator)
  ) %>%
  filter(raw_impact > 0)

# Perform the regressions on each venue and for each metric
t <- df_ranked %>%
  pivot_longer(
    cols = c(raw_impact, norm_impact),
    names_to = "impact_type",
    values_to = "impact"
  ) %>%
  group_by(venue, impact_type) %>%
  do(
    mod = rms::lrm(
      type ~ period + log(impact),
      data = .,
    )$stats[["R2"]]
  ) %>%
  mutate(mod = unlist(mod)) %>%
  mutate(mod = round(mod, 4)) %>%
  pivot_wider(names_from = "impact_type", values_from = "mod") %>%
  select(venue, raw_impact, norm_impact) %>%
  rename(
    `Venue` = venue,
    `2-Year Impact` = raw_impact,
    `Normalized 2-Year Impact` = norm_impact
  )

# Construct the latex table
latex_table <- xtable(
  t,
  align = c("llrr"),
  digits = 3
)

print(
  latex_table,
  include.rownames = FALSE,
  booktabs = TRUE,
  file = snakemake@output[[1]]
)