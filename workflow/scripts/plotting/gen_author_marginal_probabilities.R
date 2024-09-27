suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(margins))
suppressPackageStartupMessages(library(sandwich))
suppressPackageStartupMessages(library(lmtest))

authors <- read_csv(snakemake@input[[1]], col_types = cols())
histories <- read_csv(snakemake@input[[2]], col_types = cols())
matched <- read_csv(snakemake@input[[3]], col_types = cols())
elite <- read_csv(snakemake@params[[1]], col_types = cols()) %>%
  select(grid_id) %>%
  distinct(grid_id) %>%
  mutate(elite = TRUE)

# Get the authorhsip parameter
wc_authorship <- substring(snakemake@wildcards[[1]], 1, 1)


# From the author histories, identify the first year in which they published
author_first_year <- histories %>%
  filter(author_position == wc_authorship) %>%
  select(AuthorId, Year) %>%
  inner_join(authors %>% select(AuthorId), by = "AuthorId") %>%
  group_by(AuthorId) %>%
  summarize(
    first_year = min(Year),
    .groups = "drop"
  )

# Now build a dataframe of the author information for each paper...
author_df <- authors %>%
  filter(num_authors <= 20) %>% # remove papers with too many authors...
  filter(author_position == wc_authorship) %>%
  left_join(author_first_year, by = "AuthorId") %>%
  left_join(elite, by = c("GridId" = "grid_id")) %>%
  filter(!is.na(Gender)) %>%
  group_by(PaperId) %>%
  summarize(
    authorship = "first",
    elite = any(!is.na(elite)),
    gender = ifelse(Gender == 1, "Female", "Male"),
    first_year = first(first_year),
    num_authors = first(num_authors),
    .groups = "drop"
  ) %>%
  rename(id = "PaperId")

# Now move towards the targeted dataframe...
targeted_df <- matched %>%
  inner_join(author_df, by = c("id")) %>%
  filter(!is.na(gender)) %>%
  mutate(
    gender = factor(gender, levels = c("Male", "Female")),
    elite = ifelse(!elite, "NonElite", "Elite"),
    elite = factor(elite, levels = c("NonElite", "Elite")),
    impact = impact_norm,
    type = factor(type, levels = c("Control", "Treatment")),
    career_age = year - first_year,
    seniority = ifelse(career_age <= 10, "Junior", "Senior"),
    seniority = factor(seniority, levels = c("Junior", "Senior")),
    year = year - min(year),
    num_authors = cut(num_authors, c(0, 1, 5, 20), include_lowest = TRUE),
    venue_group = ifelse(
      venue %in% c("PR-A", "PR-B", "PR-C", "PR-D", "PR-E"),
      "Other APS",
      venue
    )
  )

venues <- c("Nature", "Science", "PNAS", "PRL", "Other APS")
margins_dfs <- lapply(venues, function(v) {
  model <- glm(
    type ~
      gender + seniority + elite + num_authors + impact + year,
    data = targeted_df %>% filter(venue_group == v),
    family = binomial
  )

  num_observations <- nobs(model)

  marginal_effects_gender <- margins(
    model,
    variable = "gender"
  )

  marginal_effects_seniority <- margins(
    model,
    variable = "seniority"
  )

  marginal_effects_elite <- margins(
    model,
    variable = "elite"
  )

  # coef_se_clustered <- coeftest(
  #   model,
  #   vcov = vcovCL(model, cluster = ~match.group)
  # )

  # Extract the coefficients and standard errors from the coeftest output
  # coef_estimates <- coef_se_clustered[, 1]  # Coefficients
  # coef_se <- coef_se_clustered[, 2]         # Clustered standard errors

  # se <- coef_se[which(names(coef_estimates) == "genderFemale")]
  margins_df_gender <- as.data.frame(summary(marginal_effects_gender)) %>%
    #select(factor, AME) %>%
    mutate(
      venue = v,
      #lower_ci = AME - 1.96 * se,
      #upper_ci = AME + 1.96 * se
    )

  # se <- coef_se[which(names(coef_estimates) == "senioritySenior")]
  margins_df_seniority <- as.data.frame(summary(marginal_effects_seniority)) %>%
    #select(factor, AME) %>%
    mutate(
      venue = v,
      #lower_ci = AME - 1.96 * se,
     # upper_ci = AME + 1.96 * se
    )

  # se <- coef_se[which(names(coef_estimates) == "eliteElite")]
  margins_df_elite <- as.data.frame(summary(marginal_effects_elite)) %>%
    #select(factor, AME) %>%
    mutate(
      venue = v,
      #lower_ci = AME - 1.96 * se,
      #upper_ci = AME + 1.96 * se
    )

  df_part <- rbind(
    margins_df_gender,
    margins_df_seniority,
    margins_df_elite
  )

  df_part$n <- num_observations
  return(df_part)
})

# Combine all the results into one dataframe...
margins_df <- data.table::rbindlist(margins_dfs) %>%
  arrange(factor, venue)

write_csv(margins_df, snakemake@output[[1]])
