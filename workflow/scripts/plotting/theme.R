# The order of labels to use
venue_levels <- function() {
  c(
    "Nature",
    "Science",
    "PNAS",
    "PRL",
    "Other APS"
  )
}

venue_levels_all <- function() {
  c(
    "Nature",
    "Science",
    "PNAS",
    "PRL",
    "PR-A",
    "PR-B",
    "PR-C",
    "PR-D",
    "PR-E"
  )
}

field_values <- function() {
  c(
    `142362112` = "art",
    `185592680` = "chemistry",
    `71924100` = "medicine",
    `192562407` = "materials science",
    `144024400` = "sociology",
    `17744445` = "political science",
    `39432304` = "environmental science",
    `127413603` = "engineering",
    `41008148` = "computer science",
    `205649164` = "geography",
    `15744967` = "psychology",
    `127313418` = "geology",
    `144133560` = "business",
    `162324750` = "economics",
    `138885662` = "philosophy",
    `33923547` = "mathematics",
    `86803240` = "biology",
    `121332964` = "physics",
    `95457728` = "history"
  )
}

venue_colors <- function() {
  c(
    "Nature" = "forestgreen",
    "Science" = "firebrick",
    "PNAS" = "goldenrod",
    "PRL" = "steelblue",
    "Other APS" = "lightblue"
  )
}


theme_criticism <- function() {
  require(ggplot2)
  theme_minimal() +
  theme(
    panel.background = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.5
    ),
    text = element_text(family = "Helvetica", size = 12),
    legend.background = element_rect(fill = "white", linewidth = 0.6),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    legend.position = "inside",
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )
}