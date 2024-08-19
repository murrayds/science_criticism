# The order of labels to use
venue_levels <- function() {
  c(
    "Nature",
    "Science",
    "PNAS",
    "PRL"
  )
}

venue_colors <- function() {
  c(
      "Nature" = "forestgreen",
      "Science" = "firebrick",
      "PNAS" = "goldenrod",
      "PRL" = "steelblue"
  ) 
}


theme_criticism <- function() {
  require(ggplot2)
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = "black", linewidth = 0.5),
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