plot_boxplot <- function(data, metric, legend_nrow = 2) {
  #metric <- rlang::as_name(rlang::ensym(metric))
  
  ggplot(
    data,
    aes(
      x = x_group,
      y = .data[[metric]],
      fill = method_code
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.9,
      width = 0.7,
      position = position_dodge2(width = 0.75, preserve = "single")
    ) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # panel.grid = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(fill = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 3)
    ) +
    facet_grid(
      model ~ delta_label,
      labeller = label_parsed,
      scales = "free_y"
    )
}
