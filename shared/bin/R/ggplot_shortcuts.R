gg_smaller_legend <- function(p, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    guides(
        shape = guide_legend(override.aes = list(size = pointSize)), 
        color = guide_legend(override.aes = list(size = pointSize))
    ) +
    theme(
        legend.title = element_text(size = textSize), 
        legend.text  = element_text(size = textSize),
        legend.key.size = unit(spaceLegend, "lines")
    )
}

gg_legend_bottom <- function(p) {
    theme(legend.position = 'bottom') + 
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))
}

