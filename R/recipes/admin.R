# Adapted from the user's original admin final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
make_panel <- function(metric, xlab, mode = c("color", "grayscale"), show_y = TRUE) {
    mode <- match.arg(mode)
    if (mode == "color") {
        fill_scale <- scale_fill_manual(values = REGION_COLORS)
        col_scale <- scale_colour_manual(values = REGION_COLORS)
    }
    else {
        gray_cols <- style$region_grays
        fill_scale <- scale_fill_manual(values = gray_cols)
        col_scale <- scale_colour_manual(values = gray_cols)
    }
    p <- ggplot(dat, aes(x = .data[[metric]], y = Administrative_division, shape = Region, fill = Region, colour = Region)) + 
        geom_segment(aes(x = 0, xend = .data[[metric]], yend = Administrative_division), linewidth = 0.3, colour = "grey80") + 
        geom_point(size = 2.8, stroke = 0.55) + scale_shape_manual(values = REGION_SHAPES) + fill_scale + col_scale + labs(x = xlab, 
        y = if (show_y) 
            "Administrative division"
        else NULL) + theme_pub
    if (!show_y) {
        p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    p
}
build_combined <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    p1 <- make_panel("DC", "Discrimination capacity", mode, TRUE) + labs(tag = "(a)")
    p2 <- make_panel("Fraction_unique_haplotypes", "Fraction of unique haplotypes", mode, FALSE) + labs(tag = "(b)")
    p3 <- make_panel("HMP", "Haplotype match probability", mode, FALSE) + labs(tag = "(c)")
    (p1 | p2 | p3) + plot_layout(guides = "collect", widths = c(1.18, 1, 1)) & theme(legend.position = "top", plot.tag = element_text(family = style$font, 
        size = 9.5, face = "plain"))
}
theme_pub <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), 
    axis.text = element_text(size = 8, colour = "black"), legend.position = "top", legend.title = element_blank(), plot.title = element_blank(), 
    plot.subtitle = element_blank())
