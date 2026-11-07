# Adapted from the user's original spatial final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_ibd_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    if (mode == "color") {
        point_fill <- "#9ECAE1"
        point_edge <- "#3F3F3F"
        fit_color <- "#2C7FB8"
        zero_color <- "grey60"
    }
    else {
        point_fill <- "grey80"
        point_edge <- "grey25"
        fit_color <- "black"
        zero_color <- "grey60"
    }
    ggplot(ibd_df, aes(x = geographic_km, y = rst)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, colour = zero_color) + 
        geom_point(shape = 21, size = 2.05, stroke = 0.52, fill = point_fill, colour = point_edge, alpha = 0.68) + geom_smooth(method = "lm", 
        formula = y ~ x, se = FALSE, linewidth = 0.78, colour = fit_color) + annotate("text", x = Inf, y = Inf, label = annotation_label, 
        hjust = 1.05, vjust = 1.2, family = style$font, size = 3.25, colour = "black", lineheight = 1.12) + labs(x = "Geographic distance (km)", 
        y = expression("Pairwise " * italic(R)[ST])) + scale_x_continuous(breaks = style$ibd_x$breaks, limits = style$ibd_x$limits, 
        expand = expansion(mult = c(0, 0))) + scale_y_continuous(breaks = style$ibd_y$breaks, limits = style$ibd_y$limits, 
        labels = scales::label_number(accuracy = 0.001), expand = expansion(mult = c(0, 0))) + coord_cartesian(clip = "off") + 
        base_theme
}
build_corr_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    if (mode == "color") {
        point_fill <- "#9ECAE1"
        point_edge <- "#3F3F3F"
        line_color <- "#2C7FB8"
        zero_color <- "grey60"
    }
    else {
        point_fill <- "grey80"
        point_edge <- "grey25"
        line_color <- "black"
        zero_color <- "grey60"
    }
    ggplot(corr_df, aes(x = distance_class, y = mantel_r, group = 1)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, 
        colour = zero_color) + geom_line(linewidth = 0.78, colour = line_color) + geom_point(shape = 21, size = 3, stroke = 0.62, 
        fill = point_fill, colour = point_edge) + labs(x = "Geographic distance class (km)", y = expression("Mantel " * italic(r))) + 
        scale_y_continuous(breaks = style$correlogram_y$breaks, limits = style$correlogram_y$limits, expand = expansion(mult = c(0, 
            0))) + base_theme + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}
make_combined_grob <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    ibd_grob <- ggplotGrob(build_ibd_plot(mode))
    corr_grob <- ggplotGrob(build_corr_plot(mode))
    lay <- grid.layout(nrow = 4, ncol = 1, heights = unit(c(0.035, 0.47, 0.035, 0.46), "npc"))
    grobTree(rectGrob(gp = gpar(fill = "white", col = NA)), grobTree(textGrob("(a)", x = unit(0.01, "npc"), y = unit(0.5, 
        "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5)), vp = viewport(layout.pos.row = 1, 
        layout.pos.col = 1)), grobTree(ibd_grob, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)), grobTree(textGrob("(b)", 
        x = unit(0.01, "npc"), y = unit(0.5, "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5)), 
        vp = viewport(layout.pos.row = 3, layout.pos.col = 1)), grobTree(corr_grob, vp = viewport(layout.pos.row = 4, layout.pos.col = 1)), 
        vp = viewport(layout = lay))
}
base_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 10, colour = "black"), 
    axis.text = element_text(size = 8.5, colour = "black"), axis.line = element_line(colour = "black", linewidth = 0.45), 
    axis.ticks = element_line(colour = "black", linewidth = 0.4), plot.title = element_blank(), plot.subtitle = element_blank(), 
    legend.position = "none", plot.margin = margin(3, 5, 3, 5, unit = "mm"))
