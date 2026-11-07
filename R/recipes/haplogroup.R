# Adapted from the user's original haplogroup final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_bar_panel <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    fill_value <- if (mode == "color") {
        BAR_COLOR
    }
    else {
        BAR_GRAY
    }
    ggplot(overall, aes(x = Percent, y = Final_Haplogroup)) + geom_col(fill = fill_value, colour = "black", linewidth = 0.23, 
        width = 0.72) + geom_text(aes(label = Label), hjust = -0.06, family = style$font, size = 2.55, colour = "black") + 
        scale_x_continuous(labels = function(x) {
            paste0(x, "%")
        }, expand = expansion(mult = c(0, 0.23))) + scale_y_discrete(limits = rev(FINAL_HG_ORDER), drop = FALSE) + labs(x = "Final Haplogroup frequency (%)", 
        y = NULL) + base_theme + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(t = 5, 
        r = 3, b = 4, l = 4, unit = "mm"))
}
build_heatmap_panel <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    max_pct <- max(region_freq$Percent)
    p <- ggplot(region_freq, aes(x = Region, y = Final_Haplogroup, fill = Percent)) + geom_tile(colour = "white", linewidth = 0.42) + 
        scale_y_discrete(limits = rev(FINAL_HG_ORDER), drop = FALSE) + labs(x = "Geographic region", y = NULL, fill = "Within-region\nfrequency (%)") + 
        base_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(size = 7.9, colour = "black"), panel.border = element_rect(colour = "black", fill = NA, 
            linewidth = 0.32), legend.position = "right", legend.title = element_text(size = 7.6, colour = "black"), legend.text = element_text(size = 7.4, 
            colour = "black"), plot.margin = margin(t = 5, r = 2, b = 4, l = 0, unit = "mm"))
    if (mode == "color") {
        p <- p + scale_fill_gradientn(colours = c(HEAT_LOW, HEAT_MID, HEAT_HIGH), values = scales::rescale(c(0, max_pct * 
            0.35, max_pct)), limits = c(0, max_pct), oob = scales::squish)
    }
    else {
        p <- p + scale_fill_gradient(low = "white", high = "grey15", limits = c(0, max_pct), oob = scales::squish)
    }
    p
}
make_combined_grob <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    bar_grob <- ggplotGrob(build_bar_panel(mode))
    heat_grob <- ggplotGrob(build_heatmap_panel(mode))
    layout <- grid.layout(nrow = 2, ncol = 2, heights = unit(c(0.05, 0.95), "npc"), widths = unit(c(1.2, 1.05), "null"))
    grobTree(rectGrob(gp = gpar(fill = "white", col = NA)), grobTree(bar_grob, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)), 
        grobTree(heat_grob, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)), grobTree(textGrob("(a)", x = unit(0.015, 
            "npc"), y = unit(0.5, "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5, 
            fontface = "plain")), vp = viewport(layout.pos.row = 1, layout.pos.col = 1)), grobTree(textGrob("(b)", x = unit(0.015, 
            "npc"), y = unit(0.5, "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5, 
            fontface = "plain")), vp = viewport(layout.pos.row = 1, layout.pos.col = 2)), vp = viewport(layout = layout))
}
base_theme <- theme_classic(base_family = style$font, base_size = 8.5) + theme(axis.title = element_text(size = 9, colour = "black"), 
    axis.text = element_text(size = 8.2, colour = "black"), axis.line = element_line(linewidth = 0.38, colour = "black"), 
    axis.ticks = element_line(linewidth = 0.32, colour = "black"), plot.title = element_blank(), plot.subtitle = element_blank(), 
    panel.grid = element_blank(), plot.margin = margin(t = 4, r = 5, b = 4, l = 4, unit = "mm"))
