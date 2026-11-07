# Adapted from the user's original rarefaction final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_admin_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    palette <- if (mode == "color") {
        REGION_COLORS
    }
    else {
        REGION_GRAYS
    }
    ggplot(admin_df, aes(x = Group, y = Mean, colour = Region, shape = Region)) + geom_errorbar(aes(ymin = Lower, ymax = Upper), 
        width = 0.16, linewidth = 0.58) + geom_point(aes(fill = Region), size = 2.8, stroke = 0.55) + scale_colour_manual(values = palette, 
        drop = FALSE) + scale_fill_manual(values = palette, drop = FALSE) + scale_shape_manual(values = REGION_SHAPES, drop = FALSE) + 
        scale_x_discrete(limits = ADMIN_ORDER, drop = FALSE) + labs(x = "Administrative division", y = paste0("Rarefied number of distinct haplotypes\n(subsample size = ", 
        ADMIN_SUBSAMPLE_N, ")")) + guides(fill = "none", colour = guide_legend(override.aes = list(shape = unname(REGION_SHAPES[REGION_ORDER]), 
        size = 3)), shape = "none") + publication_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
}
build_region_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    palette <- if (mode == "color") {
        REGION_COLORS
    }
    else {
        REGION_GRAYS
    }
    ggplot(region_df, aes(x = Group, y = Mean, colour = Group, shape = Group)) + geom_errorbar(aes(ymin = Lower, ymax = Upper), 
        width = 0.16, linewidth = 0.62) + geom_point(aes(fill = Group), size = 3.4, stroke = 0.6) + scale_colour_manual(values = palette, 
        drop = FALSE) + scale_fill_manual(values = palette, drop = FALSE) + scale_shape_manual(values = REGION_SHAPES, drop = FALSE) + 
        scale_x_discrete(limits = REGION_ORDER, drop = FALSE) + labs(x = "Geographic region", y = paste0("Rarefied number of distinct haplotypes\n(subsample size = ", 
        REGION_SUBSAMPLE_N, ")")) + guides(colour = "none", fill = "none", shape = "none") + publication_theme + theme(legend.position = "none")
}
make_combined_grob <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    admin_grob <- ggplotGrob(build_admin_plot(mode))
    region_grob <- ggplotGrob(build_region_plot(mode))
    layout <- grid.layout(nrow = 4, ncol = 1, heights = unit(c(0.035, 0.5, 0.035, 0.43), "npc"))
    grobTree(rectGrob(gp = gpar(fill = "white", col = NA)), grobTree(textGrob("(a)", x = unit(0.012, "npc"), y = unit(0.5, 
        "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5, fontface = "plain")), vp = viewport(layout.pos.row = 1, 
        layout.pos.col = 1)), grobTree(admin_grob, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)), grobTree(textGrob("(b)", 
        x = unit(0.012, "npc"), y = unit(0.5, "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5, 
            fontface = "plain")), vp = viewport(layout.pos.row = 3, layout.pos.col = 1)), grobTree(region_grob, vp = viewport(layout.pos.row = 4, 
        layout.pos.col = 1)), vp = viewport(layout = layout))
}
publication_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, 
    colour = "black"), axis.text = element_text(size = 8.3, colour = "black"), axis.line = element_line(linewidth = 0.42, 
    colour = "black"), axis.ticks = element_line(linewidth = 0.34, colour = "black"), panel.grid = element_blank(), plot.title = element_blank(), 
    plot.subtitle = element_blank(), plot.margin = margin(t = 4, r = 5, b = 5, l = 5, unit = "mm"), legend.position = "top", 
    legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 8))
