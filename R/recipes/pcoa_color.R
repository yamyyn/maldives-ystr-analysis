# Adapted from the user's original pcoa_color final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_plot <- function() {
    ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey65") + 
        geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey65") + geom_point(aes(shape = Region, 
        fill = Region, colour = Region), size = 3.5, stroke = 0.85) + ggrepel::geom_text_repel(aes(label = Admin), family = style$font, 
        size = 3.1, colour = "black", box.padding = 0.2, point.padding = 0.12, min.segment.length = 0.1, segment.colour = "grey70", 
        segment.size = 0.22, max.overlaps = Inf, seed = 20261010) + scale_fill_manual(values = REGION_COLOURS, drop = FALSE) + 
        scale_colour_manual(values = REGION_COLOURS, drop = FALSE) + scale_shape_manual(values = REGION_SHAPES, drop = FALSE) + 
        labs(x = paste0("PCoA 1 (", sprintf("%.2f", LOCKED_AXIS1_PERCENT), "%)"), y = paste0("PCoA 2 (", sprintf("%.2f", 
            LOCKED_AXIS2_PERCENT), "%)"), fill = "Region", colour = "Region", shape = "Region") + coord_cartesian(clip = "off") + 
        guides(fill = "none", colour = "none", shape = guide_legend(title.position = "left", title.hjust = 0.5, override.aes = list(size = 3.8, 
            fill = unname(REGION_COLOURS[REGION_LEVELS]), colour = unname(REGION_COLOURS[REGION_LEVELS])))) + theme_classic(base_family = style$font, 
        base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), axis.text = element_text(size = 8.5, 
        colour = "black"), axis.line = element_line(linewidth = 0.4, colour = "black"), axis.ticks = element_line(linewidth = 0.35, 
        colour = "black"), legend.position = "top", legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_text(size = 8.5, 
        colour = "black"), legend.text = element_text(size = 8, colour = "black"), legend.key.width = grid::unit(0.7, "cm"), 
        legend.spacing.x = grid::unit(0.08, "cm"), plot.title = element_blank(), plot.subtitle = element_blank(), plot.margin = margin(t = 4, 
            r = 8, b = 5, l = 5, unit = "mm"))
}
