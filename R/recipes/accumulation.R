# Adapted from the user's original accumulation final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_plot <- function() {
    ggplot(plot_data, aes(x = Sample_Size, y = Mean, group = Series)) + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Series), 
        alpha = 0.22, colour = NA) + geom_line(aes(linetype = Series), colour = "black", linewidth = 0.8) + scale_linetype_manual(values = line_types, 
        name = NULL) + scale_fill_manual(values = ribbon_fills, name = NULL) + scale_x_continuous(breaks = style$accumulation_x$breaks, 
        limits = style$accumulation_x$limits, expand = expansion(mult = c(0, 0))) + scale_y_continuous(breaks = style$accumulation_y$breaks, 
        limits = style$accumulation_y$limits, expand = expansion(mult = c(0, 0))) + labs(x = "Sample size", y = "Number of haplotypes") + 
        guides(fill = "none", linetype = guide_legend(override.aes = list(linewidth = 0.9))) + theme_classic(base_family = style$font, 
        base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), axis.text = element_text(size = 8.5, 
        colour = "black"), axis.line = element_line(linewidth = 0.35, colour = "black"), axis.ticks = element_line(linewidth = 0.3, 
        colour = "black"), panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25), panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = c(0.22, 0.88), legend.background = element_rect(fill = "white", 
            colour = "grey78", linewidth = 0.2), legend.key.width = grid::unit(0.75, "cm"), legend.key.height = grid::unit(0.4, 
            "cm"), legend.margin = margin(t = 2, r = 3, b = 2, l = 3, unit = "mm"), legend.text = element_text(size = 7.5), 
        plot.margin = margin(t = 5, r = 6, b = 5, l = 6, unit = "mm"))
}
