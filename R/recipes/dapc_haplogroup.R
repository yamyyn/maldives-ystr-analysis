# Adapted from the user's original dapc_haplogroup final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_hg_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    p <- ggplot(hg_df, aes(x = DA1, y = DA2)) + geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.3) + geom_vline(xintercept = 0, 
        colour = "grey82", linewidth = 0.3)
    if (mode == "color") {
        p <- p + geom_point(aes(colour = Group), shape = 16, size = 1.55, alpha = 0.75, show.legend = FALSE) + geom_point(data = hg_centroids, 
            aes(colour = Group), shape = 21, fill = "white", size = 3.6, stroke = 0.9, show.legend = FALSE) + scale_colour_manual(values = HG_COLOR_PALETTE, 
            drop = FALSE)
    }
    else {
        p <- p + geom_point(colour = "grey80", shape = 1, size = 1.1, alpha = 0.7, stroke = 0.35) + geom_point(data = hg_centroids, 
            colour = "black", shape = 21, fill = "white", size = 3.8, stroke = 0.95)
    }
    p + ggrepel::geom_label_repel(data = hg_centroids, aes(label = Group), colour = "black", fill = "white", family = style$font, 
        size = 2.85, label.size = 0.18, box.padding = 0.22, point.padding = 0.15, min.segment.length = 0, segment.colour = "grey65", 
        segment.size = 0.2, max.overlaps = Inf, seed = 20261013, show.legend = FALSE) + labs(x = sprintf("DA 1 (%.1f%%)", 
        HG_DA1), y = sprintf("DA 2 (%.1f%%)", HG_DA2)) + guides(colour = "none") + base_theme + theme(legend.position = "none")
}
base_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), 
    axis.text = element_text(size = 8.5, colour = "black"), axis.line = element_line(linewidth = 0.4, colour = "black"), 
    axis.ticks = element_line(linewidth = 0.35, colour = "black"), plot.title = element_blank(), plot.subtitle = element_blank(), 
    panel.grid = element_blank(), legend.position = "top", legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_text(size = 8.5), 
    legend.text = element_text(size = 8), plot.margin = margin(t = 4, r = 5, b = 5, l = 5, unit = "mm"))
