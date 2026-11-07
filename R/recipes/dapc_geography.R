# Adapted from the user's original dapc_geography final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_region_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    palette <- if (mode == "color") {
        REGION_COLORS
    }
    else {
        REGION_GRAYS
    }
    p <- ggplot(region_df, aes(x = DA1, y = DA2)) + geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.3) + geom_vline(xintercept = 0, 
        colour = "grey82", linewidth = 0.3) + geom_point(aes(shape = Group), colour = "grey55", alpha = 0.62, size = 1.55, 
        stroke = 0.45) + geom_point(data = region_centroids, aes(shape = Group, fill = Group, colour = Group), size = 5.2, 
        stroke = 1.05) + ggrepel::geom_label_repel(data = region_centroids, aes(label = Group), colour = "black", fill = "white", 
        family = style$font, size = 3, label.size = 0.2, box.padding = 0.28, point.padding = 0.25, min.segment.length = 0, 
        segment.colour = "grey65", segment.size = 0.22, seed = 20261011, max.overlaps = Inf, show.legend = FALSE) + scale_shape_manual(values = REGION_SHAPES, 
        drop = FALSE) + scale_colour_manual(values = palette, drop = FALSE) + scale_fill_manual(values = palette, drop = FALSE) + 
        labs(x = paste0("DA 1 (", sprintf("%.1f", REGION_DA1), "%)"), y = paste0("DA 2 (", sprintf("%.1f", REGION_DA2), "%)"), 
            shape = "Region", colour = "Region", fill = "Region") + guides(colour = "none", fill = "none", shape = guide_legend(override.aes = list(size = 4.2, 
        colour = if (mode == "color") {
            unname(palette[REGION_GROUPS])
        } else {
            unname(palette[REGION_GROUPS])
        }, fill = unname(palette[REGION_GROUPS])))) + base_theme
    p
}
build_admin_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    palette <- if (mode == "color") {
        REGION_COLORS
    }
    else {
        REGION_GRAYS
    }
    p <- ggplot(admin_df, aes(x = DA1, y = DA2)) + geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.3) + geom_vline(xintercept = 0, 
        colour = "grey82", linewidth = 0.3) + geom_point(colour = "grey82", alpha = 0.62, size = 1.15, shape = 1, stroke = 0.38) + 
        geom_point(data = admin_centroids, aes(shape = Region, fill = Region, colour = Region), size = 3.6, stroke = 0.95) + 
        ggrepel::geom_text_repel(data = admin_centroids, aes(label = Group), family = style$font, colour = "black", size = 2.75, 
            box.padding = 0.2, point.padding = 0.15, min.segment.length = 0.05, segment.colour = "grey70", segment.size = 0.2, 
            max.overlaps = Inf, seed = 20261012, show.legend = FALSE) + scale_shape_manual(values = REGION_SHAPES, drop = FALSE) + 
        scale_colour_manual(values = palette, drop = FALSE) + scale_fill_manual(values = palette, drop = FALSE) + labs(x = paste0("DA 1 (", 
        sprintf("%.1f", ADMIN_DA1), "%)"), y = paste0("DA 2 (", sprintf("%.1f", ADMIN_DA2), "%)"), shape = "Region", colour = "Region", 
        fill = "Region") + guides(colour = "none", fill = "none", shape = guide_legend(override.aes = list(size = 4.2, colour = unname(palette[REGION_GROUPS]), 
        fill = unname(palette[REGION_GROUPS])))) + base_theme
    p
}
base_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), 
    axis.text = element_text(size = 8.5, colour = "black"), axis.line = element_line(linewidth = 0.4, colour = "black"), 
    axis.ticks = element_line(linewidth = 0.35, colour = "black"), plot.title = element_blank(), plot.subtitle = element_blank(), 
    panel.grid = element_blank(), legend.position = "top", legend.direction = "horizontal", legend.box = "horizontal", legend.title = element_text(size = 8.5), 
    legend.text = element_text(size = 8), plot.margin = margin(t = 4, r = 5, b = 5, l = 5, unit = "mm"))
