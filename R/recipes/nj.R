# Adapted from the user's original nj final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
make_nj_plot <- function(filename, device_type = c("pdf", "svg", "tiff", "png"), mode = c("color", "grayscale")) {
    device_type <- match.arg(device_type)
    mode <- match.arg(mode)
    WIDTH_MM <- 174
    HEIGHT_MM <- 150
    if (device_type == "pdf") {
        grDevices::cairo_pdf(filename = filename, width = WIDTH_MM/25.4, height = HEIGHT_MM/25.4, family = style$font, bg = "white")
    }
    if (device_type == "svg") {
        grDevices::svg(filename, width = WIDTH_MM/25.4, height = HEIGHT_MM/25.4, family = style$font, bg = "white")
    }
    if (device_type == "tiff") {
        grDevices::tiff(filename = filename, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", res = style$dpi, compression = "lzw", 
            bg = "white")
    }
    if (device_type == "png") {
        grDevices::png(filename = filename, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", res = style$dpi, bg = "white")
    }
    oldpar <- par(no.readonly = TRUE)
    on.exit({
        par(oldpar)
        grDevices::dev.off()
    }, add = TRUE)
    par(family = style$font, mar = c(5, 5.2, 1, 4), oma = c(1.2, 0, 0, 0), mgp = c(2.3, 0.75, 0), tcl = -0.25, xpd = NA)
    if (mode == "color") {
        tip_text_col <- tip_meta$Color
        tip_symbol_col <- tip_meta$Color
        tip_symbol_bg <- tip_meta$Color
        legend_col <- unname(COLOR_PALETTE[REGION_LEVELS])
    }
    else {
        tip_text_col <- rep("black", nrow(tip_meta))
        tip_symbol_col <- tip_meta$Gray
        tip_symbol_bg <- tip_meta$Gray
        legend_col <- unname(GRAYSCALE_PALETTE[REGION_LEVELS])
    }
    ape::plot.phylo(x = display_tree, type = "phylogram", use.edge.length = TRUE, show.tip.label = TRUE, tip.color = tip_text_col, 
        cex = 0.82, font = 1, edge.color = "black", edge.width = 0.75, direction = "leftwards", root.edge = FALSE, no.margin = FALSE, 
        align.tip.label = TRUE, label.offset = max(display_tree$edge.length, na.rm = TRUE) * 0.3)
    ape::tiplabels(pch = tip_meta$Shape, col = tip_symbol_col, bg = tip_symbol_bg, cex = 0.68)
    ape::axisPhylo(side = 1, cex.axis = 0.72, las = 1)
    mtext(expression(italic(R)[ST] * "-derived NJ branch length"), side = 1, line = 2.2, cex = 0.82, family = style$font)
    legend("bottom", inset = c(0, -0.19), xpd = NA, horiz = TRUE, bty = "n", legend = REGION_LEVELS, pch = unname(REGION_SHAPES[REGION_LEVELS]), 
        col = legend_col, pt.bg = legend_col, pt.cex = 1, cex = 0.76, x.intersp = 0.65, y.intersp = 0.8)
    invisible(TRUE)
}
