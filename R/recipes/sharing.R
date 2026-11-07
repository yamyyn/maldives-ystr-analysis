# Adapted from the user's original sharing final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
matrix_to_long <- function(m, group_order) {
    tibble(Row = rep(rownames(m), times = ncol(m)), Col = rep(colnames(m), each = nrow(m)), Shared = as.vector(m)) %>% mutate(RowIndex = match(Row, 
        group_order), ColIndex = match(Col, group_order))
}
build_admin_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    p <- ggplot(admin_long, aes(x = Col, y = Row, fill = Shared)) + geom_tile(colour = "white", linewidth = 0.35, width = 0.98, 
        height = 0.98) + geom_vline(xintercept = style$sharing_vlines, colour = "grey65", linewidth = 0.35) + geom_hline(yintercept = style$sharing_hlines, 
        colour = "grey65", linewidth = 0.35) + scale_x_discrete(limits = ADMIN_ORDER, drop = FALSE) + scale_y_discrete(limits = rev(ADMIN_ORDER), 
        drop = FALSE) + coord_fixed(clip = "off") + labs(x = "Administrative division", y = "Administrative division", fill = "Shared distinct\nhaplotypes") + 
        publication_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.border = element_rect(colour = "black", 
        fill = NA, linewidth = 0.35), legend.position = "right")
    if (mode == "color") {
        p <- p + scale_fill_gradientn(colours = c(HEAT_LOW, HEAT_MID, HEAT_HIGH), values = scales::rescale(c(0, 1, 2)), limits = c(0, 
            2), breaks = c(0, 1, 2), oob = scales::squish)
    }
    else {
        p <- p + scale_fill_gradient(low = "grey92", high = "grey15", limits = c(0, 2), breaks = c(0, 1, 2), oob = scales::squish)
    }
    p
}
build_region_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    p <- ggplot(region_long, aes(x = Col, y = Row, fill = Shared)) + geom_tile(colour = "white", linewidth = 0.55, width = 0.98, 
        height = 0.98) + geom_text(aes(label = Label, colour = LabelColor), family = style$font, size = 4) + scale_colour_identity() + 
        scale_x_discrete(limits = REGION_ORDER, drop = FALSE) + scale_y_discrete(limits = rev(REGION_ORDER), drop = FALSE) + 
        coord_fixed(clip = "off") + labs(x = "Geographic region", y = "Geographic region", fill = "Shared distinct\nhaplotypes") + 
        publication_theme + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4), legend.position = "right")
    if (mode == "color") {
        p <- p + scale_fill_gradientn(colours = c(HEAT_LOW, HEAT_MID, HEAT_HIGH), values = scales::rescale(c(0, 3, 6)), limits = c(0, 
            6), breaks = 0:6, oob = scales::squish)
    }
    else {
        p <- p + scale_fill_gradient(low = "grey92", high = "grey15", limits = c(0, 6), breaks = 0:6, oob = scales::squish)
    }
    p
}
make_combined_grob <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    admin_grob <- ggplotGrob(build_admin_plot(mode))
    region_grob <- ggplotGrob(build_region_plot(mode))
    lay <- grid.layout(nrow = 4, ncol = 1, heights = unit(c(0.035, 0.55, 0.035, 0.38), "npc"))
    grobTree(rectGrob(gp = gpar(fill = "white", col = NA)), grobTree(textGrob("(a)", x = unit(0.01, "npc"), y = unit(0.5, 
        "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5)), vp = viewport(layout.pos.row = 1, 
        layout.pos.col = 1)), grobTree(admin_grob, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)), grobTree(textGrob("(b)", 
        x = unit(0.01, "npc"), y = unit(0.5, "npc"), just = c("left", "center"), gp = gpar(fontfamily = style$font, fontsize = 9.5)), 
        vp = viewport(layout.pos.row = 3, layout.pos.col = 1)), grobTree(region_grob, vp = viewport(layout.pos.row = 4, layout.pos.col = 1)), 
        vp = viewport(layout = lay))
}
publication_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, 
    colour = "black"), axis.text = element_text(size = 8.2, colour = "black"), axis.line = element_blank(), axis.ticks = element_blank(), 
    panel.grid = element_blank(), plot.title = element_blank(), plot.subtitle = element_blank(), plot.margin = margin(t = 4, 
        r = 5, b = 5, l = 5, unit = "mm"), legend.title = element_text(size = 8.2, colour = "black"), legend.text = element_text(size = 7.8, 
        colour = "black"))
