# Adapted from the user's original panels final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
build_metric_panel <- function(metric, panel_title, tag, mode = c("color", "grayscale"), y_limits, y_breaks, label_digits = 3, 
    show_x_text = FALSE) {
    mode <- match.arg(mode)
    if (mode == "color") {
        fills <- c(`Reference / assay progression` = "#4C78A8", `Yfiler Plus subsets` = "#8FB6D9")
    }
    else {
        fills <- c(`Reference / assay progression` = "grey52", `Yfiler Plus subsets` = "grey72")
    }
    df <- mutate(panel_df, plot_value = .data[[metric]])
    p <- ggplot(df, aes(x = x_pos, y = plot_value, fill = group)) + geom_col(width = 0.72, colour = "black", linewidth = 0.35) + 
        geom_vline(xintercept = GROUP_SEPARATOR_X, linewidth = 0.35, linetype = "dashed", colour = "grey55") + geom_text(aes(label = sprintf(paste0("%.", 
        label_digits, "f"), plot_value)), vjust = -0.38, size = 2.7, family = style$font, colour = "black") + scale_fill_manual(values = fills, 
        drop = FALSE) + scale_x_continuous(breaks = PANEL_BREAKS, labels = PANEL_LABELS, limits = style$panel_xlim, expand = c(0, 
        0)) + scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = expansion(mult = c(0, 0.07))) + labs(x = NULL, 
        y = NULL, title = panel_title, tag = tag) + guides(fill = guide_legend(title = NULL, nrow = 1, byrow = TRUE)) + theme_pub
    if (show_x_text) {
        p <- p + theme(axis.text.x = element_text(size = 7.2, angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.9, margin = margin(t = 2)))
    }
    else {
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    p
}
build_figure <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    p_a <- build_metric_panel(metric = "dc", panel_title = "Discrimination capacity (DC)", tag = "(a)", mode = mode, y_limits = style$panel_axes$dc$limits, 
        y_breaks = style$panel_axes$dc$breaks, label_digits = 3, show_x_text = FALSE)
    p_b <- build_metric_panel(metric = "fraction_unique", panel_title = "Fraction of unique haplotypes", tag = "(b)", mode = mode, 
        y_limits = style$panel_axes$fraction_unique$limits, y_breaks = style$panel_axes$fraction_unique$breaks, label_digits = 3, 
        show_x_text = FALSE)
    p_c <- build_metric_panel(metric = "hd", panel_title = "Haplotype diversity (HD)", tag = "(c)", mode = mode, y_limits = style$panel_axes$hd$limits, 
        y_breaks = style$panel_axes$hd$breaks, label_digits = 3, show_x_text = TRUE)
    p_d <- build_metric_panel(metric = "hmp", panel_title = "Haplotype match probability (HMP)", tag = "(d)", mode = mode, 
        y_limits = style$panel_axes$hmp$limits, y_breaks = style$panel_axes$hmp$breaks, label_digits = 4, show_x_text = TRUE)
    combined <- ((p_a | p_b)/(p_c | p_d)) + plot_layout(widths = c(1, 1), heights = c(1, 1), guides = "collect") + plot_annotation(caption = "Y-STR panel configuration", 
        theme = theme(legend.position = "top", legend.direction = "horizontal", plot.caption = element_text(family = style$font, 
            size = 10, hjust = 0.5, colour = "black", margin = margin(t = 5)), plot.margin = margin(4, 4, 4, 4, unit = "mm")))
    combined
}
theme_pub <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 9.5, colour = "black"), 
    axis.text = element_text(size = 8, colour = "black"), axis.line = element_line(colour = "black", linewidth = 0.4), axis.ticks = element_line(colour = "black", 
        linewidth = 0.35), panel.grid = element_blank(), plot.title = element_text(size = 9.5, face = "plain", hjust = 0.5, 
        margin = margin(b = 5)), plot.tag = element_text(family = style$font, size = 9.5, face = "plain"), plot.tag.position = c(0.01, 
        0.99), legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 8.2), legend.key.width = grid::unit(8, 
        "mm"), plot.margin = margin(4, 5, 4, 5, unit = "mm"))
