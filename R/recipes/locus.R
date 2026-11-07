# Adapted from the user's original locus final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
make_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    if (mode == "color") {
        cols <- c(`Single-copy locus` = "#2C7FB8", `Multicopy pattern state` = "#D95F02")
    }
    else {
        cols <- c(`Single-copy locus` = "grey30", `Multicopy pattern state` = "grey65")
    }
    ggplot(plot_df, aes(x = diversity, y = locus)) + geom_segment(aes(x = 0, xend = diversity, yend = locus), linewidth = 0.35, 
        colour = "grey75") + geom_point(aes(fill = system_type, shape = system_type), size = 3, stroke = 0.55, colour = "black") + 
        scale_fill_manual(values = cols) + scale_shape_manual(values = c(`Single-copy locus` = 21, `Multicopy pattern state` = 22)) + 
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) + labs(x = "Diversity (PD = 1 - MP)", 
        y = "Y-STR system") + theme_pub
}
theme_pub <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 10, colour = "black"), 
    axis.text = element_text(size = 8.5, colour = "black"), legend.position = "top", legend.title = element_blank(), plot.title = element_blank(), 
    plot.subtitle = element_blank())
