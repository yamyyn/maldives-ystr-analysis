# Adapted from the user's original nocall final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
make_panel_a <- function(color_mode = TRUE) {
    fill_vals <- if (color_mode) {
        c(`Allele observed` = "#f2f2f2", `No allele call` = "#d95f02")
    }
    else {
        c(`Allele observed` = "#f2f2f2", `No allele call` = "#4d4d4d")
    }
    ggplot(shared_pattern, aes(x = locus, y = pattern, fill = locus_call_state)) + geom_tile(color = "white", linewidth = 0.45, 
        height = 0.7) + scale_fill_manual(values = fill_vals, name = "Locus state") + labs(x = "Y-STR locus", y = NULL) + 
        base_theme + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8), axis.text.y = element_text(size = 9), 
        axis.ticks.y = element_blank(), axis.line.y = element_blank(), legend.position = "top", plot.margin = margin(5, 5, 
            5, 5))
}
make_panel_b <- function(color_mode = TRUE) {
    hg_levels <- sort(unique(bar_dat$final_haplogroup))
    if (color_mode) {
        base_cols <- style$haplogroup_colors
        fill_vals <- base_cols[intersect(names(base_cols), hg_levels)]
        missing_hg <- setdiff(hg_levels, names(fill_vals))
        if (length(missing_hg) > 0) {
            extras <- hue_pal()(length(missing_hg))
            names(extras) <- missing_hg
            fill_vals <- c(fill_vals, extras)
        }
    }
    else {
        fill_vals <- setNames(c("#4d4d4d", "#969696", "#bdbdbd", "#d9d9d9")[seq_along(hg_levels)], hg_levels)
    }
    ggplot(bar_dat, aes(x = atoll_region, y = n, fill = final_haplogroup)) + geom_col(width = 0.75, color = "black", linewidth = 0.25) + 
        geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.3, color = "white") + scale_fill_manual(values = fill_vals, 
        name = "Final haplogroup") + scale_y_continuous(breaks = style$nocall_breaks, expand = expansion(mult = c(0, 0.03))) + 
        labs(x = "Administrative division (region)", y = "Number of recurrent-null profiles") + coord_flip(clip = "off") + 
        base_theme + theme(legend.position = "top", axis.text.y = element_text(size = 9))
}
make_combined <- function(color_mode = TRUE) {
    p1 <- make_panel_a(color_mode)
    p2 <- make_panel_b(color_mode)
    (p1/p2) + plot_layout(heights = c(0.42, 0.58)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", 
        size = 13))
}
base_theme <- theme_classic(base_size = 12) + theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black", 
    face = "bold"), legend.title = element_text(face = "bold"))
