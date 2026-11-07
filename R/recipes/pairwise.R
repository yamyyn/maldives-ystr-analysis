# Adapted from the user's original pairwise final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
make_lower_long <- function(mat, order_vec) {
    n <- length(order_vec)
    grid <- expand.grid(row_index = seq_len(n), col_index = seq_len(n)) %>% filter(row_index > col_index) %>% mutate(row_group = order_vec[row_index], 
        col_group = order_vec[col_index], rst = mapply(function(r, c) {
            mat[r, c]
        }, row_group, col_group))
    grid
}
build_admin_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    min_rst <- min(admin_df$rst, na.rm = TRUE)
    max_rst <- max(admin_df$rst, na.rm = TRUE)
    admin_legend_breaks <- c(-0.025, 0, 0.05, 0.1, 0.135)
    p <- ggplot(admin_df, aes(x = col_index, y = row_index, fill = rst)) + geom_tile(width = 0.96, height = 0.96) + geom_segment(data = tibble::tibble(b = ADMIN_BLOCK_ENDS), 
        aes(x = b, xend = b, y = b, yend = length(ADMIN_ORDER) + 0.5), inherit.aes = FALSE, linewidth = 0.35, colour = "grey55") + 
        geom_segment(data = tibble::tibble(b = ADMIN_BLOCK_ENDS), aes(x = 0.5, xend = b, y = b, yend = b), inherit.aes = FALSE, 
            linewidth = 0.35, colour = "grey55") + geom_text(data = admin_df %>% filter(bonf_significant), aes(label = "*"), 
        colour = "black", size = 4.3, fontface = "bold") + scale_x_continuous(breaks = seq_along(ADMIN_ORDER), labels = ADMIN_ORDER, 
        limits = c(0.5, length(ADMIN_ORDER) + 0.5), expand = c(0, 0)) + scale_y_continuous(breaks = seq_along(ADMIN_ORDER), 
        labels = ADMIN_ORDER, limits = c(length(ADMIN_ORDER) + 0.5, 0.5), expand = c(0, 0)) + coord_fixed(clip = "off") + 
        labs(x = "Administrative division", y = "Administrative division", fill = expression(R[ST])) + heatmap_theme
    if (mode == "color") {
        p <- p + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(min_rst, 
            max_rst), oob = scales::squish, breaks = admin_legend_breaks, labels = scales::label_number(accuracy = 0.001))
    }
    else {
        p <- p + scale_fill_gradient2(low = "grey70", mid = "white", high = "grey15", midpoint = 0, limits = c(min_rst, max_rst), 
            oob = scales::squish, breaks = admin_legend_breaks, labels = scales::label_number(accuracy = 0.001))
    }
    p
}
build_region_plot <- function(mode = c("color", "grayscale")) {
    mode <- match.arg(mode)
    min_rst <- min(region_df$rst, na.rm = TRUE)
    max_rst <- max(region_df$rst, na.rm = TRUE)
    max_abs_rst <- max(abs(c(min_rst, max_rst)))
    region_plot_df <- region_df %>% mutate(text_colour = if (mode == "grayscale") {
        ifelse(rst > 0.0045, "white", "black")
    }
    else {
        ifelse(abs(rst) > 0.65 * max_abs_rst, "white", "black")
    }, label = sprintf("%.4f", rst))
    p <- ggplot(region_plot_df, aes(x = col_index, y = row_index, fill = rst)) + geom_tile(width = 0.96, height = 0.96) + 
        geom_text(aes(label = label, colour = text_colour), size = 4) + scale_colour_identity() + scale_x_continuous(breaks = seq_along(REGION_ORDER), 
        labels = REGION_ORDER, limits = c(0.5, length(REGION_ORDER) + 0.5), expand = c(0, 0)) + scale_y_continuous(breaks = seq_along(REGION_ORDER), 
        labels = REGION_ORDER, limits = c(length(REGION_ORDER) + 0.5, 0.5), expand = c(0, 0)) + coord_fixed(clip = "off") + 
        labs(x = "Geographic region", y = "Geographic region", fill = expression(R[ST])) + heatmap_theme + theme(axis.text.x = element_text(angle = 0, 
        hjust = 0.5))
    if (mode == "color") {
        p <- p + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(min_rst, 
            max_rst), oob = scales::squish, labels = scales::label_number(accuracy = 0.001))
    }
    else {
        p <- p + scale_fill_gradient2(low = "grey70", mid = "white", high = "grey15", midpoint = 0, limits = c(min_rst, max_rst), 
            oob = scales::squish, labels = scales::label_number(accuracy = 0.001))
    }
    p
}
heatmap_theme <- theme_classic(base_family = style$font, base_size = 9) + theme(axis.title = element_text(size = 10, colour = "black"), 
    axis.text = element_text(size = 8.5, colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
    axis.line = element_line(colour = "black", linewidth = 0.4), axis.ticks = element_blank(), panel.grid = element_blank(), 
    legend.title = element_text(size = 9, colour = "black"), legend.text = element_text(size = 8, colour = "black"), plot.title = element_blank(), 
    plot.subtitle = element_blank(), plot.margin = margin(4, 5, 4, 5, unit = "mm"))
