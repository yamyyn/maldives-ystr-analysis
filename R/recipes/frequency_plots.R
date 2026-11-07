# Original final allele-frequency figure functions.
numeric_allele_sort_key <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
}
multicopy_component_sort_key <- function(x, component) {
    vapply(strsplit(as.character(x), "-", fixed = TRUE), function(parts) {
        nums <- suppressWarnings(as.numeric(parts))
        nums <- nums[!is.na(nums)]
        nums <- sort(nums)
        if (length(nums) < component) {
            return(Inf)
        }
        nums[component]
    }, numeric(1))
}
calculate_y_max <- function(values) {
    observed_max <- max(values, na.rm = TRUE)
    if (!is.finite(observed_max)) {
        return(Y_ROUND_STEP)
    }
    y_max <- ceiling((observed_max * Y_HEADROOM_MULTIPLIER)/Y_ROUND_STEP) * Y_ROUND_STEP
    y_max <- max(y_max, Y_ROUND_STEP)
    min(y_max, 1)
}
calculate_y_breaks <- function(y_max) {
    by_value <- if (y_max <= 0.25) {
        0.05
    }
    else if (y_max <= 0.5) {
        0.1
    }
    else {
        0.2
    }
    seq(0, y_max + 1e-12, by = by_value)
}
strip_plot_key <- function(x) {
    sub("^.*___", "", x)
}
prepare_single_plot_data <- function(tbl) {
    out <- tbl %>% mutate(Locus = factor(Locus, levels = LOCUS_ORDER), Allele = as.character(Allele), allele_numeric = numeric_allele_sort_key(Allele)) %>% 
        arrange(Locus, allele_numeric, Allele) %>% mutate(Allele_Plot_Key = paste(as.character(Locus), Allele, sep = "___"))
    out$Allele_Plot <- factor(out$Allele_Plot_Key, levels = unique(out$Allele_Plot_Key))
    out
}
prepare_multicopy_plot_data <- function(tbl) {
    out <- tbl %>% mutate(Marker = factor(Marker, levels = MULTICOPY_ORDER), Allele_Combination = as.character(Allele_Combination), 
        sort_1 = multicopy_component_sort_key(Allele_Combination, 1), sort_2 = multicopy_component_sort_key(Allele_Combination, 
            2), sort_3 = multicopy_component_sort_key(Allele_Combination, 3)) %>% arrange(Marker, sort_1, sort_2, sort_3, 
        Allele_Combination) %>% mutate(Combination_Plot_Key = paste(as.character(Marker), Allele_Combination, sep = "___"))
    out$Combination_Plot <- factor(out$Combination_Plot_Key, levels = unique(out$Combination_Plot_Key))
    out
}
theme_ystr_frequency <- function(base_size = 8.5) {
    theme_bw(base_family = FIGURE_FONT, base_size = base_size) + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = base_size), axis.title = element_text(size = base_size + 
            0.5), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size - 1.5), axis.text.y = element_text(size = base_size - 
            1), plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2), plot.subtitle = element_text(hjust = 0.5, 
            size = base_size), plot.margin = margin(5, 5, 5, 5))
}
make_single_copy_plot <- function(tbl, descriptive_title = NULL, subtitle = NULL, publication_mode = FALSE, style = c("grayscale", 
    "colour")) {
    style <- match.arg(style)
    plot_data <- prepare_single_plot_data(tbl)
    y_max <- calculate_y_max(plot_data$Allele_Frequency)
    y_breaks <- calculate_y_breaks(y_max)
    single_fill <- if (style == "colour") {
        COLOUR_SINGLE
    }
    else {
        "black"
    }
    p <- ggplot(plot_data, aes(x = Allele_Plot, y = Allele_Frequency)) + geom_col(width = 0.78, fill = single_fill) + facet_wrap(~Locus, 
        ncol = 4, scales = "free_x") + scale_x_discrete(labels = strip_plot_key) + scale_y_continuous(limits = c(0, y_max), 
        breaks = y_breaks, labels = label_number(accuracy = 0.01, trim = TRUE), expand = expansion(mult = c(0, 0.01))) + 
        labs(x = "Allele", y = "Allele frequency") + theme_ystr_frequency(base_size = 8.3)
    if (!publication_mode) {
        p <- p + labs(title = descriptive_title, subtitle = subtitle)
    }
    p
}
make_multicopy_plot <- function(tbl, descriptive_title = NULL, subtitle = NULL, publication_mode = FALSE, style = c("grayscale", 
    "colour")) {
    style <- match.arg(style)
    plot_data <- prepare_multicopy_plot_data(tbl)
    y_max <- calculate_y_max(plot_data$Combination_Frequency)
    y_breaks <- calculate_y_breaks(y_max)
    if (style == "colour") {
        p <- ggplot(plot_data, aes(x = Combination_Plot, y = Combination_Frequency, fill = Marker)) + geom_col(width = 0.78, 
            show.legend = FALSE) + scale_fill_manual(values = c(DYS385 = COLOUR_DYS385, DYF387S1 = COLOUR_DYF387S1), drop = FALSE)
    }
    else {
        p <- ggplot(plot_data, aes(x = Combination_Plot, y = Combination_Frequency)) + geom_col(width = 0.78, fill = "black")
    }
    p <- p + facet_wrap(~Marker, ncol = 1, scales = "free_x") + scale_x_discrete(labels = strip_plot_key) + scale_y_continuous(limits = c(0, 
        y_max), breaks = y_breaks, labels = label_number(accuracy = 0.01, trim = TRUE), expand = expansion(mult = c(0, 0.01))) + 
        labs(x = "Unordered allele pattern", y = "Pattern frequency") + theme_ystr_frequency(base_size = 8.3) + theme(axis.text.x = element_text(angle = 90, 
        hjust = 1, vjust = 0.5, size = 6.4))
    if (!publication_mode) {
        p <- p + labs(title = descriptive_title, subtitle = subtitle)
    }
    p
}
make_combined_plot <- function(single_tbl, multi_tbl, style = c("grayscale", "colour")) {
    style <- match.arg(style)
    p_single <- make_single_copy_plot(single_tbl, publication_mode = TRUE, style = style)
    p_multi <- make_multicopy_plot(multi_tbl, publication_mode = TRUE, style = style)
    combined <- p_single/p_multi + plot_layout(heights = c(5, 1.55)) + plot_annotation(tag_levels = "a", tag_prefix = "(", 
        tag_suffix = ")", theme = theme(plot.tag = element_text(family = FIGURE_FONT, face = "bold", size = 10)))
    combined
}
