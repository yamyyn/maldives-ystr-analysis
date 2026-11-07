# Adapted from the user's original table_alleles final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
format_allele_state <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & !is.finite(x_num))) {
        stop("Non-numeric allele state encountered.", call. = FALSE)
    }
    out <- sprintf("%.10f", x_num)
    out <- sub("0+$", "", out)
    out <- sub("\\.$", "", out)
    out[is.na(x_num)] <- NA_character_
    out
}
make_pattern <- function(data, columns) {
    apply(data[, columns, drop = FALSE], 1, function(z) {
        z <- z[!is.na(z) & trimws(z) != ""]
        if (length(z) == 0L) {
            return(NA_character_)
        }
        z <- format_allele_state(z)
        observed_multicopy_state(z,collapse_identical=identical(columns,DYS385_COLUMNS))
    })
}
frequency_vector <- function(x) {
    x <- x[!is.na(x) & trimws(as.character(x)) != ""]
    n <- length(x)
    if (n == 0L) {
        return(list(frequency = numeric(0), n = 0L))
    }
    counts <- table(x, useNA = "no")
    freq <- as.numeric(counts)/n
    names(freq) <- names(counts)
    list(frequency = freq, n = n)
}
metrics_from_frequency <- function(freq, n) {
    if (length(freq) == 0L || n == 0L) {
        return(tibble(DC = NA_real_, GD = NA_real_, MCAF = NA_real_, MP = NA_real_, N = 0L, Na = 0L, Ne = NA_real_))
    }
    mp <- sum(freq^2)
    dc <- 1 - mp
    gd <- if (n > 1L) {
        (n/(n - 1)) * dc
    }
    else {
        NA_real_
    }
    tibble(DC = dc, GD = gd, MCAF = max(freq), MP = mp, N = as.integer(n), Na = as.integer(length(freq)), Ne = if (mp > 0) {
        1/mp
    }
    else {
        NA_real_
    })
}
pattern_sort_key <- function(x) {
    split <- strsplit(x, "-", fixed = TRUE)
    max_len <- max(lengths(split))
    mat <- do.call(rbind, lapply(split, function(z) {
        c(as.numeric(z), rep(Inf, max_len - length(z)))
    }))
    do.call(order, as.data.frame(mat))
}
build_population_table <- function(data) {
    single_freq <- vector("list", length(SINGLE_COPY_LOCI))
    names(single_freq) <- SINGLE_COPY_LOCI
    single_metrics <- vector("list", length(SINGLE_COPY_LOCI))
    names(single_metrics) <- SINGLE_COPY_LOCI
    allele_union <- character(0)
    for (locus in SINGLE_COPY_LOCI) {
        states <- format_allele_state(data[[locus]])
        fv <- frequency_vector(states)
        single_freq[[locus]] <- fv$frequency
        single_metrics[[locus]] <- metrics_from_frequency(fv$frequency, fv$n)
        allele_union <- union(allele_union, names(fv$frequency))
    }
    allele_union <- allele_union[order(as.numeric(allele_union))]
    dyf_pattern <- make_pattern(data, DYF387S1_COLUMNS)
    dys385_pattern <- make_pattern(data, DYS385_COLUMNS)
    dyf_fv <- frequency_vector(dyf_pattern)
    dys385_fv <- frequency_vector(dys385_pattern)
    dyf_metrics <- metrics_from_frequency(dyf_fv$frequency, dyf_fv$n)
    dys385_metrics <- metrics_from_frequency(dys385_fv$frequency, dys385_fv$n)
    dyf_patterns <- names(dyf_fv$frequency)
    dys385_patterns <- names(dys385_fv$frequency)
    if (length(dyf_patterns) > 1L) {
        dyf_patterns <- dyf_patterns[pattern_sort_key(dyf_patterns)]
    }
    if (length(dys385_patterns) > 1L) {
        dys385_patterns <- dys385_patterns[pattern_sort_key(dys385_patterns)]
    }
    n_body <- max(length(allele_union), length(dyf_patterns), length(dys385_patterns))
    body <- matrix(NA_character_, nrow = n_body, ncol = 28)
    colnames(body) <- c("Allele", SINGLE_COPY_LABELS, "Pattern_DYF387S1", "DYF387S1", "Pattern_DYS385", "DYS385")
    for (i in seq_along(allele_union)) {
        allele <- allele_union[i]
        body[i, "Allele"] <- allele
        for (j in seq_along(SINGLE_COPY_LOCI)) {
            locus <- SINGLE_COPY_LOCI[j]
            if (allele %in% names(single_freq[[locus]])) {
                body[i, SINGLE_COPY_LABELS[j]] <- sprintf("%.4f", single_freq[[locus]][[allele]])
            }
        }
    }
    if (length(dyf_patterns) > 0L) {
        for (i in seq_along(dyf_patterns)) {
            pattern <- dyf_patterns[i]
            body[i, "Pattern_DYF387S1"] <- pattern
            body[i, "DYF387S1"] <- sprintf("%.4f", dyf_fv$frequency[[pattern]])
        }
    }
    if (length(dys385_patterns) > 0L) {
        for (i in seq_along(dys385_patterns)) {
            pattern <- dys385_patterns[i]
            body[i, "Pattern_DYS385"] <- pattern
            body[i, "DYS385"] <- sprintf("%.4f", dys385_fv$frequency[[pattern]])
        }
    }
    body <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
    summary_metric_order <- c("DC", "GD", "MCAF", "MP", "N", "Na", "Ne")
    summary_rows <- vector("list", length(summary_metric_order))
    names(summary_rows) <- summary_metric_order
    for (metric in summary_metric_order) {
        row <- rep(NA_character_, 28)
        row[1] <- metric
        for (j in seq_along(SINGLE_COPY_LOCI)) {
            value <- single_metrics[[SINGLE_COPY_LOCI[j]]][[metric]]
            row[1 + j] <- if (metric %in% c("N", "Na")) {
                as.character(as.integer(value))
            }
            else {
                sprintf("%.4f", as.numeric(value))
            }
        }
        row[25] <- metric
        row[26] <- if (metric %in% c("N", "Na")) {
            as.character(as.integer(dyf_metrics[[metric]]))
        }
        else {
            sprintf("%.4f", as.numeric(dyf_metrics[[metric]]))
        }
        row[27] <- metric
        row[28] <- if (metric %in% c("N", "Na")) {
            as.character(as.integer(dys385_metrics[[metric]]))
        }
        else {
            sprintf("%.4f", as.numeric(dys385_metrics[[metric]]))
        }
        summary_rows[[metric]] <- row
    }
    summary_df <- as.data.frame(do.call(rbind, summary_rows), stringsAsFactors = FALSE, check.names = FALSE)
    names(summary_df) <- names(body)
    display_header <- c("Allele", SINGLE_COPY_LABELS, "Pattern", "DYF387S1", "Pattern", "DYS385")
    list(body = body, summary = summary_df, header = display_header, single_metrics = single_metrics, dyf_metrics = dyf_metrics, 
        dys385_metrics = dys385_metrics)
}
write_population_sheet <- function(wb, sheet_name, population_name, population_data, table_object) {
    addWorksheet(wb, sheet_name, gridLines = TRUE)
    n_pop <- nrow(population_data)
    display_name <- if (population_name == study$population_label) {
        study$population_label
    }
    else {
        paste0("Administrative Division ", population_name)
    }
    mergeCells(wb, sheet_name, cols = 1:28, rows = 1)
    writeData(wb, sheet_name, x = paste0("Allele Frequency Distribution \342\200\224 ", display_name, " (N=", n_pop, ")"), 
        startCol = 1, startRow = 1, colNames = FALSE)
    addStyle(wb, sheet_name, style = style_title, rows = 1, cols = 1:28, gridExpand = TRUE)
    setRowHeights(wb, sheet_name, rows = 1, heights = 24)
    mergeCells(wb, sheet_name, cols = 1:28, rows = 2)
    writeData(wb, sheet_name, x = "Single-copy allele frequencies and multicopy state-pattern frequencies; frequencies shown to four decimals", 
        startCol = 1, startRow = 2, colNames = FALSE)
    addStyle(wb, sheet_name, style = style_subtitle, rows = 2, cols = 1:28, gridExpand = TRUE)
    writeData(wb, sheet_name, x = as.data.frame(t(table_object$header), stringsAsFactors = FALSE, check.names = FALSE), startCol = 1, 
        startRow = 4, colNames = FALSE, rowNames = FALSE)
    addStyle(wb, sheet_name, style = style_header, rows = 4, cols = 1:28, gridExpand = TRUE)
    body_start <- 5L
    body_end <- body_start + nrow(table_object$body) - 1L
    writeData(wb, sheet_name, x = table_object$body, startCol = 1, startRow = body_start, colNames = FALSE, rowNames = FALSE, 
        keepNA = FALSE)
    addStyle(wb, sheet_name, style = style_body, rows = body_start:body_end, cols = 1:28, gridExpand = TRUE)
    addStyle(wb, sheet_name, style = style_text, rows = body_start:body_end, cols = c(25, 27), gridExpand = TRUE, stack = TRUE)
    summary_start <- body_end + 1L
    summary_end <- summary_start + nrow(table_object$summary) - 1L
    writeData(wb, sheet_name, x = table_object$summary, startCol = 1, startRow = summary_start, colNames = FALSE, rowNames = FALSE, 
        keepNA = FALSE)
    addStyle(wb, sheet_name, style = style_summary, rows = summary_start:summary_end, cols = 1:28, gridExpand = TRUE)
    addStyle(wb, sheet_name, style = style_text, rows = summary_start:summary_end, cols = c(25, 27), gridExpand = TRUE, stack = TRUE)
    note_row <- summary_end + 2L
    mergeCells(wb, sheet_name, cols = 1:28, rows = note_row)
    writeData(wb, sheet_name, x = study$allele_note, startCol = 1, startRow = note_row, colNames = FALSE)
    addStyle(wb, sheet_name, style = style_note, rows = note_row, cols = 1:28, gridExpand = TRUE)
    setRowHeights(wb, sheet_name, rows = note_row, heights = 42)
    setColWidths(wb, sheet_name, cols = 1, widths = 9)
    setColWidths(wb, sheet_name, cols = 2:24, widths = 10)
    setColWidths(wb, sheet_name, cols = 25, widths = 14)
    setColWidths(wb, sheet_name, cols = 26, widths = 11)
    setColWidths(wb, sheet_name, cols = 27, widths = 14)
    setColWidths(wb, sheet_name, cols = 28, widths = 11)
    freezePane(wb, sheet_name, firstActiveRow = 5, firstActiveCol = 2)
    invisible(list(body_start = body_start, body_end = body_end, summary_start = summary_start, summary_end = summary_end, 
        note_row = note_row))
}
style_title <- createStyle(fontName = style$font, fontSize = 12, fontColour = "#FFFFFF", fgFill = "#1F4E78", textDecoration = "bold", 
    halign = "center", valign = "center")
style_subtitle <- createStyle(fontName = style$font, fontSize = 9, fgFill = "#D9EAF7", textDecoration = "italic", halign = "center", 
    valign = "center")
style_header <- createStyle(fontName = style$font, fontSize = 9, fontColour = "#FFFFFF", fgFill = "#5B9BD5", textDecoration = "bold", 
    halign = "center", valign = "center", border = "TopBottomLeftRight", borderColour = "#D9E1F2")
style_body <- createStyle(fontName = style$font, fontSize = 8, halign = "center", valign = "center")
style_text <- createStyle(numFmt = "TEXT")
style_summary <- createStyle(fontName = style$font, fontSize = 8, fgFill = "#E2F0D9", textDecoration = "bold", halign = "center", 
    valign = "center", border = "TopBottom", borderColour = "#A9D18E")
style_note <- createStyle(fontName = style$font, fontSize = 8, fgFill = "#FFF2CC", textDecoration = "italic", wrapText = TRUE, 
    valign = "top")
