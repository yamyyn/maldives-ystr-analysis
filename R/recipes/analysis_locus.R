# Adapted from the user's original analysis_locus final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
clean_state <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NA", "NaN", ".", "NULL")] <- NA_character_
    x
}
format_numeric_state <- function(x) {
    x <- clean_state(x)
    out <- rep(NA_character_, length(x))
    keep <- !is.na(x)
    if (any(keep)) {
        num <- suppressWarnings(as.numeric(x[keep]))
        if (any(is.na(num))) {
            stop(paste("Non-numeric allele state(s):", paste(unique(x[keep][is.na(num)]), collapse = ", ")), call. = FALSE)
        }
        formatted <- format(num, scientific = FALSE, trim = TRUE, digits = 15)
        formatted <- sub("\\.0+$", "", formatted)
        formatted <- sub("(\\.[0-9]*?)0+$", "\\1", formatted)
        out[keep] <- formatted
    }
    out
}
make_unordered_pattern <- function(data, columns) {
    apply(data[, columns, drop = FALSE], 1, function(z) {
        z <- format_numeric_state(z)
        z <- z[!is.na(z)]
        if (length(z) == 0L) {
            return(NA_character_)
        }
        observed_multicopy_state(z,collapse_identical=all(grepl("^dys385",columns,ignore.case=TRUE)))
    })
}
calc_single_locus_stats <- function(x) {
    x <- format_numeric_state(x)
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2L) {
        stop("Insufficient typed observations at single-copy locus.", call. = FALSE)
    }
    numeric_x <- as.numeric(x)
    tab <- table(x)
    p <- as.numeric(tab)/n
    mp <- sum(p^2)
    gd <- (n/(n - 1L)) * (1 - mp)
    pd <- 1 - mp
    mcaf <- max(p)
    ne <- 1/mp
    tibble(N = n, Na_Np = length(tab), Observed_Range = paste0(format(min(numeric_x), scientific = FALSE, trim = TRUE), "\342\200\223", 
        format(max(numeric_x), scientific = FALSE, trim = TRUE)), MCAF = mcaf, GD = gd, MP = mp, PD = pd, Ne = ne)
}
calc_multicopy_stats <- function(pattern, allele_columns) {
    pattern <- as.character(pattern)
    keep <- !is.na(pattern) & pattern != ""
    pattern <- pattern[keep]
    n <- length(pattern)
    if (n < 2L) {
        stop("Insufficient typed observations at multicopy system.", call. = FALSE)
    }
    tab <- table(pattern)
    p <- as.numeric(tab)/n
    mp <- sum(p^2)
    gd <- (n/(n - 1L)) * (1 - mp)
    pd <- 1 - mp
    mcaf <- max(p)
    ne <- 1/mp
    allele_values <- unlist(lapply(allele_columns, function(x) {
        as.numeric(format_numeric_state(x))
    }), use.names = FALSE)
    allele_values <- allele_values[is.finite(allele_values)]
    tibble(N = n, Na_Np = length(tab), Observed_Range = paste0(format(min(allele_values), scientific = FALSE, trim = TRUE), 
        "\342\200\223", format(max(allele_values), scientific = FALSE, trim = TRUE)), MCAF = mcaf, GD = gd, MP = mp, PD = pd, 
        Ne = ne)
}
