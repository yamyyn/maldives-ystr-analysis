# Original mathematical/parsing helpers; top-level study scripts were not executed.
normalise_label <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x[is.na(x) | x == ""] <- NA_character_
    x
}
strict_optional_numeric <- function(x, field_name, allow_percent = FALSE) {
    original <- as.character(x)
    trimmed <- stringr::str_trim(original)
    blank <- is.na(trimmed) | trimmed == ""
    parse_string <- trimmed
    if (allow_percent) {
        parse_string <- stringr::str_replace_all(parse_string, "%", "")
        parse_string <- stringr::str_trim(parse_string)
    }
    parsed <- suppressWarnings(as.numeric(parse_string))
    failed <- !blank & is.na(parsed)
    if (any(failed)) {
        fail_with_log(paste("Unparseable nonblank", field_name, "value(s):", paste(unique(trimmed[failed]), collapse = ", ")))
    }
    parsed
}
numeric_summary_one <- function(x, field_name) {
    nonmissing <- x[!is.na(x)]
    if (length(nonmissing) == 0L) {
        return(tibble(Field = field_name, Nonmissing_N = 0L, Missing_N = length(x), Minimum = NA_real_, Q1 = NA_real_, Median = NA_real_, 
            Mean = NA_real_, Q3 = NA_real_, Maximum = NA_real_))
    }
    tibble(Field = field_name, Nonmissing_N = length(nonmissing), Missing_N = sum(is.na(x)), Minimum = min(nonmissing), Q1 = unname(stats::quantile(nonmissing, 
        0.25)), Median = stats::median(nonmissing), Mean = mean(nonmissing), Q3 = unname(stats::quantile(nonmissing, 0.75)), 
        Maximum = max(nonmissing))
}
