# Adapted from the user's original inventory_rules final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
required_formats_vector <- function(x) {
    if (is.na(x) || trimws(x) == "") {
        return(character(0))
    }
    unlist(strsplit(tolower(trimws(x)), "\\|"))
}
format_present <- function(candidates, fmt) {
    if (nrow(candidates) == 0L) {
        return(FALSE)
    }
    fmt <- tolower(fmt)
    if (fmt == "tiff") {
        return(any(candidates$Extension %in% c("tif", "tiff")))
    }
    any(candidates$Extension == fmt)
}
check_format_rule <- function(candidates, required_formats, rule) {
    req <- required_formats_vector(required_formats)
    if (length(req) == 0L) {
        return(TRUE)
    }
    present <- vapply(req, function(fmt) {
        format_present(candidates, fmt)
    }, logical(1))
    if (toupper(rule) == "ANY") {
        any(present)
    }
    else {
        all(present)
    }
}
