# Original mathematical/parsing helpers; top-level study scripts were not executed.
make_safe_token <- function(x) {
    out <- as.character(x)
    out <- str_replace_all(out, "[^A-Za-z0-9_-]", "_")
    out <- str_replace_all(out, "_+", "_")
    out <- str_replace_all(out, "^_+|_+$", "")
    if_else(out == "", "UNNAMED", out)
}
make_sample_name <- function(prefix, x) {
    paste0(prefix, "_", make_safe_token(x))
}
format_arlequin_integer <- function(x) {
    if (any(is.na(x))) {
        stop("Missing value encountered while writing Arlequin allele state.", call. = FALSE)
    }
    as.character(as.integer(round(x)))
}
build_sample_block <- function(sample_data, sample_name) {
    if (nrow(sample_data) == 0) {
        stop(paste("Cannot create empty Arlequin sample:", sample_name), call. = FALSE)
    }
    header <- c(paste0("SampleName=\"", sample_name, "\""), paste0("SampleSize=", nrow(sample_data)), "SampleData= {")
    sample_rows <- vapply(seq_len(nrow(sample_data)), function(i) {
        sample_id_safe <- make_safe_token(sample_data$sample_id[i])
        allele_values <- vapply(STRUCTURE_LOCI, function(locus) {
            format_arlequin_integer(sample_data[[locus]][i])
        }, character(1))
        paste(c(sample_id_safe, "1", allele_values), collapse = "\t")
    }, character(1))
    c(header, sample_rows, "}", "")
}
build_profile_header <- function(title, nb_samples) {
    c("# ================================================================", "# Arlequin 3.5 project file", "# Configured Y-STR population-structure analysis", 
        "# Complete-profile dataset", "# 23 single-copy Y-STR repeat states", "# DYS389II represented as DYS389IIadj = DYS389II - DYS389I", 
        "# Integer repeat states are used only for the Arlequin RST/AMOVA branch", "# ================================================================", 
        "", "[Profile]", paste0("Title=\"", title, "\""), paste0("NbSamples=", nb_samples), "DataType=MICROSAT", "GenotypicData=0", 
        "GameticPhase=0", "RecessiveData=0", "LocusSeparator=WHITESPACE", "MissingData='?'", "Frequency=ABS", "", "[Data]", 
        "[[Samples]]", "", paste0("# Number of loci: ", length(STRUCTURE_LOCI)), paste0("# Locus order: ", paste(unname(STRUCTURE_LOCI_DISPLAY[STRUCTURE_LOCI]), 
            collapse = " ")), "")
}
validate_common_arp <- function(file, expected_samples, expected_profiles) {
    lines <- readLines(file, warn = FALSE)
    checks <- tibble(Check = c("Profile section present", "Data section present", "Samples section present", "DataType MICROSAT", 
        "GenotypicData 0", "Locus separator whitespace", "NbSamples correct", "SampleName block count correct", "SampleSize total correct", 
        "No decimal point in sample data lines"), Pass = c(any(lines == "[Profile]"), any(lines == "[Data]"), any(lines == 
        "[[Samples]]"), any(lines == "DataType=MICROSAT"), any(lines == "GenotypicData=0"), any(lines == "LocusSeparator=WHITESPACE"), 
        any(lines == paste0("NbSamples=", expected_samples)), sum(str_detect(lines, "^SampleName=\"")) == expected_samples, 
        sum(as.integer(str_extract(lines[str_detect(lines, "^SampleSize=")], "[0-9]+$"))) == expected_profiles, !any(str_detect(lines[str_detect(lines, 
            "^[A-Za-z0-9_-]+\\t1\\t")], "\\."))))
    if (any(!checks$Pass)) {
        stop(paste("Arlequin project validation failed for:", basename(file), "| Failed check(s):", paste(checks$Check[!checks$Pass], 
            collapse = "; ")), call. = FALSE)
    }
    checks
}
