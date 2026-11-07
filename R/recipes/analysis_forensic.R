# Adapted from the user's original analysis_forensic final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
validate_forensic_level <- function(level) {
    allowed <- c("overall", "region", "atoll")
    if (length(level) != 1 || !level %in% allowed) {
        stop(paste0("level must be one of: ", paste(allowed, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
validate_forensic_required_columns <- function(data, required_columns) {
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0) {
        stop(paste("Missing required column(s):", paste(missing_columns, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
validate_forensic_unique_ids <- function(data) {
    validate_forensic_required_columns(data, "sample_id")
    if (any(is.na(data$sample_id) | str_trim(as.character(data$sample_id)) == "")) {
        stop("Blank Sample_ID values detected.", call. = FALSE)
    }
    if (n_distinct(data$sample_id) != nrow(data)) {
        stop("Sample_ID values are not unique.", call. = FALSE)
    }
    invisible(TRUE)
}
get_forensic_grouping_vars <- function(level) {
    validate_forensic_level(level)
    switch(level, overall = character(0), region = "region", atoll = c("region", "atoll"))
}
add_forensic_population_labels <- function(data, level) {
    validate_forensic_level(level)
    if (level == "overall") {
        data %>% mutate(Population = study$population_label, Region = "All", Administrative_Division = "All")
    }
    else if (level == "region") {
        data %>% mutate(Population = study$population_label, Region = as.character(region), Administrative_Division = "All")
    }
    else {
        data %>% mutate(Population = study$population_label, Region = as.character(region), Administrative_Division = as.character(atoll))
    }
}
validate_forensic_marker_completeness <- function(data, markers) {
    validate_forensic_unique_ids(data)
    validate_forensic_required_columns(data, c("sample_id", "region", "atoll", markers))
    blank_matrix <- sapply(data[markers], function(x) {
        is.na(x) | str_trim(as.character(x)) == ""
    })
    if (is.null(dim(blank_matrix))) {
        blank_matrix <- matrix(blank_matrix, ncol = 1)
    }
    blank_n_by_profile <- rowSums(blank_matrix)
    if (any(blank_n_by_profile > 0)) {
        affected <- data$sample_id[blank_n_by_profile > 0]
        stop(paste0("Forensic haplotype construction requires complete observed states ", "for every selected marker system. Blank state(s) detected in ", 
            length(affected), " profile(s). Example Sample_ID(s): ", paste(head(affected, 10), collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
marker_set_signature <- function(markers) {
    validate_forensic_required_columns(tibble(dummy = 1), "dummy")
    if (length(markers) == 0) {
        stop("At least one marker system must be selected.", call. = FALSE)
    }
    paste(markers, collapse = ";")
}
construct_haplotype_strings <- function(data, markers, include_marker_names = TRUE) {
    validate_forensic_marker_completeness(data, markers)
    marker_values <- data %>% select(all_of(markers)) %>% mutate(across(everything(), ~str_trim(as.character(.x))))
    if (include_marker_names) {
        haplotype_strings <- pmap_chr(marker_values, function(...) {
            vals <- c(...)
            paste(paste0(markers, "=", vals), collapse = "|")
        })
    }
    else {
        haplotype_strings <- pmap_chr(marker_values, function(...) {
            vals <- c(...)
            paste(vals, collapse = "|")
        })
    }
    data %>% mutate(Haplotype = haplotype_strings)
}
calculate_haplotype_frequencies <- function(data, markers, level = c("overall", "region", "atoll"), panel_name = "Full Yfiler Plus") {
    level <- match.arg(level)
    grouping_vars <- get_forensic_grouping_vars(level)
    hap_data <- construct_haplotype_strings(data, markers = markers, include_marker_names = TRUE)
    counts <- hap_data %>% group_by(across(all_of(c(grouping_vars, "Haplotype")))) %>% summarise(Haplotype_Count = n(), .groups = "drop")
    if (length(grouping_vars) == 0) {
        denominators <- tibble(Population_N = nrow(hap_data))
        counts <- counts %>% mutate(Population_N = denominators$Population_N)
    }
    else {
        denominators <- hap_data %>% group_by(across(all_of(grouping_vars))) %>% summarise(Population_N = n(), .groups = "drop")
        counts <- counts %>% left_join(denominators, by = grouping_vars)
    }
    result <- counts %>% mutate(Haplotype_Frequency = Haplotype_Count/Population_N, Singleton = Haplotype_Count == 1)
    result <- add_forensic_population_labels(result, level)
    result %>% mutate(Panel = panel_name, Marker_System_Count = length(markers), Marker_Set = marker_set_signature(markers)) %>% 
        arrange(Population, Region, Administrative_Division, desc(Haplotype_Count), Haplotype) %>% select(Population, Region, 
        Administrative_Division, Panel, Marker_System_Count, Marker_Set, Haplotype, Haplotype_Count, Population_N, Haplotype_Frequency, 
        Singleton)
}
calculate_forensic_parameters_from_frequencies <- function(haplotype_frequency_table) {
    required <- c("Population", "Region", "Administrative_Division", "Panel", "Marker_System_Count", "Marker_Set", "Haplotype", 
        "Haplotype_Count", "Population_N", "Haplotype_Frequency", "Singleton")
    validate_forensic_required_columns(haplotype_frequency_table, required)
    result <- haplotype_frequency_table %>% group_by(Population, Region, Administrative_Division, Panel, Marker_System_Count, 
        Marker_Set) %>% summarise(N = first(Population_N), Distinct_Haplotypes = n(), Unique_Haplotypes = sum(Singleton), 
        Nonunique_Distinct_Haplotypes = sum(!Singleton), Profiles_in_Repeated_Haplotypes = sum(Haplotype_Count[Haplotype_Count > 
            1]), HMP = sum(Haplotype_Frequency^2), HD = if_else(first(Population_N) > 1, first(Population_N)/(first(Population_N) - 
            1) * (1 - sum(Haplotype_Frequency^2)), NA_real_), DC = n()/first(Population_N), Unique_Profile_Fraction = sum(Singleton)/first(Population_N), 
        Unique_Haplotype_Fraction = sum(Singleton)/n(), Most_Common_Haplotype_Count = max(Haplotype_Count), Most_Common_Haplotype_Frequency = max(Haplotype_Frequency), 
        .groups = "drop")
    result
}
calculate_forensic_parameters <- function(data, markers, level = c("overall", "region", "atoll"), panel_name = "Full Yfiler Plus") {
    level <- match.arg(level)
    frequency_table <- calculate_haplotype_frequencies(data = data, markers = markers, level = level, panel_name = panel_name)
    calculate_forensic_parameters_from_frequencies(frequency_table)
}
calculate_forensic_all_levels <- function(data, markers, panel_name = "Full Yfiler Plus") {
    list(Overall_Haplotype_Frequencies = calculate_haplotype_frequencies(data = data, markers = markers, level = "overall", 
        panel_name = panel_name), Region_Haplotype_Frequencies = calculate_haplotype_frequencies(data = data, markers = markers, 
        level = "region", panel_name = panel_name), Atoll_Haplotype_Frequencies = calculate_haplotype_frequencies(data = data, 
        markers = markers, level = "atoll", panel_name = panel_name), Overall_Forensic_Parameters = calculate_forensic_parameters(data = data, 
        markers = markers, level = "overall", panel_name = panel_name), Region_Forensic_Parameters = calculate_forensic_parameters(data = data, 
        markers = markers, level = "region", panel_name = panel_name), Atoll_Forensic_Parameters = calculate_forensic_parameters(data = data, 
        markers = markers, level = "atoll", panel_name = panel_name))
}
validate_haplotype_frequency_table <- function(frequency_table, tolerance = 1e-12) {
    required <- c("Population", "Region", "Administrative_Division", "Panel", "Haplotype", "Haplotype_Count", "Population_N", 
        "Haplotype_Frequency")
    validate_forensic_required_columns(frequency_table, required)
    check <- frequency_table %>% group_by(Population, Region, Administrative_Division, Panel) %>% summarise(Count_Sum = sum(Haplotype_Count), 
        Population_N = first(Population_N), Population_N_Values = n_distinct(Population_N), Frequency_Sum = sum(Haplotype_Frequency), 
        .groups = "drop") %>% mutate(Count_Check = Count_Sum == Population_N, Population_N_Consistency = Population_N_Values == 
        1, Frequency_Check = abs(Frequency_Sum - 1) <= tolerance)
    if (any(!check$Count_Check | !check$Population_N_Consistency | !check$Frequency_Check)) {
        stop("Haplotype-frequency validation failed.", call. = FALSE)
    }
    invisible(check)
}
validate_forensic_parameter_table <- function(parameter_table, tolerance = 1e-12) {
    required <- c("N", "Distinct_Haplotypes", "Unique_Haplotypes", "HMP", "HD", "DC", "Unique_Profile_Fraction", "Unique_Haplotype_Fraction", 
        "Most_Common_Haplotype_Count", "Most_Common_Haplotype_Frequency")
    validate_forensic_required_columns(parameter_table, required)
    invalid <- parameter_table %>% mutate(distinct_valid = Distinct_Haplotypes >= 1 & Distinct_Haplotypes <= N, unique_valid = Unique_Haplotypes >= 
        0 & Unique_Haplotypes <= Distinct_Haplotypes, hmp_valid = HMP >= -tolerance & HMP <= 1 + tolerance, hd_valid = is.na(HD) | 
        (HD >= -tolerance & HD <= 1 + tolerance), dc_valid = DC >= -tolerance & DC <= 1 + tolerance, unique_profile_valid = Unique_Profile_Fraction >= 
        -tolerance & Unique_Profile_Fraction <= 1 + tolerance, unique_haplotype_valid = Unique_Haplotype_Fraction >= -tolerance & 
        Unique_Haplotype_Fraction <= 1 + tolerance, max_count_valid = Most_Common_Haplotype_Count >= 1 & Most_Common_Haplotype_Count <= 
        N, max_frequency_valid = Most_Common_Haplotype_Frequency >= -tolerance & Most_Common_Haplotype_Frequency <= 1 + tolerance, 
        dc_formula_valid = abs(DC - Distinct_Haplotypes/N) <= tolerance, unique_profile_formula_valid = abs(Unique_Profile_Fraction - 
            Unique_Haplotypes/N) <= tolerance, unique_haplotype_formula_valid = abs(Unique_Haplotype_Fraction - Unique_Haplotypes/Distinct_Haplotypes) <= 
            tolerance, max_frequency_formula_valid = abs(Most_Common_Haplotype_Frequency - Most_Common_Haplotype_Count/N) <= 
            tolerance) %>% filter(!distinct_valid | !unique_valid | !hmp_valid | !hd_valid | !dc_valid | !unique_profile_valid | 
        !unique_haplotype_valid | !max_count_valid | !max_frequency_valid | !dc_formula_valid | !unique_profile_formula_valid | 
        !unique_haplotype_formula_valid | !max_frequency_formula_valid)
    if (nrow(invalid) > 0) {
        stop("Forensic parameter validation failed.", call. = FALSE)
    }
    invisible(TRUE)
}
remove_one_marker <- function(markers, marker_to_remove) {
    if (!marker_to_remove %in% markers) {
        stop(paste(marker_to_remove, "is not present in the supplied marker set."), call. = FALSE)
    }
    remaining <- setdiff(markers, marker_to_remove)
    if (length(remaining) == 0) {
        stop("Marker removal would leave an empty marker set.", call. = FALSE)
    }
    remaining
}
calculate_single_marker_forensics <- function(data, markers, level = "overall") {
    map_dfr(markers, function(marker) {
        display_name <- if (marker %in% names(FORENSIC_MARKER_DISPLAY_NAMES)) {
            unname(FORENSIC_MARKER_DISPLAY_NAMES[marker])
        }
        else {
            marker
        }
        calculate_forensic_parameters(data = data, markers = marker, level = level, panel_name = paste0("Single marker: ", 
            display_name)) %>% mutate(Marker = display_name, Marker_Internal = marker)
    })
}
format_forensic_parameters <- function(parameter_table, digits = 6) {
    parameter_table %>% mutate(HMP = round(HMP, digits), HD = round(HD, digits), DC = round(DC, digits), Unique_Profile_Fraction = round(Unique_Profile_Fraction, 
        digits), Unique_Haplotype_Fraction = round(Unique_Haplotype_Fraction, digits), Most_Common_Haplotype_Frequency = round(Most_Common_Haplotype_Frequency, 
        digits))
}
