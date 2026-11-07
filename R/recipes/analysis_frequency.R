# Adapted from the user's original analysis_frequency final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
validate_level <- function(level) {
    allowed <- c("overall", "region", "atoll")
    if (length(level) != 1 || !level %in% allowed) {
        stop(paste0("level must be one of: ", paste(allowed, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
validate_required_columns <- function(data, required_columns) {
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0) {
        stop(paste("Missing required column(s):", paste(missing_columns, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
validate_unique_sample_ids <- function(data) {
    validate_required_columns(data, "sample_id")
    if (any(is.na(data$sample_id) | data$sample_id == "")) {
        stop("Blank Sample_ID values detected.", call. = FALSE)
    }
    if (n_distinct(data$sample_id) != nrow(data)) {
        stop("Sample_ID values are not unique.", call. = FALSE)
    }
    invisible(TRUE)
}
is_blank_value <- function(x) {
    is.na(x) | str_trim(as.character(x)) == ""
}
is_deletion_value <- function(x) {
    !is.na(x) & toupper(str_trim(as.character(x))) == "DEL"
}
is_successfully_typed_single_copy <- function(x) {
    !is_blank_value(x) & !is_deletion_value(x)
}
is_numeric_ystr_allele <- function(x) {
    x_chr <- str_trim(as.character(x))
    is_successfully_typed_single_copy(x_chr) & str_detect(x_chr, "^[0-9]+(?:\\.[0-9]+)?$")
}
add_population_labels <- function(data, level) {
    validate_level(level)
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
get_grouping_vars <- function(level) {
    validate_level(level)
    switch(level, overall = character(0), region = "region", atoll = c("region", "atoll"))
}
prepare_single_copy_long <- function(data, loci = SINGLE_COPY_LOCI) {
    validate_unique_sample_ids(data)
    validate_required_columns(data, c("sample_id", "region", "atoll", loci))
    long <- data %>% select(sample_id, region, atoll, all_of(loci)) %>% pivot_longer(cols = all_of(loci), names_to = "locus_internal", 
        values_to = "allele_raw") %>% mutate(allele_raw = str_trim(as.character(allele_raw)), Locus = unname(LOCUS_DISPLAY_NAMES[locus_internal]), 
        blank = is_blank_value(allele_raw), structural_deletion = is_deletion_value(allele_raw), successfully_typed = is_successfully_typed_single_copy(allele_raw), 
        valid_numeric_allele = is_numeric_ystr_allele(allele_raw))
    invalid_typed <- long %>% filter(successfully_typed, !valid_numeric_allele)
    if (nrow(invalid_typed) > 0) {
        stop(paste0("Found ", nrow(invalid_typed), " successfully typed single-copy value(s) that are not valid ", "integer or microvariant Y-STR allele strings."), 
            call. = FALSE)
    }
    long
}
calculate_single_copy_frequencies <- function(data, level = c("overall", "region", "atoll"), loci = SINGLE_COPY_LOCI) {
    level <- match.arg(level)
    grouping_vars <- get_grouping_vars(level)
    long <- prepare_single_copy_long(data, loci = loci)
    denominators <- long %>% group_by(across(all_of(c(grouping_vars, "locus_internal", "Locus")))) %>% summarise(Successfully_Typed_N = sum(successfully_typed), 
        .groups = "drop")
    counts <- long %>% filter(successfully_typed) %>% group_by(across(all_of(c(grouping_vars, "locus_internal", "Locus"))), 
        Allele = allele_raw) %>% summarise(Allele_Count = n(), .groups = "drop")
    result <- counts %>% left_join(denominators, by = c(grouping_vars, "locus_internal", "Locus")) %>% mutate(Allele_Frequency = Allele_Count/Successfully_Typed_N, 
        allele_sort_numeric = suppressWarnings(as.numeric(Allele)), locus_order = match(locus_internal, loci))
    result <- add_population_labels(result, level)
    result %>% arrange(Population, Region, Administrative_Division, locus_order, allele_sort_numeric, Allele) %>% select(Population, 
        Region, Administrative_Division, Locus, Allele, Allele_Count, Successfully_Typed_N, Allele_Frequency)
}
calculate_single_copy_typing_summary <- function(data, level = c("overall", "region", "atoll"), loci = SINGLE_COPY_LOCI) {
    level <- match.arg(level)
    grouping_vars <- get_grouping_vars(level)
    long <- prepare_single_copy_long(data, loci = loci)
    result <- long %>% group_by(across(all_of(c(grouping_vars, "locus_internal", "Locus")))) %>% summarise(Total_Profiles_N = n(), 
        Successfully_Typed_N = sum(successfully_typed), Blank_N = sum(blank), Structural_Deletion_N = sum(structural_deletion), 
        Not_Successfully_Typed_N = sum(!successfully_typed), .groups = "drop") %>% mutate(Typing_Rate = Successfully_Typed_N/Total_Profiles_N, 
        Count_Check = Total_Profiles_N - Successfully_Typed_N - Not_Successfully_Typed_N, locus_order = match(locus_internal, 
            loci))
    result <- add_population_labels(result, level)
    result %>% arrange(Population, Region, Administrative_Division, locus_order) %>% select(Population, Region, Administrative_Division, 
        Locus, Total_Profiles_N, Successfully_Typed_N, Blank_N, Structural_Deletion_N, Not_Successfully_Typed_N, Typing_Rate, 
        Count_Check)
}
prepare_multicopy_long <- function(data, state_columns = MULTICOPY_STATE_COLUMNS) {
    validate_unique_sample_ids(data)
    validate_required_columns(data, c("sample_id", "region", "atoll", names(state_columns)))
    long <- data %>% select(sample_id, region, atoll, all_of(names(state_columns))) %>% pivot_longer(cols = all_of(names(state_columns)), 
        names_to = "state_column", values_to = "allele_combination_raw") %>% mutate(allele_combination_raw = str_trim(as.character(allele_combination_raw)), 
        Marker = unname(state_columns[state_column]), successfully_typed = !is_blank_value(allele_combination_raw))
    long
}
multicopy_first_allele_numeric <- function(x) {
    vapply(strsplit(as.character(x), "-", fixed = TRUE), function(parts) {
        values <- suppressWarnings(as.numeric(parts))
        if (length(values) == 0 || all(is.na(values))) {
            return(Inf)
        }
        min(values, na.rm = TRUE)
    }, numeric(1))
}
multicopy_copy_count <- function(x) {
    vapply(strsplit(as.character(x), "-", fixed = TRUE), length, integer(1))
}
calculate_multicopy_frequencies <- function(data, level = c("overall", "region", "atoll"), state_columns = MULTICOPY_STATE_COLUMNS) {
    level <- match.arg(level)
    grouping_vars <- get_grouping_vars(level)
    long <- prepare_multicopy_long(data, state_columns = state_columns)
    denominators <- long %>% group_by(across(all_of(c(grouping_vars, "state_column", "Marker")))) %>% summarise(Successfully_Typed_N = sum(successfully_typed), 
        .groups = "drop")
    counts <- long %>% filter(successfully_typed) %>% group_by(across(all_of(c(grouping_vars, "state_column", "Marker"))), 
        Allele_Combination = allele_combination_raw) %>% summarise(Combination_Count = n(), .groups = "drop")
    result <- counts %>% left_join(denominators, by = c(grouping_vars, "state_column", "Marker")) %>% mutate(Combination_Frequency = Combination_Count/Successfully_Typed_N, 
        marker_order = match(state_column, names(state_columns)), first_allele_sort = multicopy_first_allele_numeric(Allele_Combination), 
        observed_copy_count = multicopy_copy_count(Allele_Combination))
    result <- add_population_labels(result, level)
    result %>% arrange(Population, Region, Administrative_Division, marker_order, first_allele_sort, observed_copy_count, 
        Allele_Combination) %>% select(Population, Region, Administrative_Division, Marker, Allele_Combination, Combination_Count, 
        Successfully_Typed_N, Combination_Frequency)
}
calculate_multicopy_typing_summary <- function(data, level = c("overall", "region", "atoll"), state_columns = MULTICOPY_STATE_COLUMNS) {
    level <- match.arg(level)
    grouping_vars <- get_grouping_vars(level)
    long <- prepare_multicopy_long(data, state_columns = state_columns)
    result <- long %>% group_by(across(all_of(c(grouping_vars, "state_column", "Marker")))) %>% summarise(Total_Profiles_N = n(), 
        Successfully_Typed_N = sum(successfully_typed), Not_Successfully_Typed_N = sum(!successfully_typed), .groups = "drop") %>% 
        mutate(Typing_Rate = Successfully_Typed_N/Total_Profiles_N, Count_Check = Total_Profiles_N - Successfully_Typed_N - 
            Not_Successfully_Typed_N, marker_order = match(state_column, names(state_columns)))
    result <- add_population_labels(result, level)
    result %>% arrange(Population, Region, Administrative_Division, marker_order) %>% select(Population, Region, Administrative_Division, 
        Marker, Total_Profiles_N, Successfully_Typed_N, Not_Successfully_Typed_N, Typing_Rate, Count_Check)
}
validate_single_copy_frequency_table <- function(frequency_table, tolerance = 1e-12) {
    required <- c("Population", "Region", "Administrative_Division", "Locus", "Allele", "Allele_Count", "Successfully_Typed_N", 
        "Allele_Frequency")
    validate_required_columns(frequency_table, required)
    check <- frequency_table %>% group_by(Population, Region, Administrative_Division, Locus) %>% summarise(Count_Sum = sum(Allele_Count), 
        Denominator = first(Successfully_Typed_N), Frequency_Sum = sum(Allele_Frequency), Denominator_Values = n_distinct(Successfully_Typed_N), 
        .groups = "drop") %>% mutate(Count_Check = Count_Sum == Denominator, Frequency_Check = abs(Frequency_Sum - 1) <= 
        tolerance, Denominator_Consistency = Denominator_Values == 1)
    if (any(!check$Count_Check | !check$Frequency_Check | !check$Denominator_Consistency)) {
        stop("Single-copy allele-frequency validation failed.", call. = FALSE)
    }
    invisible(check)
}
validate_multicopy_frequency_table <- function(frequency_table, tolerance = 1e-12) {
    required <- c("Population", "Region", "Administrative_Division", "Marker", "Allele_Combination", "Combination_Count", 
        "Successfully_Typed_N", "Combination_Frequency")
    validate_required_columns(frequency_table, required)
    check <- frequency_table %>% group_by(Population, Region, Administrative_Division, Marker) %>% summarise(Count_Sum = sum(Combination_Count), 
        Denominator = first(Successfully_Typed_N), Frequency_Sum = sum(Combination_Frequency), Denominator_Values = n_distinct(Successfully_Typed_N), 
        .groups = "drop") %>% mutate(Count_Check = Count_Sum == Denominator, Frequency_Check = abs(Frequency_Sum - 1) <= 
        tolerance, Denominator_Consistency = Denominator_Values == 1)
    if (any(!check$Count_Check | !check$Frequency_Check | !check$Denominator_Consistency)) {
        stop("Multicopy frequency validation failed.", call. = FALSE)
    }
    invisible(check)
}
calculate_all_allele_frequency_tables <- function(data, single_copy_loci = SINGLE_COPY_LOCI, multicopy_state_columns = MULTICOPY_STATE_COLUMNS) {
    output <- list(Overall_SingleCopy = calculate_single_copy_frequencies(data, level = "overall", loci = single_copy_loci), 
        Overall_Multicopy = calculate_multicopy_frequencies(data, level = "overall", state_columns = multicopy_state_columns), 
        Regions_SingleCopy = calculate_single_copy_frequencies(data, level = "region", loci = single_copy_loci), Regions_Multicopy = calculate_multicopy_frequencies(data, 
            level = "region", state_columns = multicopy_state_columns), Atolls_SingleCopy = calculate_single_copy_frequencies(data, 
            level = "atoll", loci = single_copy_loci), Atolls_Multicopy = calculate_multicopy_frequencies(data, level = "atoll", 
            state_columns = multicopy_state_columns), Overall_SingleCopy_Typing = calculate_single_copy_typing_summary(data, 
            level = "overall", loci = single_copy_loci), Regions_SingleCopy_Typing = calculate_single_copy_typing_summary(data, 
            level = "region", loci = single_copy_loci), Atolls_SingleCopy_Typing = calculate_single_copy_typing_summary(data, 
            level = "atoll", loci = single_copy_loci), Overall_Multicopy_Typing = calculate_multicopy_typing_summary(data, 
            level = "overall", state_columns = multicopy_state_columns), Regions_Multicopy_Typing = calculate_multicopy_typing_summary(data, 
            level = "region", state_columns = multicopy_state_columns), Atolls_Multicopy_Typing = calculate_multicopy_typing_summary(data, 
            level = "atoll", state_columns = multicopy_state_columns))
    validate_single_copy_frequency_table(output$Overall_SingleCopy)
    validate_single_copy_frequency_table(output$Regions_SingleCopy)
    validate_single_copy_frequency_table(output$Atolls_SingleCopy)
    validate_multicopy_frequency_table(output$Overall_Multicopy)
    validate_multicopy_frequency_table(output$Regions_Multicopy)
    validate_multicopy_frequency_table(output$Atolls_Multicopy)
    output
}
single_copy_frequency_wide <- function(long_frequency_table) {
    validate_single_copy_frequency_table(long_frequency_table)
    long_frequency_table %>% select(Population, Region, Administrative_Division, Locus, Allele, Allele_Frequency) %>% pivot_wider(names_from = Allele, 
        values_from = Allele_Frequency)
}
multicopy_frequency_wide <- function(long_frequency_table) {
    validate_multicopy_frequency_table(long_frequency_table)
    long_frequency_table %>% select(Population, Region, Administrative_Division, Marker, Allele_Combination, Combination_Frequency) %>% 
        pivot_wider(names_from = Allele_Combination, values_from = Combination_Frequency)
}
