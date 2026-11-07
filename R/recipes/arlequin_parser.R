# Original mathematical/parsing helpers; top-level study scripts were not executed.
read_arlequin_text <- function(file) {
    lines <- readLines(file, warn = FALSE, encoding = "latin1")
    list(lines = lines, text = paste(lines, collapse = "\n"))
}
extract_first_number <- function(text, pattern, label) {
    match <- str_match(text, regex(pattern, dotall = TRUE))
    if (nrow(match) == 0 || is.na(match[1, 2])) {
        stop(paste("Could not extract:", label), call. = FALSE)
    }
    as.numeric(match[1, 2])
}
extract_first_integer <- function(text, pattern, label) {
    as.integer(extract_first_number(text, pattern, label))
}
extract_first_string <- function(text, pattern, label) {
    match <- str_match(text, regex(pattern, dotall = TRUE))
    if (nrow(match) == 0 || is.na(match[1, 2])) {
        stop(paste("Could not extract:", label), call. = FALSE)
    }
    str_trim(match[1, 2])
}
extract_xml_block <- function(text, tag) {
    pattern <- paste0("<", tag, "[^>]*>(.*?)</", tag, ">")
    match <- str_match(text, regex(pattern, dotall = TRUE))
    if (nrow(match) == 0 || is.na(match[1, 2])) {
        stop(paste("Could not find XML block:", tag), call. = FALSE)
    }
    match[1, 2]
}
parse_common_settings <- function(arlequin_text, expected_samples) {
    text <- arlequin_text$text
    nb_samples <- extract_first_integer(text, "NbSamples\\s*=\\s*([0-9]+)", "NbSamples")
    data_type <- extract_first_string(text, "DataType\\s*=\\s*([^\\r\\n<]+)", "DataType")
    genotypic_data <- extract_first_integer(text, "GenotypicData\\s*=\\s*([0-9]+)", "GenotypicData")
    requested_permutations <- extract_first_integer(text, "No\\. of\\s+(?:[Pp]ermutations for significance|[Pp]ermutations)\\s*=\\s*([0-9]+)", 
        "requested significance permutations")
    distance_method <- extract_first_string(text, "Molecular distance\\s*:\\s*([^\\r\\n<]+)", "molecular distance")
    usable_loci <- extract_first_integer(text, "Number of usable loci for distance computation\\s*:\\s*([0-9]+)", "usable loci")
    tibble(NbSamples = nb_samples, DataType = data_type, GenotypicData = genotypic_data, Requested_Permutations = requested_permutations, 
        Distance_Method = distance_method, Usable_Loci = usable_loci, Expected_Samples = expected_samples, Sample_Count_Pass = nb_samples == 
            expected_samples, DataType_Pass = data_type == "MICROSAT", Haplotypic_Pass = genotypic_data == 0L, Permutation_Setting_Pass = requested_permutations == 
            REQUESTED_PERMUTATIONS, RST_Distance_Pass = str_detect(distance_method, fixed("Sum of squared size difference (RST)")), 
        Locus_Count_Pass = usable_loci == EXPECTED_LOCI)
}
parse_hierarchical_amova <- function(file) {
    obj <- read_arlequin_text(file)
    text <- obj$text
    common_settings <- parse_common_settings(obj, expected_samples = EXPECTED_ATOLLS)
    df_among_groups <- extract_first_integer(text, "Among\\s+groups\\s+([0-9]+)\\s+[-0-9.]+\\s+[-0-9.]+\\s+Va\\s+[-0-9.]+", 
        "AMOVA among-group df")
    ss_among_groups <- extract_first_number(text, "Among\\s+groups\\s+[0-9]+\\s+([-0-9.]+)\\s+[-0-9.]+\\s+Va\\s+[-0-9.]+", 
        "AMOVA among-group sum of squares")
    var_among_groups <- extract_first_number(text, "Among\\s+groups\\s+[0-9]+\\s+[-0-9.]+\\s+([-0-9.]+)\\s+Va\\s+[-0-9.]+", 
        "AMOVA Va")
    pct_among_groups <- extract_first_number(text, "Among\\s+groups\\s+[0-9]+\\s+[-0-9.]+\\s+[-0-9.]+\\s+Va\\s+([-0-9.]+)", 
        "AMOVA among-group percentage")
    df_among_pop_within <- extract_first_integer(text, "Among\\s+populations\\s+within\\s+groups\\s+([0-9]+)\\s+[-0-9.]+\\s+[-0-9.]+\\s+Vb\\s+[-0-9.]+", 
        "AMOVA among-populations-within-groups df")
    ss_among_pop_within <- extract_first_number(text, "Among\\s+populations\\s+within\\s+groups\\s+[0-9]+\\s+([-0-9.]+)\\s+[-0-9.]+\\s+Vb\\s+[-0-9.]+", 
        "AMOVA among-populations-within-groups sum of squares")
    var_among_pop_within <- extract_first_number(text, "Among\\s+populations\\s+within\\s+groups\\s+[0-9]+\\s+[-0-9.]+\\s+([-0-9.]+)\\s+Vb\\s+[-0-9.]+", 
        "AMOVA Vb")
    pct_among_pop_within <- extract_first_number(text, "Among\\s+populations\\s+within\\s+groups\\s+[0-9]+\\s+[-0-9.]+\\s+[-0-9.]+\\s+Vb\\s+([-0-9.]+)", 
        "AMOVA among-populations-within-groups percentage")
    df_within <- extract_first_integer(text, "Within\\s+populations\\s+([0-9]+)\\s+[-0-9.]+\\s+[-0-9.]+\\s+Vc\\s+[-0-9.]+", 
        "AMOVA within-populations df")
    ss_within <- extract_first_number(text, "Within\\s+populations\\s+[0-9]+\\s+([-0-9.]+)\\s+[-0-9.]+\\s+Vc\\s+[-0-9.]+", 
        "AMOVA within-populations sum of squares")
    var_within <- extract_first_number(text, "Within\\s+populations\\s+[0-9]+\\s+[-0-9.]+\\s+([-0-9.]+)\\s+Vc\\s+[-0-9.]+", 
        "AMOVA Vc")
    pct_within <- extract_first_number(text, "Within\\s+populations\\s+[0-9]+\\s+[-0-9.]+\\s+[-0-9.]+\\s+Vc\\s+([-0-9.]+)", 
        "AMOVA within-populations percentage")
    df_total <- extract_first_integer(text, "Total\\s+([0-9]+)\\s+[-0-9.]+\\s+[-0-9.]+", "AMOVA total df")
    ss_total <- extract_first_number(text, "Total\\s+[0-9]+\\s+([-0-9.]+)\\s+[-0-9.]+", "AMOVA total sum of squares")
    var_total <- extract_first_number(text, "Total\\s+[0-9]+\\s+[-0-9.]+\\s+([-0-9.]+)", "AMOVA total variance")
    fsc <- extract_first_number(text, "FSC\\s*:\\s*([-0-9.]+)", "FSC")
    fst <- extract_first_number(text, "FST\\s*:\\s*([-0-9.]+)", "FST")
    fct <- extract_first_number(text, "FCT\\s*:\\s*([-0-9.]+)", "FCT")
    p_fst <- extract_first_number(text, "Vc and FST.*?P-value\\s*=\\s*([0-9.]+)", "FST P value")
    p_fsc <- extract_first_number(text, "Vb and FSC.*?P-value\\s*=\\s*([0-9.]+)", "FSC P value")
    p_fct <- extract_first_number(text, "Va and FCT.*?P-value\\s*=\\s*([0-9.]+)", "FCT P value")
    reported_permutations <- extract_first_integer(text, "Significance tests \\(([0-9]+) permutations\\)", "AMOVA reported permutations")
    amova_table <- tibble(Source = c("Among regions", "Among administrative divisions within regions", "Within administrative divisions", 
        "Total"), df = c(df_among_groups, df_among_pop_within, df_within, df_total), Sum_of_Squares = c(ss_among_groups, 
        ss_among_pop_within, ss_within, ss_total), Variance_Component = c(var_among_groups, var_among_pop_within, var_within, 
        var_total), Variance_Component_Label = c("Va", "Vb", "Vc", NA_character_), Percentage_of_Variation = c(pct_among_groups, 
        pct_among_pop_within, pct_within, NA_real_), Fixation_Index_Label = c("FCT", "FSC", "FST", NA_character_), Fixation_Index = c(fct, 
        fsc, fst, NA_real_), P_Raw = c(p_fct, p_fsc, p_fst, NA_real_)) %>% mutate(P_Display = case_when(is.na(P_Raw) ~ "", 
        P_Raw == 0 ~ "<0.0001", TRUE ~ format(P_Raw, digits = 5, nsmall = 5, trim = TRUE)))
    list(Settings = common_settings, AMOVA = amova_table, Reported_Permutations = reported_permutations)
}
parse_population_labels <- function(text) {
    block <- extract_xml_block(text, "pairDistPopLabels")
    lines <- str_split(block, "\n", simplify = FALSE)[[1]]
    matches <- str_match(lines, "^\\s*([0-9]+):\\s*([A-Za-z0-9_-]+)")
    labels <- tibble(Index = suppressWarnings(as.integer(matches[, 2])), Arlequin_Name = matches[, 3]) %>% filter(!is.na(Index), 
        !is.na(Arlequin_Name)) %>% arrange(Index)
    labels
}
parse_lower_triangle_numeric <- function(block, n_populations, p_value_mode = FALSE) {
    lines <- str_split(block, "\n", simplify = FALSE)[[1]]
    output <- matrix(NA_real_, nrow = n_populations, ncol = n_populations)
    for (line in lines) {
        tokens <- str_split(str_trim(line), "\\s+", simplify = FALSE)[[1]]
        if (length(tokens) == 0 || !str_detect(tokens[1], "^[0-9]+$")) {
            next
        }
        row_index <- as.integer(tokens[1])
        if (row_index < 1L || row_index > n_populations) {
            next
        }
        values <- tokens[-1]
        if (p_value_mode) {
            parsed_values <- vapply(values, function(x) {
                if (x == "*") {
                  return(NA_real_)
                }
                as.numeric(str_split(x, fixed("+-"), simplify = TRUE)[1])
            }, numeric(1))
        }
        else {
            parsed_values <- suppressWarnings(as.numeric(values))
        }
        expected_value_count <- row_index
        if (length(parsed_values) != expected_value_count) {
            next
        }
        for (column_index in seq_len(row_index)) {
            value <- parsed_values[column_index]
            output[row_index, column_index] <- value
            output[column_index, row_index] <- value
        }
    }
    output
}
parse_pairwise_arlequin <- function(file, expected_populations) {
    obj <- read_arlequin_text(file)
    text <- obj$text
    settings <- parse_common_settings(obj, expected_samples = expected_populations)
    labels <- parse_population_labels(text)
    if (nrow(labels) != expected_populations) {
        stop(paste("Expected", expected_populations, "population labels but parsed", nrow(labels)), call. = FALSE)
    }
    rst_block <- extract_xml_block(text, "PairFstMat")
    p_block <- extract_xml_block(text, "PairFstPvalMat")
    rst_matrix <- parse_lower_triangle_numeric(rst_block, n_populations = expected_populations, p_value_mode = FALSE)
    p_matrix <- parse_lower_triangle_numeric(p_block, n_populations = expected_populations, p_value_mode = TRUE)
    if (any(is.na(rst_matrix))) {
        stop("Pairwise RST matrix contains unparsed cells.", call. = FALSE)
    }
    if (any(is.na(p_matrix[lower.tri(p_matrix)]))) {
        stop("Pairwise P-value matrix contains unparsed off-diagonal cells.", call. = FALSE)
    }
    population_names <- labels$Arlequin_Name
    rownames(rst_matrix) <- population_names
    colnames(rst_matrix) <- population_names
    rownames(p_matrix) <- population_names
    colnames(p_matrix) <- population_names
    if (any(abs(diag(rst_matrix)) > 1e-12)) {
        stop("Pairwise RST matrix diagonal is not zero.", call. = FALSE)
    }
    if (max(abs(rst_matrix - t(rst_matrix)), na.rm = TRUE) > 1e-12) {
        stop("Pairwise RST matrix is not symmetric.", call. = FALSE)
    }
    pair_rows <- vector("list", choose(expected_populations, 2))
    counter <- 1L
    for (i in seq_len(expected_populations)) {
        if (i == 1L) {
            next
        }
        for (j in seq_len(i - 1L)) {
            pair_rows[[counter]] <- tibble(Index_1 = j, Index_2 = i, Population_1 = population_names[j], Population_2 = population_names[i], 
                RST = rst_matrix[i, j], P_Raw = p_matrix[i, j])
            counter <- counter + 1L
        }
    }
    pairwise_long <- bind_rows(pair_rows)
    if (any(is.na(pairwise_long$P_Raw))) {
        stop("At least one off-diagonal pairwise P value could not be parsed.", call. = FALSE)
    }
    pairwise_long <- pairwise_long %>% mutate(P_Bonferroni = p.adjust(P_Raw, method = "bonferroni"), P_BH = p.adjust(P_Raw, 
        method = "BH"), Significant_Raw = P_Raw < ALPHA, Significant_Bonferroni = P_Bonferroni < ALPHA, Significant_BH = P_BH < 
        ALPHA)
    reported_permutation_matches <- str_match_all(text, "Number of permutations\\s*:\\s*([0-9]+)")[[1]]
    reported_permutations <- unique(suppressWarnings(as.integer(reported_permutation_matches[, 2])))
    reported_permutations <- reported_permutations[!is.na(reported_permutations)]
    list(Settings = settings, Labels = labels, RST_Matrix = rst_matrix, P_Raw_Matrix = p_matrix, Pairwise_Long = pairwise_long, 
        Reported_Permutations = reported_permutations)
}
pairwise_to_matrix <- function(pairwise_data, labels, label_1_column, label_2_column, value_column, diagonal_value = NA_real_) {
    mat <- matrix(NA_real_, nrow = length(labels), ncol = length(labels), dimnames = list(labels, labels))
    diag(mat) <- diagonal_value
    for (i in seq_len(nrow(pairwise_data))) {
        label_1 <- pairwise_data[[label_1_column]][i]
        label_2 <- pairwise_data[[label_2_column]][i]
        value <- pairwise_data[[value_column]][i]
        mat[label_1, label_2] <- value
        mat[label_2, label_1] <- value
    }
    mat
}
