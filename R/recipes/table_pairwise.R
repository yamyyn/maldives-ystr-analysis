# Adapted from the user's original table_pairwise final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
matrix_to_export_df <- function(mat, row_name_column = "Population") {
    out <- as.data.frame(mat, check.names = FALSE)
    out <- tibble::rownames_to_column(out, var = row_name_column)
    out
}
write_workbook <- function(path) {
    openxlsx::write.xlsx(list(AMOVA = amova_table, AMOVA_Manuscript = amova_manuscript, Atoll_Pairs = atoll_pairs, Atoll_RST_Matrix = matrix_to_export_df(atoll_rst_matrix, 
        "Administrative_Division"), Atoll_P_Raw = matrix_to_export_df(atoll_p_raw_matrix, "Administrative_Division"), Atoll_P_Bonf = matrix_to_export_df(atoll_p_bonf_matrix, 
        "Administrative_Division"), Atoll_P_BH = matrix_to_export_df(atoll_p_bh_matrix, "Administrative_Division"), Region_Pairs = region_pairs, 
        Region_RST_Matrix = matrix_to_export_df(region_rst_matrix, "Region"), Region_P_Raw = matrix_to_export_df(region_p_raw_matrix, 
            "Region"), Region_P_Bonf = matrix_to_export_df(region_p_bonf_matrix, "Region"), Region_P_BH = matrix_to_export_df(region_p_bh_matrix, 
            "Region"), Significance_Summary = pairwise_significance_summary, Negative_RST_Audit = negative_rst_audit, Settings = all_settings, 
        Permutation_Documentation = permutation_documentation, Source_File_Index = source_file_index, Figure_Captions = figure_captions), 
        path, overwrite = TRUE)
}
