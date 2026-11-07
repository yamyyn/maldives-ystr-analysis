# Adapted from the user's original table_catalogue final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_workbook <- function(path) {
    {
        wb <- createWorkbook(creator = study$population_label)
        header_style <- createStyle(textDecoration = "bold", halign = "center", valign = "center", border = "Bottom", borderStyle = "thin", 
            fgFill = "#EDEDED")
        freq_style <- createStyle(numFmt = "0.000000")
        hmp_style <- createStyle(numFmt = "0.0000000000")
        add_sheet <- function(wb, sheet, x) {
            addWorksheet(wb, sheet)
            writeData(wb, sheet, x, headerStyle = header_style, withFilter = TRUE)
            freezePane(wb, sheet, firstActiveRow = 2)
            setColWidths(wb, sheet, 1:ncol(x), "auto")
        }
        add_sheet(wb, "Distinct_Haplotypes", catalogue)
        add_sheet(wb, "Recurrent_Haplotypes", recurrent_catalogue)
        add_sheet(wb, "Haplotype_Admin_Distribution", admin_distribution)
        add_sheet(wb, "Haplotype_Region_Distribution", region_distribution)
        add_sheet(wb, "Predictor_Consistency", predictor_consistency)
        add_sheet(wb, "Catalogue_Summary", summary_table)
        add_sheet(wb, "Notes_Definitions", notes)
        add_sheet(wb, "QC_Locked_Parameters", qc_locked)
        add_sheet(wb, "Source_Audit", source_audit)
        for (sh in c("Distinct_Haplotypes", "Recurrent_Haplotypes")) {
            x <- if (sh == "Distinct_Haplotypes") {
                catalogue
            }
            else {
                recurrent_catalogue
            }
            fcol <- which(names(x) == "Empirical_haplotype_frequency")
            hcol <- which(names(x) == "HMP_contribution")
            addStyle(wb, sh, freq_style, rows = 2:(nrow(x) + 1), cols = fcol, gridExpand = TRUE)
            addStyle(wb, sh, hmp_style, rows = 2:(nrow(x) + 1), cols = hcol, gridExpand = TRUE)
            long_cols <- which(names(x) %in% c("Recurrence_scope", "Administrative_divisions", "Regions", "Predictor_Haplogroup", 
                "Predictor_probability_percent", "Predictor_fitness_score", "Deletion_pattern"))
            if (length(long_cols) > 0L) {
                setColWidths(wb, sh, long_cols, 27)
            }
            locus_cols <- which(names(x) %in% HAPLOTYPE_SYSTEMS)
            if (length(locus_cols) > 0L) {
                setColWidths(wb, sh, locus_cols, 11)
            }
        }
        saveWorkbook(wb, path, overwrite = TRUE)
    }
}
