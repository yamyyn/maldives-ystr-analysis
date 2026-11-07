# Adapted from the user's original table_recoding final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_workbook <- function(path) {
    {
        wb <- createWorkbook()
        addWorksheet(wb, "Recoding_Audit")
        addWorksheet(wb, "Predictor_to_Final_Mapping")
        addWorksheet(wb, "Final_Haplogroup_Summary")
        addWorksheet(wb, "Notes")
        writeData(wb, "Recoding_Audit", audit_summary)
        writeData(wb, "Predictor_to_Final_Mapping", mapping)
        writeData(wb, "Final_Haplogroup_Summary", final_summary)
        writeData(wb, "Notes", tibble::tibble(Item = c("Predictor_Haplogroup", "Final_Haplogroup", "Exact label", "Parent-category recoding", 
            "Interpretation limitation", "Source"), Description = c("Detailed Y-STR predictor output.", "User-defined parent reporting category used for population-level analyses.", 
            "Normalized detailed predictor label equals the Final Haplogroup label.", "Detailed predictor output was recoded to a broader parent reporting category.", 
            "Assignments are Y-STR-derived and were not confirmed by Y-SNP genotyping.", best$source)))
        for (s in names(wb)) {
            freezePane(wb, s, firstActiveRow = 2)
            setColWidths(wb, s, 1:10, "auto")
        }
        saveWorkbook(wb, path, overwrite = TRUE)
    }
}
