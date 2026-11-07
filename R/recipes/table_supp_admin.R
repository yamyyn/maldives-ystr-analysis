# Adapted from the user's original table_supp_admin final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_workbook <- function(path) {
    {
        wb <- createWorkbook()
        addWorksheet(wb, "Admin_Forensic_Parameters")
        addWorksheet(wb, "Notes")
        writeData(wb, "Admin_Forensic_Parameters", tab)
        writeData(wb, "Notes", notes)
        freezePane(wb, "Admin_Forensic_Parameters", firstActiveRow = 2)
        setColWidths(wb, "Admin_Forensic_Parameters", 1:ncol(tab), "auto")
        setColWidths(wb, "Notes", 1:2, "auto")
        saveWorkbook(wb, path, overwrite = TRUE)
    }
}
