# Adapted from the user's original table_forensic final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_table_workbook <- function(path, table_data, footnotes, qc, source_audit) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, table_meta$sheet)
    openxlsx::addWorksheet(wb, "Footnotes")
    openxlsx::addWorksheet(wb, "QC")
    openxlsx::addWorksheet(wb, "Source_Audit")
    title_style <- openxlsx::createStyle(fontName = style$font, fontSize = 11, textDecoration = "bold", fontColour = "#000000", 
        halign = "left", valign = "center")
    header_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, textDecoration = "bold", fontColour = "#FFFFFF", 
        fgFill = "#1F4E78", halign = "center", valign = "center", wrapText = TRUE, border = "Bottom", borderColour = "#D9E1F2")
    body_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, valign = "center")
    integer_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0", halign = "center", valign = "center")
    statistic_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0.0000000000", halign = "center", 
        valign = "center")
    note_style <- openxlsx::createStyle(fontName = style$font, fontSize = 9, wrapText = TRUE, valign = "top")
    openxlsx::writeData(wb, table_meta$sheet, table_meta$title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, title_style, rows = 1, cols = 1)
    openxlsx::writeData(wb, table_meta$sheet, table_data, startRow = 3, startCol = 1, headerStyle = header_style)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = 4, cols = 1:ncol(table_data), gridExpand = TRUE, stack = TRUE)
    integer_cols <- which(names(table_data) %in% c("N", "Analytical systems", "Distinct haplotypes", "Singleton haplotypes"))
    openxlsx::addStyle(wb, table_meta$sheet, integer_style, rows = 4, cols = integer_cols, gridExpand = TRUE, stack = TRUE)
    statistic_cols <- which(names(table_data) %in% c("HMP", "HD", "DC"))
    openxlsx::addStyle(wb, table_meta$sheet, statistic_style, rows = 4, cols = statistic_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, table_meta$sheet, cols = 1:ncol(table_data), widths = c(22, 10, 18, 20, 20, 16, 16, 16))
    openxlsx::freezePane(wb, table_meta$sheet, firstActiveRow = 4)
    openxlsx::writeData(wb, "Footnotes", footnotes, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, "Footnotes", note_style, rows = 1:max(1L, nrow(footnotes) + 1L), cols = 1:max(1L, ncol(footnotes)), 
        gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, "Footnotes", cols = 1:2, widths = c(18, 100))
    openxlsx::writeData(wb, "QC", qc, startRow = 1, startCol = 1)
    openxlsx::setColWidths(wb, "QC", cols = 1:ncol(qc), widths = "auto")
    openxlsx::writeData(wb, "Source_Audit", source_audit, startRow = 1, startCol = 1)
    openxlsx::setColWidths(wb, "Source_Audit", cols = 1:ncol(source_audit), widths = "auto")
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}
