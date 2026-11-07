# Adapted from the user's original table_sensitivity final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_workbook <- function(path, table_data, footnotes, qc, source_audit) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, table_meta$sheet)
    openxlsx::addWorksheet(wb, "Footnotes")
    openxlsx::addWorksheet(wb, "QC")
    openxlsx::addWorksheet(wb, "Source_Audit")
    title_style <- openxlsx::createStyle(fontName = style$font, fontSize = 11, textDecoration = "bold", fontColour = "#000000", 
        halign = "left", valign = "center")
    header_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, textDecoration = "bold", fontColour = "#FFFFFF", 
        fgFill = "#1F4E78", halign = "center", valign = "center", wrapText = TRUE)
    body_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, valign = "center", wrapText = TRUE)
    note_style <- openxlsx::createStyle(fontName = style$font, fontSize = 9, wrapText = TRUE, valign = "top")
    openxlsx::writeData(wb, table_meta$sheet, table_meta$title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, title_style, rows = 1, cols = 1)
    openxlsx::writeData(wb, table_meta$sheet, table_data, startRow = 3, startCol = 1, headerStyle = header_style)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = 4:(3 + nrow(table_data)), cols = 1:ncol(table_data), gridExpand = TRUE, 
        stack = TRUE)
    openxlsx::setColWidths(wb, table_meta$sheet, cols = 1:ncol(table_data), widths = c(30, 21, 24, 22, 18, 54))
    openxlsx::freezePane(wb, table_meta$sheet, firstActiveRow = 4)
    openxlsx::writeData(wb, "Footnotes", footnotes, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, "Footnotes", note_style, rows = 1:max(1L, nrow(footnotes) + 1L), cols = 1:max(1L, ncol(footnotes)), 
        gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, "Footnotes", cols = 1:2, widths = c(18, 110))
    openxlsx::writeData(wb, "QC", qc, startRow = 1, startCol = 1)
    openxlsx::setColWidths(wb, "QC", cols = 1:ncol(qc), widths = "auto")
    openxlsx::writeData(wb, "Source_Audit", source_audit, startRow = 1, startCol = 1)
    openxlsx::setColWidths(wb, "Source_Audit", cols = 1:ncol(source_audit), widths = "auto")
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}
