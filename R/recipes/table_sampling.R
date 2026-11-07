# Adapted from the user's original table_sampling final builder.
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
    number_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0", halign = "center", valign = "center")
    coord_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0.000000", halign = "center", valign = "center")
    total_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, textDecoration = "bold", fgFill = "#D9EAF7", 
        border = "Top", borderColour = "#808080")
    note_style <- openxlsx::createStyle(fontName = style$font, fontSize = 9, wrapText = TRUE, valign = "top")
    openxlsx::writeData(wb, table_meta$sheet, table_meta$title, startRow = 1, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, title_style, rows = 1, cols = 1)
    openxlsx::writeData(wb, table_meta$sheet, table_data, startRow = 3, startCol = 1, headerStyle = header_style)
    n_rows <- nrow(table_data)
    n_cols <- ncol(table_data)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = 4:(3 + n_rows), cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
    integer_cols <- which(names(table_data) %in% c("N_full_forensic", "N_complete_structure", "Excluded_structural_anomaly_profiles"))
    if (length(integer_cols) > 0L) {
        openxlsx::addStyle(wb, table_meta$sheet, number_style, rows = 4:(3 + n_rows), cols = integer_cols, gridExpand = TRUE, 
            stack = TRUE)
    }
    coord_cols <- which(names(table_data) %in% c("Weighted_latitude", "Weighted_longitude"))
    if (length(coord_cols) > 0L) {
        openxlsx::addStyle(wb, table_meta$sheet, coord_style, rows = 4:(3 + n_rows), cols = coord_cols, gridExpand = TRUE, 
            stack = TRUE)
    }
    total_row <- 3 + n_rows
    openxlsx::addStyle(wb, table_meta$sheet, total_style, rows = total_row, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, table_meta$sheet, cols = 1:n_cols, widths = c(16, 24, 14, 15, 18, 18, 17, 17))
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
