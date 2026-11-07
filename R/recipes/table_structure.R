# Adapted from the user's original table_structure final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
write_table_workbook <- function(path, amova, pairwise, mantel, footnotes, qc, source_audit) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, table_meta$sheet)
    openxlsx::addWorksheet(wb, "Footnotes")
    openxlsx::addWorksheet(wb, "QC")
    openxlsx::addWorksheet(wb, "Source_Audit")
    title_style <- openxlsx::createStyle(fontName = style$font, fontSize = 11, textDecoration = "bold", fontColour = "#000000", 
        halign = "left", valign = "center")
    section_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, textDecoration = "bold", fontColour = "#FFFFFF", 
        fgFill = "#1F4E78", halign = "left", valign = "center")
    header_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, textDecoration = "bold", fontColour = "#000000", 
        fgFill = "#D9EAF7", halign = "center", valign = "center", wrapText = TRUE, border = "Bottom", borderColour = "#808080")
    body_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, valign = "center", wrapText = TRUE)
    numeric_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0.00000", halign = "center", valign = "center")
    percent_style <- openxlsx::createStyle(fontName = style$font, fontSize = 10, numFmt = "0.00", halign = "center", valign = "center")
    note_style <- openxlsx::createStyle(fontName = style$font, fontSize = 9, wrapText = TRUE, valign = "top")
    row <- 1L
    openxlsx::writeData(wb, table_meta$sheet, table_meta$title, startRow = row, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, title_style, rows = row, cols = 1)
    row <- row + 2L
    openxlsx::mergeCells(wb, table_meta$sheet, cols = 1:5, rows = row)
    openxlsx::writeData(wb, table_meta$sheet, "A. Hierarchical AMOVA", startRow = row, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, section_style, rows = row, cols = 1:5, gridExpand = TRUE)
    row <- row + 1L
    amova_start <- row
    openxlsx::writeData(wb, table_meta$sheet, amova, startRow = row, startCol = 1, headerStyle = header_style)
    amova_data_start <- row + 1L
    amova_data_end <- row + nrow(amova)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = amova_data_start:amova_data_end, cols = 1:5, gridExpand = TRUE, 
        stack = TRUE)
    openxlsx::addStyle(wb, table_meta$sheet, numeric_style, rows = amova_data_start:amova_data_end, cols = c(2, 4), gridExpand = TRUE, 
        stack = TRUE)
    openxlsx::addStyle(wb, table_meta$sheet, percent_style, rows = amova_data_start:amova_data_end, cols = 3, gridExpand = TRUE, 
        stack = TRUE)
    row <- amova_data_end + 2L
    openxlsx::mergeCells(wb, table_meta$sheet, cols = 1:7, rows = row)
    openxlsx::writeData(wb, table_meta$sheet, "B. Pairwise RST summary", startRow = row, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, section_style, rows = row, cols = 1:7, gridExpand = TRUE)
    row <- row + 1L
    pairwise_start <- row
    openxlsx::writeData(wb, table_meta$sheet, pairwise, startRow = row, startCol = 1, headerStyle = header_style)
    pairwise_data_start <- row + 1L
    pairwise_data_end <- row + nrow(pairwise)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = pairwise_data_start:pairwise_data_end, cols = 1:7, gridExpand = TRUE, 
        stack = TRUE)
    row <- pairwise_data_end + 2L
    openxlsx::mergeCells(wb, table_meta$sheet, cols = 1:6, rows = row)
    openxlsx::writeData(wb, table_meta$sheet, "C. Isolation-by-distance Mantel test", startRow = row, startCol = 1)
    openxlsx::addStyle(wb, table_meta$sheet, section_style, rows = row, cols = 1:6, gridExpand = TRUE)
    row <- row + 1L
    mantel_start <- row
    openxlsx::writeData(wb, table_meta$sheet, mantel, startRow = row, startCol = 1, headerStyle = header_style)
    mantel_data_start <- row + 1L
    mantel_data_end <- row + nrow(mantel)
    openxlsx::addStyle(wb, table_meta$sheet, body_style, rows = mantel_data_start:mantel_data_end, cols = 1:6, gridExpand = TRUE, 
        stack = TRUE)
    openxlsx::setColWidths(wb, table_meta$sheet, cols = 1:7, widths = c(34, 18, 18, 18, 20, 22, 58))
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
