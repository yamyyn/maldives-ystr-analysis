# Adapted from the user's original table_numbering final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
safe_md5 <- function(path) {
    if (!file.exists(path)) {
        return(NA_character_)
    }
    unname(tools::md5sum(path))
}
read_sheet_raw <- function(path, sheet) {
    openxlsx::read.xlsx(path, sheet = sheet, colNames = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE, check.names = FALSE)
}
compare_workbook_content_except_first_title <- function(source_path, final_path) {
    source_sheets <- openxlsx::getSheetNames(source_path)
    final_sheets <- openxlsx::getSheetNames(final_path)
    if (length(source_sheets) != length(final_sheets)) {
        return(FALSE)
    }
    for (i in seq_along(source_sheets)) {
        a <- read_sheet_raw(source_path, source_sheets[i])
        b <- read_sheet_raw(final_path, final_sheets[i])
        nr <- max(nrow(a), nrow(b))
        nc <- max(ncol(a), ncol(b))
        pad_matrix <- function(x, nr, nc) {
            out <- matrix(NA_character_, nrow = nr, ncol = nc)
            if (nrow(x) > 0L && ncol(x) > 0L) {
                vals <- as.matrix(data.frame(lapply(x, as.character), check.names = FALSE))
                out[seq_len(nrow(vals)), seq_len(ncol(vals))] <- vals
            }
            out
        }
        aa <- pad_matrix(a, nr, nc)
        bb <- pad_matrix(b, nr, nc)
        if (i == 1L) {
            aa[1, 1] <- NA_character_
            bb[1, 1] <- NA_character_
        }
        aa[is.na(aa)] <- ""
        bb[is.na(bb)] <- ""
        if (!identical(aa, bb)) {
            return(FALSE)
        }
    }
    TRUE
}
standardize_workbook <- function(source_path, destination_path, final_sheet_name, final_title) {
    wb <- openxlsx::loadWorkbook(source_path)
    current_sheets <- names(wb)
    if (length(current_sheets) < 1L) {
        stop(paste("Workbook has no worksheets:", source_path), call. = FALSE)
    }
    current_sheets[1] <- final_sheet_name
    names(wb) <- current_sheets
    openxlsx::writeData(wb, sheet = final_sheet_name, x = final_title, startRow = 1, startCol = 1, colNames = FALSE, rowNames = FALSE)
    openxlsx::saveWorkbook(wb, destination_path, overwrite = TRUE)
}
