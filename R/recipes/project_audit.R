# Adapted from the user's original project_audit final builder.
# Adaptations: configurable labels, font, palettes, axes, cohort-neutral sheets and matrix boundaries; see documentation.
add_check <- function(Section, Check, Expected, Observed, Status = c("PASS", "WARN", "FAIL"), Critical = FALSE, Evidence = NA_character_, 
    Notes = NA_character_) {
    Status <- match.arg(Status)
    audit_rows[[length(audit_rows) + 1L]] <<- tibble(Section = as.character(Section), Check = as.character(Check), Expected = as.character(Expected), 
        Observed = as.character(Observed), Status = Status, Critical = as.logical(Critical), Evidence = as.character(Evidence), 
        Notes = as.character(Notes))
    invisible(NULL)
}
check_equal <- function(section, check, observed, expected, critical = TRUE, evidence = NA_character_, notes = NA_character_) {
    pass <- identical(as.character(observed), as.character(expected))
    add_check(Section = section, Check = check, Expected = expected, Observed = observed, Status = if (pass) 
        "PASS"
    else if (critical) 
        "FAIL"
    else "WARN", Critical = critical, Evidence = evidence, Notes = notes)
    pass
}
check_near <- function(section, check, observed, expected, tolerance, critical = TRUE, evidence = NA_character_, notes = NA_character_) {
    pass <- is.finite(observed) && abs(observed - expected) <= tolerance
    add_check(Section = section, Check = check, Expected = paste0(expected, " \302\261 ", tolerance), Observed = observed, 
        Status = if (pass) 
            "PASS"
        else if (critical) 
            "FAIL"
        else "WARN", Critical = critical, Evidence = evidence, Notes = notes)
    pass
}
check_file <- function(section, check, path, critical = TRUE, notes = NA_character_) {
    exists <- file.exists(path)
    add_check(Section = section, Check = check, Expected = "File exists", Observed = if (exists) 
        "EXISTS"
    else "MISSING", Status = if (exists) 
        "PASS"
    else if (critical) 
        "FAIL"
    else "WARN", Critical = critical, Evidence = path, Notes = notes)
    exists
}
write_workbook <- function(path) {
    {
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Audit")
        openxlsx::writeData(wb, "Audit", audit)
        openxlsx::addWorksheet(wb, "Status_Summary")
        openxlsx::writeData(wb, "Status_Summary", status_summary)
        openxlsx::addWorksheet(wb, "Section_Summary")
        openxlsx::writeData(wb, "Section_Summary", section_summary)
        openxlsx::addWorksheet(wb, "Locked_Values")
        openxlsx::writeData(wb, "Locked_Values", locked_values)
        if (exists("coord_table") && nrow(coord_table) > 0L) {
            openxlsx::addWorksheet(wb, "Coordinate_Evidence")
            openxlsx::writeData(wb, "Coordinate_Evidence", coord_table)
        }
        openxlsx::addWorksheet(wb, "Run_Info")
        run_info <- tibble(Item = c("Script", "Run time", "Project root", "Overall status", "PASS checks", "WARN checks", 
            "FAIL checks", "Critical FAIL checks", "Purpose"), Value = c(SCRIPT_NAME, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
            root, overall_status, sum(audit$Status == "PASS"), sum(audit$Status == "WARN"), sum(audit$Status == "FAIL"), 
            critical_fail_n, "Read-only final analytical audit before manuscript-number lock."))
        openxlsx::writeData(wb, "Run_Info", run_info)
        header_style <- openxlsx::createStyle(fgFill = "#1F4E78", fontColour = "#FFFFFF", textDecoration = "bold", halign = "center", 
            valign = "center", border = "Bottom", borderColour = "#D9E1F2")
        pass_style <- openxlsx::createStyle(fgFill = "#E2F0D9", fontColour = "#006100")
        warn_style <- openxlsx::createStyle(fgFill = "#FFF2CC", fontColour = "#9C6500")
        fail_style <- openxlsx::createStyle(fgFill = "#F4CCCC", fontColour = "#9C0006", textDecoration = "bold")
        sheet_ncols <- c(Audit = ncol(audit), Status_Summary = ncol(status_summary), Section_Summary = ncol(section_summary), 
            Locked_Values = ncol(locked_values), Run_Info = ncol(run_info))
        if (exists("coord_table") && nrow(coord_table) > 0L) {
            sheet_ncols <- c(sheet_ncols, Coordinate_Evidence = ncol(coord_table))
        }
        for (sheet in names(sheet_ncols)) {
            n_cols <- unname(sheet_ncols[[sheet]])
            if (is.na(n_cols) || n_cols < 1L) 
                next
            openxlsx::addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(n_cols), gridExpand = TRUE, stack = TRUE)
            openxlsx::freezePane(wb, sheet, firstRow = TRUE)
            openxlsx::setColWidths(wb, sheet, cols = seq_len(n_cols), widths = "auto")
        }
        status_col <- which(names(audit) == "Status")
        if (length(status_col) == 1L) {
            pass_rows <- which(audit$Status == "PASS") + 1L
            warn_rows <- which(audit$Status == "WARN") + 1L
            fail_rows <- which(audit$Status == "FAIL") + 1L
            if (length(pass_rows) > 0L) {
                openxlsx::addStyle(wb, "Audit", pass_style, rows = pass_rows, cols = status_col, gridExpand = TRUE, stack = TRUE)
            }
            if (length(warn_rows) > 0L) {
                openxlsx::addStyle(wb, "Audit", warn_style, rows = warn_rows, cols = status_col, gridExpand = TRUE, stack = TRUE)
            }
            if (length(fail_rows) > 0L) {
                openxlsx::addStyle(wb, "Audit", fail_style, rows = fail_rows, cols = status_col, gridExpand = TRUE, stack = TRUE)
            }
        }
        openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    }
}
