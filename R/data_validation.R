# Public input-validation entry point; shared with the stage runner.
validate_input_file <- function(input,settings,report_directory=NULL) {
  raw<-read_qc_csv(input)
  result<-run_initial_qc(raw,settings)
  if(!is.null(report_directory))write_initial_qc(result,safe_path(report_directory,FALSE),input)
  result
}
