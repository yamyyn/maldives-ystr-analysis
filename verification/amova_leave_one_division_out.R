# Run after 09a. All comparisons are recalculated from this run's inputs.
source("R/pipeline.R")
cfg <- load_pipeline_config(Sys.getenv("YSTR_CONFIG", "config/analysis_settings.yml"))
checkpoint <- readRDS(file.path(cfg$out, "pipeline_state.rds"))
assert(identical(checkpoint$input_hash, pipeline_input_hash(cfg)), "Inputs or code changed since the pipeline run")
structure <- checkpoint$state[["09a_prepare_structure"]]
assert(!is.null(structure), "Run the pipeline through stage 09a first")
source_data <- checkpoint$state[["02_build_datasets"]]$complete
single <- checkpoint$state[["02_build_datasets"]]$study$single
rapid <- c("dys449", "dys518", "dys570", "dys576", "dys627")
panels <- list(single_copy=single, non_rm=setdiff(single, rapid))
models <- lapply(names(panels), function(panel) {
  filtered <- prepare_structure(source_data, panels[[panel]])
  x <- filtered$integer; meta <- filtered$arlequin_metadata
  do.call(rbind, lapply(c("", sort(unique(meta$atoll))), function(excluded) {
    keep <- meta$atoll != excluded
    fit <- amova_point(x[keep,,drop=FALSE], meta$atoll[keep], meta$region[keep])
    data.frame(Panel=panel, Excluded_Division=excluded, as.list(fit), check.names=FALSE)
  }))
})
out <- do.call(rbind, models)
assert(all(abs(out$percent_region + out$percent_population + out$percent_within - 100) < 1e-7),
       "AMOVA shares do not sum to 100%")
dest <- file.path(cfg$out, "verification", "AMOVA_leave_one_division_out.csv")
dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
write.csv(out, dest, row.names=FALSE)
message("Wrote ", nrow(out), " models to ", dest)
