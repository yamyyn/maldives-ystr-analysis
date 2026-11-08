# Optional external GUI input preparation after stage 09a. No Arlequin execution.
source("R/pipeline.R")
require_packages()
cfg <- load_pipeline_config(Sys.getenv("YSTR_CONFIG", "config/analysis_settings.yml"))
checkpoint <- readRDS(file.path(cfg$out,"pipeline_state.rds"))
assert(identical(checkpoint$input_hash,pipeline_input_hash(cfg)),"Inputs/code changed; start a new run")
data <- checkpoint$state[["02_build_datasets"]]
assert(!is.null(data),"Build datasets first")
single <- data$study$single
non_rm <- setdiff(single,c("dys449","dys518","dys570","dys576","dys627"))
unaffected <- sub("^dys389iiadj$","dys389ii",unlist(cfg$haplogroup$unaffected_loci))
assert(all(unaffected %in% single),"Configure unaffected loci as typed single-copy fields")
models <- list(
  single_copy_complete=prepare_structure(data$complete,single),
  non_rm_complete=prepare_structure(data$complete,non_rm),
  unaffected_complete=prepare_structure(data$complete,unaffected),
  unaffected_full=prepare_structure(data$full,unaffected))
capital_code <- Sys.getenv("YSTR_CAPITAL_CODE",unset="MLE")
central_region <- Sys.getenv("YSTR_CENTRAL_REGION",unset="Central")
regional_available <- capital_code %in% models$single_copy_complete$arlequin_metadata$atoll &&
                      central_region %in% models$single_copy_complete$arlequin_metadata$region
if(!regional_available)message("Regional alternatives not exported: configure YSTR_CAPITAL_CODE and YSTR_CENTRAL_REGION for this dataset.")
for (panel in if(regional_available)c("single_copy_complete","non_rm_complete") else character()) {
  baseline <- models[[panel]]
  for (scheme in c("capital_to_central","capital_excluded")) {
    z <- baseline
    assert(capital_code %in% z$arlequin_metadata$atoll && central_region %in% z$arlequin_metadata$region,
           "Configure the capital code/group for this dataset before using this study-specific alternative")
    keep <- if(scheme=="capital_excluded")z$arlequin_metadata$atoll!=capital_code else rep(TRUE,nrow(z$integer))
    original_rows <- which(z$eligibility$Eligible)[keep]
    z$integer <- z$integer[keep,,drop=FALSE]
    z$arlequin_metadata <- z$arlequin_metadata[keep,,drop=FALSE]
    if(scheme=="capital_to_central")z$arlequin_metadata$region[z$arlequin_metadata$atoll==capital_code] <- central_region
    z$recoding <- z$recoding[z$recoding$Row %in% original_rows,,drop=FALSE]
    models[[paste(panel,scheme,sep="__")]] <- z
  }
}
export_id <- Sys.getenv("YSTR_SENSITIVITY_EXPORT_ID",unset=paste0(cfg$run_id,"_filtered_sensitivity"))
assert(grepl("^[A-Za-z0-9_-]+$",export_id),"Unsafe export ID")
out <- file.path("external","arlequin",export_id)
assert(!dir.exists(out),"Sensitivity export exists; preserve prior inputs or choose a new run ID")
rows <- list()
for (name in names(models)) {
  z <- models[[name]];dest <- file.path(out,name)
  export_arlequin(z,dest)
  rows[[length(rows)+1L]] <- data.frame(Model=name,Eligible_N=nrow(z$integer),Variables=ncol(z$integer),
    Populations=length(unique(z$arlequin_metadata$atoll)),Regions=length(unique(z$arlequin_metadata$region)),
    GUI_Tasks=if(grepl("capital_",name))"Hierarchical AMOVA only required; pairwise exports optional" else "Hierarchical AMOVA and division/regional pairwise RST")
}
write_csv(do.call(rbind,rows),file.path(out,"MODEL_MANIFEST.csv"))
writeLines(c("Inputs only. Primary analyses must be run separately in Arlequin GUI.",
  "Use haplotypic microsatellite squared allele-size differences and configured permutation count.",
  "The first four models provide three project files each; regional alternatives require hierarchical AMOVA only.",
  "Retain .ars settings, original logs and XML reports, including the initialized RNG values.",
  "Capital MLE/Central alternatives reproduce this workflow's grouping design; adapt these labels explicitly for other studies.",
  "Eligibility is recalculated per included marker panel. Exported integers are not truncated.",
  "Private inputs and external results stay ignored by Git; do not add genetic data to the repository."),
  file.path(out,"GUI_RUN_INSTRUCTIONS.txt"))
message("Prepared sensitivity projects in ",out)
