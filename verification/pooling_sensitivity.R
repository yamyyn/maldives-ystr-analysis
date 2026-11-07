# Run from the repository root after stage 02.
source("R/pipeline.R")
cfg <- load_pipeline_config(Sys.getenv("YSTR_CONFIG", "config/analysis_settings.yml"))
checkpoint <- readRDS(file.path(cfg$out, "pipeline_state.rds"))
assert(identical(checkpoint$input_hash, pipeline_input_hash(cfg)), "Inputs or code changed since the pipeline run")
cohort <- checkpoint$state[["02_build_datasets"]]
assert(!is.null(cohort), "Run the pipeline through stage 02 first")
full <- cohort$full
keys <- haplotype_keys(full, c(cohort$study$single, names(cohort$study$multicopy)))
counts <- table(keys, full$atoll)
assert(ncol(counts) > 1L && sum(counts) == nrow(full), "Invalid geographic counts")
hmp <- function(n) sum((n / sum(n))^2)
pooled <- hmp(rowSums(counts))
equal <- sum(rowMeans(sweep(counts, 2L, colSums(counts), "/"))^2)
excluded <- vapply(seq_len(ncol(counts)), function(i) hmp(rowSums(counts[, -i, drop=FALSE])), numeric(1))
out <- data.frame(Scenario=c("Observed pooled", "Equal division weights", paste("Excluded",colnames(counts))),
                  HMP=c(pooled, equal, excluded))
out$Relative_change_percent <- 100 * (out$HMP / pooled - 1)
dest <- file.path(cfg$out, "verification", "Pooling_sensitivity.csv")
dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
write.csv(out, dest, row.names=FALSE)
floor_pooled <- 1 / nrow(full)
floor_equal <- sum(1 / colSums(counts)) / ncol(counts)^2
decomposition <- data.frame(Weighting=c("Sample proportional","Equal divisions"),
 HMP=c(pooled,equal),Self_pair_floor=c(floor_pooled,floor_equal),
 Distinct_individual_contribution=c(pooled-floor_pooled,equal-floor_equal))
write.csv(decomposition,file.path(dirname(dest),"HMP_SelfPair_Decomposition.csv"),row.names=FALSE)
message("Wrote descriptive pooling sensitivity to ", dest)
# Subtracting the self-pair floor alone is not an unbiased population estimator.
