# Copy into config/local.R and configure for your own assay and private files.
# Source the standard panel only if it matches your assay. Do not source the
# synthetic config for a real study: its evidence files are demonstration-only.
qc_config <- initial_qc_defaults(
  single=c("replace_with_single_copy_column_names"),
  multicopy=list(replace_state_column=c("replace_slot_a","replace_slot_b")))
qc_input <- "data/private/master.csv"
qc_output <- "outputs/private_initial_qc"
# qc_config$population_mapping <- read_qc_csv("data/private/population_mapping.csv")
# qc_config$confirmed_patterns <- read_qc_csv("data/private/confirmed_patterns.csv")
# Confirmation columns: sample_id, loci (joined with |), evidence, classification.
# classification: confirmed_deletion or confirmed_no_call. An absent
# classification column never establishes deletion membership.
# Optional independent expectations: expected_n, expected_ids, id_pattern,
# expected_populations, expected_regions, allele_bounds (field/min/max).
# All metadata fields are required by default. If genuinely unavailable, set
# metadata_required=FALSE and metadata_constant=character(); document the gap.
