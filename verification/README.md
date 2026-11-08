# Standalone verification
Run from the repository root after stage 09a. Set YSTR_CONFIG to use a configuration other than the synthetic default.

    Rscript verification/amova_leave_one_division_out.R

This script calculates independent AMOVA point estimates for the single-copy and non-RM panels, with each division omitted in turn. It requires at least two remaining regions and more remaining populations than regions. The primary AMOVA and its permutation results use the Arlequin GUI handoff. This check does not calculate permutation P-values.

For the descriptive pooling sensitivity:

    Rscript verification/pooling_sensitivity.R

This compares sample-weighted pooling, equal division weights, and exclusion of each division in turn.

Run Rscript verification/local_profile_sensitivity.R after the full pipeline for full-versus-complete division rarefaction at a common feasible size and the excluded-profile haplotype counts. The pooling script also exports the self-pair floor decomposition; this is descriptive and is not an unbiased estimator.

Version 1.2.0: allele-size cohorts exclude included-locus intermediate calls. LODO reassesses eligibility independently for its full and non-RM panels. `export_filtered_sensitivity_projects.R` prepares reduced-panel and configurable capital-placement GUI inputs; it never runs Arlequin.
