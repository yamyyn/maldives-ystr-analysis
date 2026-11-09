# Version 1.2.1 — 9 November 2026

- Separated synthetic demonstration seeds from the documented study configuration; corrected the study random-removal seed and Mantel permutation count.
- Added a separately seeded correlogram, with its own iteration count and both progressive and simultaneous Holm outputs.
- Aligned DAPC group ordering, singleton RNG consumption, grouped training order, training-only PCA reuse, rank-dependent candidate filtering and highest-PC RMSE tie handling with the study scripts.
- Preserved the broad exploratory base grid. Effective grids and selected PCs remain input-derived rather than copied from manuscript results.
- Matched accumulation evaluation sizes to the study's multiples-of-ten grid and retained separate rare-group sensitivity models.
- Recorded effective configuration/seeds/RNG metadata with each run; added invented-data reproducibility regressions.
- Updated the predicted-haplogroup frequency axis label and source recipe provenance.
- Retained the external Arlequin GUI requirement and existing private-data exclusions. No private profiles or manuscript-result targets were added.

# Version 1.2.0 — 8 November 2026

- Replaced allele-size microvariant truncation with profile exclusion at included loci, assessed before DYS389II adjustment.
- Separated Arlequin-eligible metadata from the complete categorical DAPC cohort, preserving microvariant states in DAPC.
- Added data-derived eligibility audits and explicit cohort counts; corrected panel-specific filtering in leave-one-division-out verification.
- Added synthetic regression cases for alignment, omitted-marker eligibility, no-call handling and partial-repeat cancellation.
- Updated methodology and external GUI-run instructions. No private genotypes, study-specific expected results or external executable were added.

# Version 1.1.1 - 7 November 2026

- Removed an inherited non-default renv profile selector that referred to an absent profile-specific lockfile. Dependency restoration now uses the supplied root renv.lock.
- No analysis functions, scientific settings, example inputs or package versions changed.

# Version 1.1.0 — 7 November 2026

- Finalized public source archive and GPL-3.0-only licence; updated citation authors.
- Excluded cached outputs, external run files, local libraries and private material.
- Corrected random-removal eligibility to complete profiles, preserving division-matched removal and empirical tail calculations.
- Aligned full exploratory haplogroup DAPC with two-profile validation eligibility and a distinct five-profile rare-group sensitivity cohort.
- Matched regional rarefaction seed offsets to the documented method.
- Added the HMP self-pair decomposition to the separate pooling verification script.
- Updated execution, privacy, release and verification documentation. Primary AMOVA remains a separate Arlequin GUI step.

These changes prepare reusable code. They do not certify exact equality to every private manuscript output or historical figure design.

- Corrected the inherited DYS385 repeated-equal-slot bug through one shared observed-state normalizer, reused for import and locus-frequency table preparation; retained raw slots and DYF387S1 observations. Updated adapted-recipe hashes while retaining origin hashes.

- Aligned final figure/table numbering and Y-GATA-H4 typography while retaining historical labels. Renamed the derived predictor threshold flag to avoid implying it was an observed predictor automatic call.
