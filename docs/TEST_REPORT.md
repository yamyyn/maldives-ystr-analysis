# Release verification — version 1.2.1

Verified on Windows with R 4.6.1 on 9 November 2026. All 21 synthetic workflow stages completed and all 265 automated checks passed. All 102 R source files parsed. Installed dependency versions were checked against every package in renv.lock and agreed, using existing matching installations; a fresh network restore was not performed locally.

New invented-data tests independently check grouped training splits, singleton training-only handling and random-number consumption, training-only PCA prediction, input-derived PC filtering, maximum-PC selection on an RMSE tie, and independent Mantel/correlogram seeds and iteration counts. Method-setting checks distinguish the short demonstration from the study template. No private genotype, observed sample-size target, selected-PC target or manuscript P-value was used.

Separate AMOVA leave-one-division-out point estimates, filtered sensitivity GUI exports and pooling sensitivity helpers completed on invented data. The source checksum manifest passed. Scans for study identifier patterns, personal Windows paths and credential patterns found no matches. The stage log recorded no R warnings; the console contained single-observation geom_line messages and host-startup locale warnings.

Primary AMOVA remains an external Arlequin GUI step. Neither this test run nor the source manifest certifies equality to private manuscript results, original predictor records or private laboratory evidence. Local verification did not include a private GUI run or every private sensitivity analysis. Hosted CI status should be checked for the particular commit on GitHub.

## Earlier verification records

### Release verification â€” version 1.1.0

Verified locally on Windows using R 4.6.1 and installed dependencies. All 21 synthetic workflow stages completed and 235 automated checks passed. All R source files parsed successfully. Stage logs contained no analysis warnings; host startup emitted locale-setting warnings.

The checks include DYS385 single-state handling, original-call preservation, linked-haplotype arithmetic, complete-profile-only random-removal eligibility, resampling repeatability, train-only DAPC cross-validation, separate rare-group eligibility, microvariant coding, parser invariants, output formats, and the withheld synthetic publication gate. Separate pooling, local-profile sensitivity and leave-one-division-out AMOVA point-estimate verification scripts were executed successfully on invented input.

The public archive contains source code, templates, documentation and invented input only. Cached results, local libraries, external result files, private records and private verification targets are excluded. Personal-path and credential-pattern scans found no matches. The source manifest and archive SHA-256 provide integrity checks, not proof of scientific equivalence to private data.

A clean network dependency restore, another operating system, GitHub-hosted CI, a private Arlequin GUI run and numerical/visual equality to every manuscript output were not verified. The primary AMOVA remains a separate Arlequin GUI step. This repository is a workflow source package, not an installable R library package. GPL-3.0-only was selected by the study author.

## GitHub environment patch
The initial hosted run installed R successfully but stopped during renv setup because an inherited profile selector pointed away from the supplied lockfile. Version 1.1.1 removes that selector. Analysis code and dependencies remain identical to the locally tested version; the hosted rerun is assessed separately.

## Version 1.2.0 method update — 8 November 2026

All 21 synthetic pipeline stages completed locally. The full regression suite passed 244 checks, including invented-data cases for panel-specific exclusion, categorical/allele-size row alignment, unchanged retained integers, invalid calls and partial-repeat cancellation before DYS389II adjustment. Separate LODO and sensitivity-export helpers were run on synthetic data, including configured alternative geographic labels.

The tests assert data-derived invariants and invented fixture arithmetic; no private study genotypes or manuscript numerical targets are included. Synthetic R calculations remain labelled as verification, not external Arlequin results. Primary private AMOVA requires the documented GUI checkpoint.
