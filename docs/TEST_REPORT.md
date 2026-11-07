# Release verification — version 1.1.0

Verified locally on Windows using R 4.6.1 and installed dependencies. All 21 synthetic workflow stages completed and 235 automated checks passed. All R source files parsed successfully. Stage logs contained no analysis warnings; host startup emitted locale-setting warnings.

The checks include DYS385 single-state handling, original-call preservation, linked-haplotype arithmetic, complete-profile-only random-removal eligibility, resampling repeatability, train-only DAPC cross-validation, separate rare-group eligibility, microvariant coding, parser invariants, output formats, and the withheld synthetic publication gate. Separate pooling, local-profile sensitivity and leave-one-division-out AMOVA point-estimate verification scripts were executed successfully on invented input.

The public archive contains source code, templates, documentation and invented input only. Cached results, local libraries, external result files, private records and private verification targets are excluded. Personal-path and credential-pattern scans found no matches. The source manifest and archive SHA-256 provide integrity checks, not proof of scientific equivalence to private data.

A clean network dependency restore, another operating system, GitHub-hosted CI, a private Arlequin GUI run and numerical/visual equality to every manuscript output were not verified. The primary AMOVA remains a separate Arlequin GUI step. This repository is a workflow source package, not an installable R library package. GPL-3.0-only was selected by the study author.
