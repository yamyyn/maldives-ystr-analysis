# External Arlequin checkpoint

Stage 09a preserves microvariant allele states for categorical DAPC. Its separate allele-size branch excludes any profile with an intermediate allele at an included locus; no truncation or rounding of allele calls is performed. Eligibility is evaluated before DYS389II minus DYS389I, and the audit records excluded loci and cohort sizes. Stage 09b writes hierarchical AMOVA, population pairwise and region pairwise `.arp` projects with original source-compatible microsatellite settings and validated profile/sample blocks. Sanitized labels are checked for collisions and mapped back to canonical labels.

Run the three projects in Arlequin externally. Use haplotypic microsatellite data, sum of squared size differences (RST), the configured permutation setting and population nesting from the export. Preserve software version, project files, settings, raw XML/text output and checksums. Do not edit labels or transplant results from another cohort.

Configure the expected `amova`, `atoll` and `region` result paths before starting a private run. The runner stops if they are absent at 09c. Supply the files in those paths and resume without changing the frozen configuration. The source-derived parser expects the original Arlequin report tags and text conventions; incompatible output formats fail rather than silently becoming missing values.

The importer checks cohort labels, sample count, data type, haplotypic setting, usable locus count, requested permutations, matrix symmetry, diagonals and P-value ranges. Raw RST is retained; Bonferroni and BH are computed separately. Manually compare reported permutations, AMOVA degrees of freedom and allele-size-eligible cohort sizes to the exported project as part of private review. The external executable is not bundled or launched automatically.

Synthetic mode follows an explicitly different branch: R method-of-moments estimates and permutations exercise later analysis and graphics. They are labelled `not Arlequin` in tables/settings and cannot satisfy the external-result publication gate.

## Arlequin graphical interface

Run the workflow through stage 09b. Open each exported .arp in Arlequin, enable hierarchical AMOVA for hierarchical_amova.arp and pairwise population differentiation for the other two projects, and save the settings as .ars alongside the projects. Use the configured permutation count and squared microsatellite repeat differences. Run the analyses in Arlequin and put the resulting XML reports at the three paths configured before the run. Resume at 09c. The R workflow never launches Arlequin. Its synthetic smoke test uses a separately labelled R verification branch.

## Panel-specific eligibility

Reassess eligibility for every marker subset using its included typed calls. Do not build non-RM or unaffected-locus cohorts by taking columns from the already filtered full panel: a profile excluded only for an omitted marker can become eligible. Missing/unspecified and invalid calls require the initial QC rules; included microvariants are not converted to integer alleles. DYS385 and DYF387S1 are outside the single-copy population-structure panel. Categorical DAPC retains complete-profile microvariants as distinct states. The YHRD AMOVA service documents intermediate-allele exclusion: https://yhrd.org/AMOVA. This is not a claim that every YHRD function uses the same filtering.

Run `Rscript verification/export_filtered_sensitivity_projects.R` after the checkpoint to export marker-specific full, non-RM and unaffected-locus GUI projects. Regional alternatives use `YSTR_CAPITAL_CODE` (default MLE) and `YSTR_CENTRAL_REGION` (default Central); when those labels are absent, the helper explicitly reports that alternatives were not exported. Set `YSTR_SENSITIVITY_EXPORT_ID` to a new safe ID when preserving an earlier export. These files contain data and remain ignored by Git.
