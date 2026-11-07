# External Arlequin checkpoint

Stage 09a preserves adjusted decimal alleles and audits separate integer truncation. Stage 09b writes hierarchical AMOVA, population pairwise and region pairwise `.arp` projects with original source-compatible microsatellite settings and validated profile/sample blocks. Sanitized labels are checked for collisions and mapped back to canonical labels.

Run the three projects in Arlequin externally. Use haplotypic microsatellite data, sum of squared size differences (RST), the configured permutation setting and population nesting from the export. Preserve software version, project files, settings, raw XML/text output and checksums. Do not edit labels or transplant results from another cohort.

Configure the expected `amova`, `atoll` and `region` result paths before starting a private run. The runner stops if they are absent at 09c. Supply the files in those paths and resume without changing the frozen configuration. The source-derived parser expects the original Arlequin report tags and text conventions; incompatible output formats fail rather than silently becoming missing values.

The importer checks cohort labels, sample count, data type, haplotypic setting, usable locus count, requested permutations, matrix symmetry, diagonals and P-value ranges. Raw RST is retained; Bonferroni and BH are computed separately. Manually compare reported permutations, AMOVA degrees of freedom and complete-cohort sizes to the exported project as part of private review. The external executable is not bundled or launched automatically.

Synthetic mode follows an explicitly different branch: R method-of-moments estimates and permutations exercise later analysis and graphics. They are labelled `not Arlequin` in tables/settings and cannot satisfy the external-result publication gate.

## Arlequin graphical interface

Run the workflow through stage 09b. Open each exported .arp in Arlequin, enable hierarchical AMOVA for hierarchical_amova.arp and pairwise population differentiation for the other two projects, and save the settings as .ars alongside the projects. Use the configured permutation count and squared microsatellite repeat differences. Run the analyses in Arlequin and put the resulting XML reports at the three paths configured before the run. Resume at 09c. The R workflow never launches Arlequin. Its synthetic smoke test uses a separately labelled R verification branch.
