# Release scope

This is a code-only workflow release with invented data. It includes mathematical checks and synthetic integration tests, not hard-coded manuscript verification results. Methodological constants (assay panels, thresholds, seeds, iteration settings) remain explicit configuration.

The actual GitHub URL and any software DOI must be added after publication. A clean dependency restore and hosted CI require external network access. Local test success is not a claim that either has occurred.

The public code does not bundle or run a contemporary haplogroup predictor, derive confirmed deletions from a single observed allele, or estimate a casework theta/discrete-Laplace/kappa model. Supplied prediction support remains exploratory. Automatically derived threshold flags must be checked against the originating predictor batch, not represented as independently observed automatic assignments.

The `threshold_qualified` prediction field is a derived parent-support/leading-fitness rule, not an observed automatic assignment from the predictor. Preserve the predictor batch output separately if reporting actual automatic assignment counts.

The registry preserves historical locked labels and supplies final labels for locus diversity (Figure 4), predicted haplogroups (Figure 5), PCoA (Figure 6), haplogroup table (Table 5) and AMOVA table (Table 6). Main allele-frequency figures are generated as additional plate exports; the sampling map remains external.
