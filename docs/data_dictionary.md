# Data dictionary

## Master CSV

All fields are imported as text. Canonical columns are `sample_id`, `atoll`, `region`, `latitude_centroid`, `longitude_centroid`, `latitude_weighted`, `longitude_weighted`, `islands_sampled`, `total_samples`, the configured single-copy alleles and multicopy observed slots. Use an explicit source-to-canonical mapping for other headers. IDs must be nonblank and unique. Repeated population metadata must agree. Coordinates are supplied, not repaired.

Decimal allele strings represent observed microvariants. Blank optional multicopy slots do not imply absent genomic copies. The configured no-call token is permitted only with matching independent single-copy confirmation; no-call patterns are not automatically called deletions.

## QC confirmation CSV

`sample_id,loci,evidence,classification`: `loci` is a `|`-separated exact no-call set; evidence is a nonblank reference. Classification is `confirmed_deletion` or `confirmed_no_call`. Without that column, confirmation is no-call only. The mapping CSV uses `population,region`. Optional laboratory bounds use numeric `field,min,max` configuration values.

## Predictor leaf CSV

`sample_id,leaf,probability_pct,fitness`: one row per modeled leaf per sample, including zeros; unique sample/leaf pairs. Probabilities are numeric 0–100 and sum to 100 within configured rounding tolerance. All IDs must match the full cohort. Fitness must be numeric and use the documented predictor scale.

`parent_mapping.csv`: unique `leaf,parent` mapping, including every modeled class. Private studies supply their own revised mapping; synthetic labels are not biological classifications.

`reported_parents.csv`: unique `sample_id,parent` rows for the full cohort, recording the independently reviewed reported parent. These labels are preserved when posterior and fitness reporting is recalculated.

## Private validation

`original_targets.csv`: `key,comparison,expected,tolerance`, using the original audit keys and supported comparisons documented in `PRIVATE_REPRODUCTION.md`. `original_observations.csv`: `key,value,evidence`, containing fresh measurements. Do not populate observations by copying expected values. Missing configuration/observations remains incomplete. Private target values, IDs, predictor files and reports must never be committed.
