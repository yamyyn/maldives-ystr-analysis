# Study methods and demonstration settings

The invented-data demonstration and the study-scale configuration serve different purposes. `config/analysis_settings.yml` uses short runs and `config/synthetic_random_seeds.csv` for testing. `config/private_settings.template.yml` records the study's methodological settings and uses `config/random_seeds.csv`. Copy the latter template into ignored `config/private/settings.yml` and provide authorized local inputs and original external results. No observed sample sizes, selected PC counts or manuscript statistics are fixed verification targets.

| Analysis | Study setting | Seed |
|---|---|---|
| Division-stratified random removal | 100,000 draws, complete-profile removal candidates | 20261002 |
| Optional paired bootstrap | 10,000 draws, QC-only | 20261001 |
| Accumulation | 500 sample orders; sample sizes 1, multiples of 10, and final N | 20261003 |
| Arlequin primary analyses | 10,000 requested permutations, external GUI | The GUI controls its randomization |
| Synthetic structure verification | Separately labelled R verification | 20261004 |
| Mantel | 10,000 division-label permutations | 20261005 |
| Correlogram | 10,000 permutations shared across distance classes | 20261006 |
| Geographic DAPC | 100 stratified holdouts, training proportion 0.90 | 20261004 |
| Exploratory parent-category DAPC and rare-group sensitivity | 100 holdouts; singleton groups training-only; rare-group threshold 5 | 20261007 |
| Rarefaction | 10,000 subsamples per group, offsets i and 100+i | 20261008 |
| Optional bootstrap uncertainty | 10,000 draws, scope offsets | 20261009 |

Geographic DAPC evaluates the configured base grid `[1,2,3,5,10,15,20,25,30,40,50,60,70,80,90,100,120,140]` below the centered full-matrix rank. Geographic group order is alphabetical. Exploratory DAPC orders categories by descending count, breaking ties alphabetically, and starts from `[17,20,25,30,35,40,50,60,70,80,90,100,110,120,130,140]`. Its rank-minus-one candidate is added from the input matrix, candidates below the full discriminant dimension are removed, and candidates unsupported by any training-only PCA are removed. Thus a manuscript's effective candidate list need not equal the broad base grid. Neither effective candidates nor selected PCs are hardcoded from study results.

Training splits retain group-block sampling order. Singleton groups are sampled into training only, preserving the random-number stream. PCA is computed on each training subset and supplied to all candidate fits for that split. The highest PC count is selected when multiple candidates have exactly equal minimum RMSE. Assignment rates remain internal estimates from the model-selection holdouts. Both the full exploratory model and the separately filtered rare-group model follow this procedure.

The workflow records configuration, seed table, RNG kind and package versions with each run. Matching seeds alone does not establish byte-identical reproduction: exact authorized input rows and their order, group labels, allele-state encoding, dependency versions and original Arlequin outputs must also match. This public package is checked with invented data, not certified against undistributed individual profiles. Do not substitute synthetic predictor records or R verification calculations for study predictor outputs or primary GUI AMOVA.

The primary Mantel seed is separate from additional coordinate-sensitivity randomization. Additional checks may use their own documented seeds; do not replace the primary run with them or copy their P-values into the primary result. Arlequin's reported permutation count is read from its output rather than assumed to equal the requested count.

Primary structure eligibility excludes profiles containing microvariants at included typed single-copy loci before DYS389II subtraction. Categorical DAPC retains its complete-profile data and microvariant states. RST inference keeps negative values; the specified display transformations do not change those inferential values. Source data and result files are never edited by analysis functions.

Public data contain invented profiles only. Laboratory identifiers, original predictor exports, private genotype crosswalks and manuscript-result verification targets are not part of the repository. `config/original_targets.template.csv` is an unfilled optional audit schema, not evidence of private numerical validation.
