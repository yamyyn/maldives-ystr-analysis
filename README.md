# Maldives Y-STR analysis workflow

Version 1.2.1

R scripts for Y-STR diversity, forensic parameters, marker-panel comparisons, sensitivity analyses, population structure, spatial analyses, and exploratory haplogroup reporting. The repository is an RStudio and command-line workflow, rather than an installable R package.

The repository contains analysis code and synthetic example data. Individual-level study data are not included. The synthetic examples demonstrate the workflow and do not reproduce the study’s numerical results.

## Getting started

Open `maldives-ystr-analysis.Rproj` in RStudio. Install `renv` if needed, then restore the recorded dependencies:

```r
renv::restore()
```

Run the synthetic example from the repository root:

```r
source("R/pipeline.R")
run_pipeline("config/analysis_settings.yml")
```

Alternatively, use the command line:

```sh
Rscript run_pipeline.R --config config/analysis_settings.yml
```

After the workflow completes, run the automated checks:

```sh
Rscript tests/run_tests.R
```

The example uses reduced iteration counts and image resolution to shorten execution time. Publication assembly is intentionally disabled for synthetic runs. To repeat an analysis, use a new run ID or resume a valid checkpoint.

## Using your own data

Follow the [data dictionary](docs/data_dictionary.md) and copy the [configuration template](config/private_settings.template.yml) to `config/private/settings.yml`. Supply your authorized input data, quality-control evidence, geographic metadata, and predictor outputs in the designated local folders.

The template records the study’s methodological settings. The synthetic example uses a separate configuration and seed table. See [study configuration](docs/study_configuration.md) for iteration counts, random seeds, DAPC settings, and input-dependent PC selection.

Files placed in the designated private-input folders are excluded from version control by the supplied ignore rules. Review the files staged for commit before sharing a modified repository.

Haplogroup predictor outputs must be supplied separately; the workflow does not submit profiles to an external prediction service.

## Arlequin analyses

Primary AMOVA and pairwise RST analyses are performed separately in the Arlequin GUI. The R workflow prepares input projects and imports the resulting output files.

Configure the expected result paths before starting, then run through the export stage:

```r
source("R/pipeline.R")
run_pipeline(
  "config/private/settings.yml",
  until = "09b_export_arlequin"
)
```

Run the exported projects in Arlequin and save the original output files at the configured locations. Resume the workflow with:

```r
run_pipeline(
  "config/private/settings.yml",
  from = "09c_import_arlequin",
  resume = TRUE
)
```

The synthetic example uses separately labeled R calculations to test downstream processing. These do not replace the primary Arlequin analyses. See the [Arlequin instructions](docs/external_arlequin.md).

## Analytical conventions

- Single observed alleles at multicopy systems are retained as single analytical states. Missing copies are not inferred or duplicated; observed multicopy patterns are treated as unordered.
- Forensic haplotypes retain DYS389II as typed. Population-structure analyses use DYS389II minus DYS389I.
- Allele-size-based analyses exclude profiles containing microvariants at included single-copy loci before the DYS389II adjustment. Alleles are not truncated or rounded. Eligibility is reassessed separately for each marker panel.
- Categorical DAPC retains microvariants as distinct allele states in the complete-profile dataset.
- The primary random-removal analysis removes complete profiles without replacement, matching the administrative-division distribution of the excluded incomplete profiles. Incomplete profiles remain in each simulated dataset.
- DAPC preprocessing is fitted to training profiles only. Exploratory predicted-haplogroup validation includes categories represented by at least two profiles; singleton categories remain in training only. A separate sensitivity analysis uses categories represented by at least five profiles.
- Haplogroup predictions are exploratory and do not constitute diagnostic Y-SNP confirmation.
- Sample sizes, effective PC grids, selected PCs, and analysis results are calculated from the supplied inputs.

## Outputs and reproducibility

The workflow records the effective configuration, random seeds, random-number generator settings, package versions, and stage logs. Input and source fingerprints protect against resuming analyses after their inputs or code have changed.

Figures are exported as PDF, SVG, PNG, and TIFF in color and grayscale. PDF and SVG provide vector output; raster resolution is configurable. The study configuration specifies 600 dpi.

Generated outputs, local data, external Arlequin results, and installed libraries are excluded from the source archive. Sampling-map assets must be supplied separately.

Before modifying an extracted release, check its source integrity:

```sh
Rscript tools/check_source_integrity.R
```

The automated tests use synthetic data and independently specified test examples. Passing these tests verifies workflow behavior; reproducing the study’s numerical results additionally requires the corresponding study inputs, external results, and documented analysis environment.

## Repository structure

| Folder | Contents |
|---|---|
| `R/` | Analysis functions and figure/table routines |
| `scripts/` | Ordered workflow entry points |
| `config/` | Marker panels, settings, and seed tables |
| `data/example_synthetic/` | Synthetic example inputs |
| `verification/` | Additional verification and sensitivity scripts |
| `docs/` | Methods, data requirements, execution instructions, and testing information |

## Documentation

- [Data dictionary](docs/data_dictionary.md)
- [Execution guide](docs/workflow.md)
- [Methods](docs/methods.md)
- [Study configuration](docs/study_configuration.md)
- [Arlequin instructions](docs/external_arlequin.md)
- [Software environment](docs/environment.md)
- [Test report](docs/TEST_REPORT.md)
- [Change history](CHANGELOG.md)

## License and citation

The code is distributed under the GNU General Public License v3.0 only (`GPL-3.0-only`). See `LICENSE` and `NOTICE.md`.

Citation metadata are provided in `CITATION.cff`. Please cite the workflow version used and the associated study when applicable.

Repository: https://github.com/yamyyn/maldives-ystr-analysis

Versioned source: https://github.com/yamyyn/maldives-ystr-analysis/tree/v1.2.1

No software DOI has been assigned.
