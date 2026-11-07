# Maldives Y-STR analysis workflow

Version 1.1.1. Reusable R analysis scripts for Y-STR diversity, forensic parameters, marker panels, sensitivity analyses, population structure, spatial analyses, and exploratory haplogroup reporting. This is an RStudio/command-line workflow repository, not an installable R library package.

Only explicitly invented demonstration data are included. Private genetic profiles, participant identifiers, predictor records, manuscript results and private verification targets are excluded. Results are calculated from the user's input; demonstration output does not reproduce the private study.

## Quick start

Open `maldives-ystr-analysis.Rproj`. Install `renv` explicitly if needed, then run `renv::restore()` to restore the recorded dependency versions. No analysis script installs or updates packages automatically.

From the repository root:

```r
source("R/pipeline.R")
run_pipeline("config/analysis_settings.yml")
```

Or use `Rscript run_pipeline.R --config config/analysis_settings.yml`. After the synthetic pipeline completes, run `Rscript tests/run_tests.R`. The demonstration uses reduced iteration counts and raster resolution for speed. Its publication gate intentionally remains withheld. Run in a fresh extracted copy; repeated runs require a new run ID or a valid resume checkpoint.

## Own-data analysis and Arlequin GUI

Follow [the data dictionary](docs/data_dictionary.md) and [private configuration template](config/private_settings.template.yml). Supply authorized input and predictor outputs locally in ignored folders. No haplogroup predictor service is called by this repository.

```r
source("R/pipeline.R")
run_pipeline("config/private/settings.yml", until="09b_export_arlequin")
# Run the exported projects separately in the Arlequin GUI.
run_pipeline("config/private/settings.yml", from="09c_import_arlequin", resume=TRUE)
```

Configure result paths before starting. The R workflow exports and imports Arlequin projects/results; it does not run the GUI or replace primary AMOVA with R estimates. The synthetic-only R branch and standalone point-estimate verification are explicitly separate. See [Arlequin instructions](docs/external_arlequin.md).

## Scientific conventions

- Single observed multicopy allele states are not duplicated or interpreted as molecular copy-number evidence. Raw calls are preserved; unordered analytical patterns are used consistently.
- Forensic haplotypes retain typed DYS389II. Structure analysis uses DYS389II minus DYS389I, with separate audited integer coding for Arlequin and categorical microvariants for DAPC.
- Primary random removal samples complete profiles only, without replacement, matching the excluded-profile division counts; the incomplete profiles remain in each simulated dataset.
- DAPC fits preprocessing on training profiles only. Haplogroup cross-validation includes categories with at least two profiles; singleton categories remain in training. The separate rare-group analysis requires at least five profiles per category. These thresholds are configurable.
- No sample size, observed result, model-selected PC count or manuscript P-value is used as a fixed verification target.
- PDF and SVG are vector exports; PNG and TIFF use configured resolution. The private template specifies 600 dpi. Existing figure exports are separate assets and are not bundled here.

## Repository layout

`scripts/` contains ordered entry points; `R/` contains reusable methods and source-derived recipes; `config/` defines marker panels and methodological settings; `data/example_synthetic/` contains invented input. `verification/` contains separate sensitivity checks. Generated results, private inputs, external results and installed libraries are ignored and omitted from the release archive.

See [methods](docs/methods.md), [execution guide](docs/workflow.md), [release notes](CHANGELOG.md), [test report](docs/TEST_REPORT.md), and [environment](docs/environment.md). Historical source/output mappings document provenance; they do not certify private numerical equality or exact publication artwork. Sampling-map assets remain external author-supplied inputs.

## Licence and citation

Code is distributed under GPL-3.0-only, selected by the study author. See `LICENSE`, `NOTICE.md` and `CITATION.cff`. The source repository is [https://github.com/yamyyn/maldives-ystr-analysis](https://github.com/yamyyn/maldives-ystr-analysis). No software DOI has been assigned.

Before making local source edits, verify the extracted release with Rscript tools/check_source_integrity.R. The manifest covers source integrity, not private-data equivalence. Test fixtures contain independent toy arithmetic expectations, never fixed manuscript results.
