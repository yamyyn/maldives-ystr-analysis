# Workflow and execution

Run every command from the repository root. Use a fresh run ID for changed data, code or settings. Inputs and code are hashed; resume rejects changed inputs rather than mixing runs.

```text
Rscript --vanilla run_pipeline.R --list
Rscript --vanilla run_pipeline.R --config config/analysis_settings.yml
Rscript --vanilla tests/run_tests.R
```

The default is independently invented synthetic data, reduced iteration counts and real computation of the analysis stages. No cached ordination, invented spatial result or placeholder DAPC coordinate is used by this runner. The synthetic structure calculation is an explicitly labelled R method-of-moments verification, not an Arlequin run. It does not certify private manuscript results.

## Stage sequence

00 setup/environment; 01 character import and QC; 02 archival/full/complete/deletion preparation; 03 sampling summaries; 04 allele frequencies/locus diversity; 05 full and complete forensic statistics; 06 revised stratified and unstratified random-removal sensitivity plus QC-only paired bootstrap; 07 panels and leave-one-system-out contribution; 08 accumulation; 09a adjusted structure cohorts; 09b Arlequin export; 09c external import/ordination; 10 spatial Mantel/correlogram; 11 geographic DAPC/CV; 12 predicted-parent processing and DAPC/sensitivity; 13 rarefaction/sharing; 14 internal bootstrap; 15 manuscript tables; 16 manuscript and supporting figures; 17 final QC/traceability; 18 gated publication assembly.

Individual stage scripts invoke the same runner. A stage requires preceding checkpoints; it cannot silently substitute demonstration results for missing inputs.

## Private configuration

1. Copy `config/private_settings.template.yml` to ignored `config/private/settings.yml` and the QC template to `config/private/initial_qc.R`.
2. Configure private CSV input, assay slots, independent population mapping, exact confirmed no-call patterns, and their evidence. Add `classification=confirmed_deletion` only where independently established. Never source synthetic confirmation records for private data.
3. Supply the leaf posterior, leaf-to-parent mapping and independently reported parent labels described in `data_dictionary.md`. Modeled classes include zeros. The reported parent is preserved; it is not silently replaced by a new argmax call.
4. Choose and document iteration counts, PC grid and random seeds. The private template contains methodological starting values, not proof that they reproduce every historical setting; compare them to the selected original/revised analyses before inference.
5. Configure the three future external result paths before the run. These files may be absent until the Arlequin checkpoint.

```text
Rscript --vanilla run_pipeline.R --config config/private/settings.yml --until 09b_export_arlequin
```

Run Arlequin externally as described in `external_arlequin.md`. Put results in the already configured paths, then:

```text
Rscript --vanilla run_pipeline.R --config config/private/settings.yml --resume
```

Keep the config unchanged for resume. If settings/code/input change, start a new run ID. External results must match frozen cohort labels and settings. Record software version and preserve the raw external result files.

## Output locations

`data/clean/<run_id>/` contains the archival, normalized-for-review, full forensic, complete, confirmed-deletion and review/provenance exports. `results/<run_id>/` contains stage tables, serialized detailed results and iterations, logs, QC, manuscript outputs, supporting figures and captions. `external/arlequin/<run_id>/` contains projects and label mappings. All generated/private contents are excluded from the public ZIP and Git tracking by default.

## Final status

A completed synthetic run must end with `PUBLICATION_PENDING.txt`, not a locked study package. A private package requires complete independently checked original audit targets/observations, approved map files in all requested formats, approved captions, verified author/license metadata and a signed-off private numerical/visual comparison. These approvals are explicit configuration flags; the software does not infer them from a successful run. The pipeline assembles files locally and does not publish or upload them.
