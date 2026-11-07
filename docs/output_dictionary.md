# Output dictionary

| Location | Meaning |
|---|---|
| `data/clean/<run_id>/full/Full_Archival_Unchanged.*` | Preserved import; not modified for analysis |
| `data/clean/<run_id>/full/Forensic_Full_Analysis_Ready.*` | Confirmed derived no-call coding and observed multicopy states |
| `data/clean/<run_id>/complete/` | Complete analysis cohort, before structure-specific transformations |
| `data/clean/<run_id>/deletion/` | Explicitly confirmed deletion cohort with source calls retained |
| `data/clean/<run_id>/review/` | Partial/invalid profile details and independent eligibility flags |
| `results/<run_id>/qc/project_info/` | Import, metadata, genotype and cohort-creation workbooks |
| `results/<run_id>/tables/<stage>/` | Calculated stage tables and workbook exports |
| `results/<run_id>/intermediate/<stage>.rds` | Detailed computed results, raw iterations, matrices and model objects |
| `results/<run_id>/publication_outputs/tables/` | Source-derived main/supplementary writers plus inference/CV tables |
| `results/<run_id>/publication_outputs/figures/` | Selected final source-derived figures, color/grayscale PDF/PNG/TIFF |
| `results/<run_id>/figures/additional/` | Original allele-frequency plates and new supporting diagnostic plots |
| `results/<run_id>/captions/Caption_Review.csv` | Draft captions requiring author review |
| `results/<run_id>/qc/Publication_Gate.csv` | Requirements for assembling a locked study package |
| `results/<run_id>/stage_log.csv` | Completed stages, elapsed time and recorded warnings |
| `results/<run_id>/pipeline_state.rds` | Resume checkpoint tied to input/config/code fingerprint |
| `external/arlequin/<run_id>/` | Exported projects, frozen label maps and recoding audit |

All generated study outputs are private by default. The downloadable repository contains only code, documentation, settings, blank validation templates and invented example input. Optional synthetic output previews are labelled separately. Output numbering and historical counterparts are recorded in `config/items.csv` and `docs/output_traceability.csv`.
