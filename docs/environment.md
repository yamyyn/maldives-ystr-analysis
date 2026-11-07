# Environment

`renv.lock` records actual installed dependency versions used in local testing. `.Rprofile` and `renv/activate.R` provide standard activation; restore with `renv::restore()` from the project. R itself, system fonts and external Arlequin are separate requirements.

Local testing used Windows R 4.6.1. The lockfile was snapshotted from installed libraries. A full network restore into an empty library on another operating system is a separate unverified step. Retrieval may require archive access or compilation when matching binaries are unavailable.

Arial is the default font; substitutes can change line breaks. Cairo PDF/PNG/TIFF and LZW TIFF are used where supported. Snapshot intentional environment changes with `Rscript tools/snapshot_environment.R`, then repeat tests. Historical reports are reference evidence, not proof of a restored current environment.
