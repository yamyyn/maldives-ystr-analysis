# Private audit evidence

Use `config/original_audit_plan.csv` and copy the blank targets template into ignored `data/private/`. Targets use `key,comparison,expected,tolerance`; observations use `key,value,evidence` with fresh independent measurements. Supported comparisons: equal, near, file_exists, at_least, at_most.

Review legacy directory checks against documented modern equivalents rather than creating empty folders to manufacture PASS. Configure private target/observation paths before a run. Missing evidence remains incomplete. Historical statuses are not inherited. The nine optional supplemental checks do not replace the 119 original definitions. See `workflow.md` and `methods.md` for scope and source boundaries.
