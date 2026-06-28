# 0001 — Keep Nextflow as orchestrator, swap the executor

- Status: Accepted
- Date: 2026-06-28

## Context
The existing pipeline experience is in Nextflow (nf-test, modules, CI). The goal is
to run the *same* pipeline locally and on AWS Batch. Two broad options: (a) keep
Nextflow and change only the executor + work-dir, or (b) rewrite against a different
cloud-native workflow engine / direct Batch job definitions.

## Decision
Keep Nextflow as the orchestrator. Move from a local executor to the `awsbatch`
executor, and from a shared filesystem to S3 as the work-dir. Pipeline logic stays
unchanged; only the execution backend and data staging change.

## Consequences
- Local↔cloud parity becomes the core demonstrable property: one codebase, two profiles.
- Reuses existing Nextflow/nf-test/CI investment; the new learning surface narrows to
  Batch + S3 — which is the actual point of the project.
- Couples the project to Nextflow's cloud-execution model (S3/Batch integration). A
  future move to a different engine would be a larger change. Accepted, since Nextflow
  is the established tool in the target domain (nf-core).
