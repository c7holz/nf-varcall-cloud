# 0004 — Reference index: expect pre-staged, build only on explicit opt-in

- Status: Accepted
- Date: 2026-06-29

## Context
BWA_MEM requires a prebuilt index of the reference FASTA. The index can be produced
inside the pipeline or supplied as a pre-staged input. Building it on every run wastes
compute and, on Batch, rewrites the index into the S3 work-dir each time. Building it
silently when absent hides an expensive step behind a default. Requiring a pre-staged
index unconditionally would block a first run where none exists yet.

## Decision
Expect a pre-staged index by default (`--bwa_index <dir>`). Build it only when explicitly
requested via `--build_index`; the build step publishes the index to the output directory
so it can be staged as input on later runs. If neither `--build_index` is set nor
`--bwa_index` is given, the pipeline fails immediately with a message pointing at
`--build_index`.

## Consequences
- No accidental rebuilds: the expensive step runs only when explicitly requested.
- The index is an explicit, retrievable artifact (published to the output dir), not an
  implicitly cached one — the operator stages it back on later runs. Predictable and
  auditable, at the cost of one manual hand-off.
- Fail-fast with guidance keeps the default safe: a missing index stops the run early
  with a clear next step instead of failing deep in alignment.
- Rejected: unconditional in-pipeline build with `storeDir` caching. Simpler for the user
  but couples reference provenance to a work-dir cache and hides a costly step behind a
  default — weaker for reproducibility and cost control.
