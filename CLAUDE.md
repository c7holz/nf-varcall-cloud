# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow pipeline (FASTQ → variants) whose point is the **engineering, not the analysis**:
the same DAG runs locally with containers and on AWS Batch with S3 as the work-dir. The three
quality goals that drive every decision are **local↔cloud parity**, **reproducibility** (pinned
containers, identical results), and **cost-awareness** (Spot, teardown, cheap test region).

**Current state: skeleton.** Every process `script`/`stub` block only `touch`es its declared
outputs — no real bioinformatics tools run yet. Tests exercise DAG wiring via `-stub-run`, not
tool behavior. When implementing a real tool, replace the `script` block *and* add a pinned
`container` directive (none exist yet); keep the `stub` block so stub tests stay fast.

## Commands

```bash
nf-test test                              # run all tests (stub-run; what CI runs)
nf-test test tests/main.nf.test           # run one test file
nf-test test --filter '<name substring>'  # run a single named test

# Run the pipeline directly (requires --input samplesheet and --fasta):
nextflow run main.nf -profile docker,local \
  --input tests/data/samplesheet.csv --fasta tests/data/ref.fasta --build_index
nextflow run main.nf -profile docker,local --stub-run --input ... --fasta ...
```

CI (`.github/workflows/ci.yml`) runs `nf-test test` on every push and PR — keep it green.

## Architecture

- `main.nf` — the workflow DAG: FASTQC + FASTP → BWA_MEM → SAMTOOLS → BCFTOOLS, plus MULTIQC
  aggregating FASTQC/FASTP reports. This file is **executor-agnostic**; nothing here knows
  about local vs. cloud.
- `modules/*.nf` — one process each, DSL2. Convention: `tuple val(meta), path(...)` channels
  where `meta` carries `[id: sample]`; named outputs via `emit:`.
- `nextflow.config` — defines `params` and three **profiles**: `local` (Docker), `docker`,
  `awsbatch`. `conf/base.config` is always included (baseline cpus/memory); each profile pulls
  in its `conf/*.config`. Parity is achieved by keeping the DAG fixed and letting **only the
  profile** swap executor + work-dir. `conf/awsbatch.config` is still a placeholder.
- Reference index contract: the pipeline **expects a pre-staged BWA index** (`--bwa_index <dir>`)
  and fails fast if absent; building is opt-in via `--build_index` (see ADR 0004). Preserve this
  fail-fast behavior — a test asserts on the `"No BWA index provided"` message.

## Conventions

- **Decisions with a rejected alternative → an ADR** in `docs/adr/` (Nygard format); update the
  table in `docs/adr/README.md`. Stable big-picture lives in `docs/architecture.md` (arc42-light).
  A note without a rejected alternative is not an ADR.
- Samplesheet is CSV with header `sample,fastq_1,fastq_2`.
- Test data lives in `tests/data/`; keep test inputs tiny (small region / chr20) to keep runtime
  and cost minimal. Public data only — never patient data.
