# nf-varcall-cloud

A tested Nextflow pipeline (FASTQ -> variants) that runs identically on a laptop with containers and elastically on AWS Batch with S3 as the work directory.

The point is not the analysis (standard) but the engineering: local<->cloud parity, reproducibility, IaC-readiness, and cost-aware cloud execution.

Architecture & Decisions -> [`docs/architecture.md`](docs/architecture.md)
