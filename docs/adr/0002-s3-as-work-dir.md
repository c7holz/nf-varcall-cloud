# 0002 — S3 as the work-dir

- Status: Accepted
- Date: 2026-06-29

## Context
On AWS Batch, every job needs a work-dir reachable across compute nodes. On-prem this
was the shared cluster filesystem. In the cloud the options are a managed shared POSIX
filesystem (EFS, or FSx for Lustre) mounted into the compute environment, or S3 as the
work-dir via Nextflow's native S3 support.

## Decision
Use S3 as the work-dir. No shared filesystem is provisioned; Nextflow stages inputs and
outputs to and from S3.

## Consequences
- Nothing to provision, mount, scale, or tear down between runs — fits the
  teardown-on-idle cost model and the cost-transparency NFR.
- Native `-resume` works against S3; staging is handled by Nextflow.
- Higher per-task I/O latency than a Lustre filesystem. Irrelevant at this
  region-limited scope; would be revisited for large high-throughput workloads.
- FSx for Lustre would be faster but adds standing cost and provisioning complexity,
  against the cost NFR. Rejected for this scope.
