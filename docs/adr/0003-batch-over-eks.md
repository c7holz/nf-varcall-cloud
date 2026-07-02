# 0003 — AWS Batch over EKS

- Status: Accepted
- Date: 2026-06-29

## Context
Nextflow needs a cloud compute backend to dispatch jobs. Two managed options: AWS Batch
(managed job queues + compute environments, Nextflow `awsbatch` executor) or EKS (managed
Kubernetes, Nextflow `k8s` executor).

## Decision
Use AWS Batch. Jobs are dispatched to a Batch job queue backed by a compute environment
with Spot instances.

## Consequences
- Maps cleanly onto a scheduler + partition mental model: a job queue plus compute
  environment behaves like a batch scheduler with a job partition. The new operational
  surface stays small.
- Managed scaling and native Spot support in the compute environment serve the cost NFR
  with little configuration.
- No cluster to operate — no control plane, node groups, or k8s objects to maintain.
  Appropriate for a single-pipeline scope.
- EKS would give finer scheduling control and more cross-cloud portability, but carries
  standing control-plane cost and Kubernetes operational overhead that this scope does not
  justify. Would be reconsidered if requirements change — e.g. multi-cloud portability,
  shared infrastructure with non-Nextflow workloads, or scheduling needs Batch can't meet.
