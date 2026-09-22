---
code: FE
date: 2026-09-21
title: "Factoring follow-through: cache lifetimes and actual giant-step ratios"
outcome: record
claims: [C120, C121]
todo: [73]
---
# FE — Follow-through on the factoring slate

The user approved FD's recommendation: audit practical cache-construction memory
and actual arithmetic giant-step recurrences in parallel, then measure a useful
cache change. C120 owns the measured outcome and limits; C121 owns the exact
symbolic result and missing lemma. TODO73 owns the residual decisions. This
round did not revisit the preserved factor-state constructor or run challenge
factoring, Shor simulation or the speculative paired-hyperbola proposal.

## Authorship and sequence

Separate workers authored the CADO source audit/design and the arithmetic
recurrence audit. A fresh referee accepted both, with exposure to the first
review disclosed for the second. The cache author retained responsibility for
implementation and measurement; the coordinator integrated only supported
results. The existing upstream cache was the exact-format reference and direct
sparse XOR the independent semantic reference. Their shared fixture generation
and upstream loader/backend limitations remain explicit in C120; no independent
complete CADO implementation was claimed.

| Work | Frozen submission | Review |
|---|---|---|
| Cache ownership audit and experiment design | S7fd475c348924f82 | V9eeecdb5c24e4253 |
| Exact recurrence audit | S9aed30d3c93b421d | Vee36cf7331cd4f55 |
| Initial actual-CADO implementation and observations | S3d85e0672aa747a2 | Vf4eac155b56c4dd5 (timing revision required) |
| Corrected same-workload timing | S5abdf8a4b01c4834 | V06175d0782df4152 |

The core science gate passed before scientific execution. Cache validation run
R5fc0653ee9f34f90 preceded performance R85217ea3166749c8; both exited zero.
Corrected run Rb23e784ec68644ac also exited zero. Every run had its prediction
and board record before execution. The recurrence audit used hand algebra and
needed no numerical sweep. Its counterexample rejects one proposed shortcut,
not all short covers or factoring algorithms.

## Timer correction and retained limits

Author and coordinator independently noticed coarse clusters in the first
elapsed readings. Python's timeout-mode process wait polls with a sleep ceiling;
the fresh referee checked that source mechanism. Existing RSS and rounded CPU
readings remained valid, but their accompanying launcher intervals could not
resolve a modest slowdown. The completed runs were preserved with exit zero;
they were not retrospectively called failed arithmetic gates.

The focused revision replaced polling with blocking wait and a separate watchdog,
kept the same compiled programs, matrices, allocator and exact-output checks,
and reran only the prescribed sequence. The final report retains every row,
including the slower combined-arm default-policy observation. It does not pool
old coarse elapsed values with corrected measurements. A fresh second context
reviewed amendment conformance and the corrected raw records.

Core/validation, builds, fixture generation, hashing and reporting are outside
each timed builder process. All are separately identified in the frozen reports.
The initial configure failure was an uninstalled build dependency, resolved in
an attempt-local pinned environment before science. Prior interpreter setup CPU
was not measured; the correction's prior-work allowance is administrative,
not a certified bound. Child caps and the accounting ledger were respected.
The corrected watchdog's timeout branch was not exercised by these fast runs.

## Reproducibility and integration

The portable package is experiments/cado_cache/README.md, with the exact patches,
build scripts, validation driver, corrected timer and explicit input-manifest
helpers. Canonical measurements_v2.json retains all corrected observations.
The unchanged executable/fixture reports, individual logs and frozen bundles
remain local under the two experiment attempt directories and in arb archives;
they are gitignored. No rebuild or fresh reproduction at the canonical path was
performed merely for integration; installed executable code matches reviewed
bytes. The README separates replay from a new environment with new hashes.

C121 preserves the accepted worker derivation; its source is Harvey--Hittmeir's
Algorithm 4.3 and exact symbolic arithmetic, not a referee's wording. C120's
table comes from the frozen corrected measurements. The integration report
maps every review-related scope correction back to code, report, raw row or
algebra. TX44/TX45 cross-references are updated only after scientific integration
closes, preserving its transitive input identities. This is administrative
sequencing, not a replacement of the accepted scientific evidence.

The useful practical outcome is a smaller construction peak on the tested CADO
path. It does not interleave into an extra multiplicative speedup with C121;
that branch supplies an exact diagnosis and an unproved constructive lemma.
No commit, publication or upstream patch submission is included.
