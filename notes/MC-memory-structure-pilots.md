---
code: MC
title: "Memory pilots: structured Grover queries, clean-block sampling and modular branch sharing"
outcome: record
claims: [C114, C115, C116]
todo: [66, 71]
---
# MC — Exploiting structure without hiding evaluation cost (2026-09-21)

The user asked for lower simulator memory, accepting moderate extra runtime,
and authorized parallel tests. Three authors investigated clean arithmetic
blocks, structured Grover, and modular branch sharing. A provisional 2x runtime
ceiling guided engineering comparisons; it was not a universal user requirement.
The slate used construction inventory, weighted counting and an interference
obstruction. TX35–TX37 own the transfers and their hypotheses.

## What the investigation changed

[C114](../claims/C114.md) owns the structured Grover derivation and measurements.
The useful pattern is bounded-width weighted counting: it keeps interference
and charges count construction. This does not improve on the strongest
structured classical method, which is the same counter. A small-solution
representation was rejected as the explanation for this exponentially large
marked set; TX35 records the primary-source boundary.

[C115](../claims/C115.md) owns the clean-block pilot. Its result reverses the
initial engineering preference for this branch: certification cost and existing
logical baselines determine practical memory, not scratch-wire count or the
number of amplitude buffers. The code remains an experiment. TODO 66 keeps
the separate work-domain-closure question open; TODO 64 is not solved.

[C116](../claims/C116.md) owns the modular-sharing identity and collision
limitation. Arbitrary modular labels are more expressive than tensor-product
Pauli labels, but their compact notation does not evaluate coherent sums.
The optional order-discovery benchmark was retired before performance runs
after review exposed an incomplete baseline set and contaminated RSS plan.
No unmeasured win or domination is inferred.

## Evidence, review and reproduction

| Branch | Accepted submission | Review | Frozen author directory |
|---|---|---|---|
| Clean-block design | S201e0d2f9f8143ad | V7b2c6529bcbe41d5 | out/agent-board/workers/A62e7672e917644f4/ |
| Clean-block result | Sc26366f437c64436 | V63b15412bfc14257 | out/agent-board/workers/Abb5d0b5e50f146c8/ |
| Grover design | S3e65b1839e1249c7 | Vebc23ee2300242b8 | out/agent-board/workers/Afeed8b6ed1b5458c/ |
| Grover result | Sce91c456c4f9491f | V6d8d445de2cf48b2 | out/agent-board/workers/A8f0f133336a24990/ |
| Modular sharing | S11674161eb334aa0 | V0b1b235da1504034 | out/agent-board/workers/A7aad0be0fc144816/ |

Clean-block sources/results: clean_block_v1.py, experiment_v1.py,
resource_summary_v1.json and reports_exact_bundle_v1.json. The bundle contains
all reports, finished run records and sealed logs, including those beyond the
board's per-submission list limit. Sharing sources: derivation.md,
disposition_v2.md and reports_bundle_v2.json. The disposition supersedes stale
execution/status text in the preserved first derivation and retired pilot.
Grover final evidence: results_v2.md, comparison_v2.json,
result_bundle_v2.json, guard_frozen_bundle_v2.json and science_identity_v2.json.
Canonical integration is Tf77dcd85e9d54484; exact source-equivalence detectors,
validation logs and reviewer-prompted statement provenance accompany its
submission. Raw board evidence under out/ is gitignored; the canonical helper
and experiment sources below are retained in the worktree.

For a fresh bounded validation run from research/, after the core gate:

```sh
uv run python test_core.py
export LAB_GPU=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out/memory_repro
flock out/memory_repro/science.lock timeout 120s uv run python -m experiments.experiment_clean_block validate --reviewed --report out/memory_repro/clean_validation.json
flock out/memory_repro/science.lock timeout 60s uv run python -m experiments.experiment_grover_queries --stage validate --report out/memory_repro/grover_validation.json
```

Use fresh report paths; do not continue after a nonzero exit. To inspect a
single resource case, the clean runner accepts measure --reviewed --backend
candidate|logical-dense|logical-sparse --a 2|3 --width 6 --samples 128
--instrument timing|allocation --report PATH; the Grover runner accepts
--stage case --method compact|streaming|dense --mode time|allocation --n 12..20
--report PATH. Those are argument choices, not literal shell commands. Preserve
separate processes/instruments and the accepted RSS guards for a resource
replication. The archived case manifests specify the full sweep/repetitions;
validation commands alone do not reproduce its timing numbers.

The host permitted three worker contexts but rejected additional fresh-context
spawns. Non-authors cross-reviewed the designs/results, with context reuse
recorded explicitly. These were not blind or fresh-context reviews. References
used existing/direct methods under reviewed exemptions, with shared constructors
and libraries disclosed. Every science run followed the core gate. Cases used
single-thread CPU execution, sampled RSS guards and serialized measurements;
parallel authorship did not mean overlapping timed workloads.

An initial sharing diagnostic stopped on an address-space cap because eager
imports reserve far more virtual space than resident memory. A reviewed sampled
RSS guard replaced that cap; the unchanged six-check diagnostic then passed.
Grover's first resource batch stopped on a /proc process-exit race. Its failed
case and late output were excluded, seven successful cases preserved by hash,
and the guard-only repair independently reviewed before continuation. Original
failures and superseded evidence remain in the board archive.

New scoped records C114–C116 own the derived/measured findings; no existing
claim or paper is promoted, no preferred backend is changed, and no generic
Shor/Grover speedup is asserted. The next proposed step is the exact structured
sampler in TODO 71; it needs its own derivation and comparison contract.
