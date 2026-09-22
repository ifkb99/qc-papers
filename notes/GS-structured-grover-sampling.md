---
code: GS
title: "Exact structured Grover output sampling through 128 qubits"
outcome: record
claims: [C117]
todo: [71]
---
# GS — From selected probabilities to complete samples (2026-09-21)

The user asked whether C114's structured Grover memory remained exponential,
what 128 qubits would cost, and to start the proposed sampling computation.
The accepted prior authorization for parallel agents remained in force.
One independent author audited scaling and primary sources while another
derived the sampler. The coordinator proposed investigating a reversible
homogeneous suffix recurrence; the sampler author proved the explicit
Fibonacci environment and retired quadratic recomputation before testing.
C117 owns the scientific results and measured scope; TX38 owns the transfer.

The useful distinction was fixed iterations versus an amplification-length
schedule, and a compact state versus an explicitly returned probability vector.
The 128-qubit final-H distribution at three iterations is extremely concentrated
at zero. Thus ordinary samples alone would be a vacuous correctness test.
Exact rare-prefix formulas became a mandatory check before resource runs.

## Evidence and independence

| Work | Submission | Review | Author directory |
|---|---|---|---|
| Scaling/source audit | S624b8a7794754d79 | V1f6b5670d61a4a81 | out/agent-board/workers/A352551992a344206/ |
| Sampler derivation/design | Sdb1dac2311534d02 | Vcaf3cc1759ff4737 | out/agent-board/workers/A4f3325ebc5bf49f7/ |
| Execution | S76be296197df4469 | V0039c2164aaa4b59 | out/agent-board/workers/A6d1656f712ec4b6b/ |

The design referee was initially fresh, a non-author, and reviewed scaling
before the sampler design; this context reuse was disclosed. A separate fresh
results referee inspected the execution. The existing dense reference received
a reviewed direct-reference exemption: it builds the predicate clauses, applies
phase-oracle and mean-reflection updates, and uses the existing Walsh transform.
It is algorithmically independent of the sampler recurrence, but its author
had seen it, so the reference was not blind. C114 scalar checks share algebra
and are supplementary rather than independent validation. The reference report
was recorded and frozen before candidate validation, and the successful
candidate report before resource cases.

The audit's first arithmetic script printed an unconditional success line after
assertions. The shared lint flagged that display; a preserved second version
removed it, and a third added the requested exact nonzero probability. A lint
CLI misuse and the design author's 44-versus46 check-count transcription were
also preserved and explicitly corrected. These were tooling/reporting issues,
not silently replaced scientific outcomes. The frozen submissions record them.

## Reproduction and measurement contract

From research/, with fresh output paths:

```sh
uv run python test_core.py
export LAB_GPU=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p out/grover_sampling_repro
flock out/grover_sampling_repro/science.lock timeout 120s uv run python -m experiments.experiment_grover_sampling --stage reference --report out/grover_sampling_repro/reference.json
flock out/grover_sampling_repro/science.lock timeout 120s uv run python -m experiments.experiment_grover_sampling --stage validate --reference out/grover_sampling_repro/reference.json --preconditions experiments/grover_sampling_preconditions.json --report out/grover_sampling_repro/validation.json
```

Stop after any nonzero exit. These commands reproduce correctness checks,
not the resource timings. The runner's case stage accepts --n, --method
rolling|cached, --mode time|allocation, --validation and --report. A resource
replication must preserve the recorded 24-case grid, separate instruments and
accepted 1GiB sampled-RSS/120-second guard; see the frozen resource_batch_v1.py
and resource_bundle_v1.json in the execution directory. The reusable fixture
is a byte-for-byte copy of the accepted design precondition report.

Resource cases use lightweight standard-library imports; dense/reference stages
also import NumPy, Walsh and the lab harness. Whole-process RSS is therefore
not directly comparable with the earlier query pilot's eager common imports.
The memory claim concerns the measured call and the stated process, not a
ratio obtained by mixing those import environments. Calls include construction,
coefficient and suffix work, RNG and retained Python integer outputs. Imports,
launch, lock acquisition and report serialization are outside timed/traced
calls; raw process RSS includes interpreter/import overhead. All-zero integer
outputs may share objects and are not fixed-width strings in memory.

The fixed resource case order was not randomized; timing values are exploratory
host measurements. A derived dense128 payload is not an empirical benchmark,
and explicitly stored exact rationals have different precision from real64 or
complex128 arrays. The strongest baseline here is the same structured
contraction, not exhaustive enumeration. No paper or generic oracle claim is
promoted. Raw evidence under out/ remains gitignored; no commit was requested.

Canonical integration is Tc85b9015e04b40ac. The helper module docstring was
updated and the result review's inaccurate quantile-extremes comment corrected;
no executable code changed. Exact allowed-text transformations and AST identity,
including an executable-mutation detector, accompany the integration evidence.
The reference and validation stages run again under canonical package imports;
resource numbers are reused from the frozen accepted execution, not remeasured.
The coordinator authored integration prose only. A fresh integration referee
checks it separately. The original dirty worktree is preserved outside these
new records, TODO closure and generated indexes/handoff.

Integration evidence lint caught check objects labeled with name instead of id.
The mechanical checker was versioned with that schema correction, preserving
its original reports and failed lint log. All conditions were unchanged; the
scientific source and resource results were not rerun for this format repair.
