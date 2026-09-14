---
code: UG
date: 2026-09-11
title: "A uniform mixing certificate and a more useful exact checkpoint baseline"
outcome: confirmed
claims: [C64, C65, C66, C67]
todo: [14, 24, 27, 28]
---
# UG — Quantify the gap, then charge a stronger baseline

C66 owns the uniform phase-gap proof and approximation contract; C67 owns
the exact checkpoint recurrence and its cost boundary. This note records
TODO 27's bounded evidence. The qsim-research skill's stronger-baseline rule
prompted the checkpoint implementation because the sufficient mixing warm-up
from §FM was already outside the current width limit. Lower-cost agents
supplied the initial probes and algebra audits; main implemented the opt-in
production change, strengthened verifier predicates and reran the evidence.

## Uniform certificate, not a useful approximation warm-up

`experiment_uniform_phase_gap.py` freezes the repeated embedded Rx(pi/7),
b=3, NO-route family. At each of 2048 circular phase points it encloses both
alternating twelve-channel products in the existing eight-dimensional
orthonormal traceless Hermitian basis. P192 and P256 are separate checks;
the circuit builder uses period 3N to represent each exact phase j/N.

The maximum mesh Frobenius upper bound is approximately 0.7659369241.
C66's proved induced-norm perturbation adds at most 12*pi/N, bounded using
355/113, approximately 0.01840777102. The resulting uniform upper bound is
approximately 0.7843446951, below the preregistered 4/5 threshold. This is a
certificate for ALL sector phases and supplied M in the fixed family, not
an inference that a dense mesh probably missed no exceptional sector.

The no-background control retains an explicit conserved traceless mode at
a NONZERO sector phase. The undersized N=1 mesh cannot certify contraction
after its between-mesh padding is charged. Branch unitarity, basis geometry
and imaginary/zero-containment checks all pass. The experiment does not
differentiate the raw exponentially powered branch phase: C66 first removes
only the scalar phases that cancel inside forward conjugation.

Authoritative report: `out/uniform_phase_gap_20260911T055228005798Z.json`,
3/3 checks; log `out/uniform_phase_gap_audited.log`. The reported run takes
about 8.5 seconds and includes both precision meshes. A 32 MiB structural
preflight covers their retained records plus serialization; it is not RSS.
The actual conservative warm-up consequence remains the one in §FM, now
with the numerical constant justified uniformly. It does not produce a
useful memory saving within the current API. No width limit was increased.

## Exact checkpoint/recompute implementation

`VerifiedFiniteWork(..., checkpoint_spacing=k)` is opt-in; None preserves
the default. Shared exact branch and forward-update helpers keep the existing
target and finite-TV certificate. The cursor retains block starts and one
reconstructed descending block, rebuilding branches on demand. Precision
refinement preserves selected bits and replays the backward effect. Lazy
helpers use weak references so discarded cursors do not retain native
matrices in reference cycles. C67 gives the proof and complete cost scope.

The comparison's tiny fixture is r9,b3,t4 with Rx(pi/7) after EVERY control
and two reflection insertions. This is a deliberately specified fixture,
not an assertion that it is identical to §OB. It covers all four histories,
three initial sectors and sixteen outputs at target TV 1e-6. Full-storage
and k=1,2,4 laws match as exact Fractions; all are enclosed within the existing
prefix target's accuracy budget. Every target component mass encloses one.
This file uses that distinct prefix contraction, not a second full-r circuit
simulation; the separate b=2/3 full-law regressions cover that comparison.

The reset control changes NORMALIZED conditional weights at a fixed prefix,
not merely their scale. Correct weights are (228058542726,325014642565);
incorrect reset weights are (575577455539,523934172237). Widening E at fixed
p=40, P=64 after one bit forces a genuine rebuild to P128 and one effect
replay; it agrees with a fresh same-prefix high-precision cursor. The test
does not change the declared bit-accuracy contract after sampling begins.

No-background t4 complete laws normalize and agree, but do not exercise
an all-zero rounded block. A separately bounded t32 forced rare path does:
all compared spacings record five declared fallbacks and the same exact
forced probability. That is not an exhaustive t32 law or a sampled frequency.

The supplied wide case uses r=3*(2^60-1), b3, t63, the repeated mixer and no
routes. At target TV 1e-6, output-zero forced probabilities match exactly and
three matched seeded samples return the same output and final sector. Full
storage versus k8 gives the following measured counters for the forced path:

| Counter | Full storage | k8 |
|---|---:|---:|
| Peak retained forward matrices | 64 | 15 |
| Persistent branch matrices | 126 | 0 |
| Conservative structural matrix bound | 205 | 40 |
| Total forward updates across refinements | 126 | 187 |
| Recomputed forward updates | 0 | 61 |
| Branch-pair constructions | 126 | 251 |

Both paths rebuild once and end at P202. The recomputation count includes
the abandoned lower-precision block as well as the final pass. A single
worker-construction-plus-forced-path timing is about 0.00434 seconds full
versus 0.00692 seconds checkpointed. This is not a median, end-to-end order
discovery benchmark, speedup claim or native-memory measurement. Input gate
descriptions are supplied once outside that timer and remain storage costs.
The full output law is not enumerated at width63.

Authoritative report: `out/checkpoint_comparison_20260911T055820817140Z.json`,
4/4 checks; log `out/checkpoint_comparison_audited2.log`.

## Audit trail and preserved failures

The initial mesh attempt failed with an Acb API call error:
`out/uniform_phase_gap_failure_20260911T054933655247Z.json`. The agent's first
passing report is `out/uniform_phase_gap_20260911T054948489378Z.json`.
Main repaired the written Lipschitz constant/norm, a null eigenphase
denominator hidden by sector zero, the missing serialization allowance and
the undersized-grid control's actual charged bound. All original thresholds
and mesh/product lengths were retained.

The checkpoint invocation before its module existed is preserved in
`out/checkpoint_comparison_test_20260911T055213Z.log`. The first scientific
run, `out/checkpoint_comparison_20260911T055239729562Z.json`, failed its
zero-fallback prediction. The explicit rare-prefix probe repaired coverage;
the prior agent pass is `out/checkpoint_comparison_20260911T055442675026Z.json`.
Main then strengthened the savings/recomputation predicates, fixed replay to
hold the accuracy contract constant, charged construction in the timing and
bounded retained full laws before allocation. A main syntax error is retained
in `out/checkpoint_comparison_audited.log`; the successful rerun is separate.

An agent summary mistakenly called the branch-matrix count the forward-state
count. The authoritative fields and table above distinguish them. A passing
report never licenses silently treating structural matrix counts as RSS.

## Validation and scope

Core passed before science (`out/uniform_gap_core.log`). Backend-present lab
and claims pass (`out/uniform_gap_lab.log`, `out/uniform_gap_claims.log`). Lab
also passes with the optional backend absent and verified arithmetic explicitly
skipped (`out/uniform_gap_lab_no_flint.log`). Lab covers checkpoint bounds,
boundary routes, widths zero through five, precision replay and cursor lifetime.
Existing binary and odd-block full-law regressions pass:
`out/uniform_gap_binary_regression.log`, `out/uniform_gap_odd_regression.log`.
The other six science suites were not rerun. The documentation-gate log is
`out/uniform_gap_docs.log`; no manuscripts, abstract workshops or commits
were changed in this follow-up.

Reproduce the new probes with this environment, changing only the module name:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_uniform_phase_gap
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_checkpoint_comparison
```

The mathematical gap is useful structural progress, while checkpointing is a
known scheduling technique with demonstrated local savings. Neither establishes
a general simulation breakthrough. A stronger reverse-instrument lead emerged
from this baseline comparison; C68/§RI record its initial audit and TODO 28
alone specifies its next bounded certification question. TODO 24 keeps the
broader backend bit-cost issue deferred.
