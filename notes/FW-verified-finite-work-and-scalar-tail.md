---
code: FW
date: 2026-09-11
title: "Verified finite-work proposals work, but a late-work invariance makes the binary benchmark simpler still"
outcome: mixed
claims: [C56, C63, C64, C65]
todo: [14, 24, 25]
---
# FW — Certification survives; the benchmark needs a stronger baseline

C64 owns the unnormalized-mass error argument and verified forward/backward
implementation. C65 owns the collective late-work invariance and stronger
binary scalar proposal. Neither result establishes general simulation hardness,
a new sampling principle, or a factoring breakthrough.

The qsim-research skill materially changed this investigation: checking the
strongest simple baseline exposed that the fixed binary benchmark's late
background rotations are invisible to its requested output. Their presence
still exercises the numerical implementation, but cannot demonstrate persistent
observable fine-work mixing. We retained that fixture and added an observable
initial rotation instead of treating a successful certification as the end of
the investigation. The established finite-work identity is reused; no second
generic propagator was written.

## Implementation and evidence

Three lower-cost agents supplied an exact-rational tree test, complete-law
verifiers and a read-only proof/code audit. Main implemented and audited the
verified density/effect and scalar helpers, corrected verifier weaknesses,
and reran the affected checks. All code remains uncommitted. Both new proposal
methods are opt-in in `VerifiedRejectionSampler`; the default prefix component
and the uncertified floating sampler retain their contracts.

The mass-tree experiment fixes depth 5 and finite-bit resolution 12, varying
only requested weight bits through 0..12. Exact Fractions cover tiny masses,
genuine zero internal prefixes and uniform fallbacks. Its rare deterministic
node has local conditional TV 1/2 but tiny true-weighted error; a separate
untracked-rescaling control violates the claimed absolute-mass budget.

The component experiment enumerates every conditional and joint law for
r=10, b=2, t=4, all four route histories and five initial sectors, at requested
TV 1e-3 and 1e-6. It tests both the original fixture and the same fixture with
W0=Rx(pi/4), using both verified implementations. Every actual finite-bit law
normalizes exactly. Outward ball bounds against the ideal component law obey
the respective plans. Independent full-r complex128 calculations agree with
ideal midpoints within 5.6e-17; those are independent numerical diagnostics,
not interval certificates. Forced paths also exercise exact-zero fallbacks;
their count is not a claim that those prefixes occur in actual ideal samples.

The full-r late-background experiment finds a joint-output difference at
roundoff when removing W1 and W3 together. In contrast, adding W0 changes
that joint law by TV about 0.2067. Keeping the fine work-basis label reveals a
maximum joint-cell difference about 0.02151 under removal of the late gates.
These controls distinguish the specified traced output from the full state.

End-to-end rejection tests enumerate the complete accepted laws for both
fixtures and requested accuracies with both new component methods. All eight
rows normalize exactly and satisfy the proposal, success and final-law plans.
At target 1e-6, their outward accepted-law TV upper bounds are between
1.777e-8 and 1.959e-8. The acceptance computation still includes all coherent
histories; using scalar component proposals does not replace interference by
an incoherent mixture.

A valid deliberately loosened backward effect forces refinement after a bit
has already been selected. Rebuilding at higher precision and replaying that
bit reproduces the fresh high-precision weights exactly. Resetting the effect
instead gives wrong UNNORMALIZED weights. In this binary scalar-tail family,
resetting may leave conditional ratios unchanged; the control establishes
violation of the mass-oracle promise, not necessarily a changed conditional
ratio on this particular fixture.

Authoritative main-reviewed reports and logs:

| Experiment | Checks | Report | Log |
|---|---:|---|---|
| `experiment_prefix_mass_budget` | 3/3 | `out/prefix_mass_budget_20260911T044405797995Z.json` | `out/prefix_mass_budget_audited.log` |
| `experiment_verified_finite_work` | 3/3 | `out/verified_finite_work_20260911T045224738535Z.json` | `out/verified_finite_work_audited.log` |
| `experiment_late_backgrounds` | 3/3 | `out/late_backgrounds_20260911T045224840527Z.json` | `out/late_backgrounds_audited.log` |
| `experiment_finite_work_comparison` | 12/12 | `out/finite_work_comparison_20260911T045058422401Z.json` | `out/finite_work_comparison.log` |
| `experiment_odd_block_tail` | 4/4 | `out/odd_block_tail_20260911T050316685319Z.json` | `out/odd_block_tail_audited.log` |

These reports store complete tiny-law evidence, exact bounds, resource
categories, precision/refinement counters and controls, not just histograms.

## Initial odd-block boundary test

After the binary comparison, a lower-cost agent tested a fixed M=3,t=4
coherent fixture with W1=embedded Rx(pi/7), W3=embedded Rz(pi/5), and
reflections K2(q=0,pi/5), K3(q=1,pi/5). Within each fixture only the joint
removal of W1/W3 changes. The paired b=2,r=6 and b=3,r=9 cases also change
the period when changing block size: this is a boundary comparison, not a
one-parameter scaling law.

For b=2 the complete joint law stays unchanged to roundoff. For b=3 the
removal changes it by TV about 0.05095, with maximum cell difference about
0.00801. Existing route contractions agree with independent full-r products
within 1.39e-16; all laws normalize. A final-only work rotation is invisible
in both cases, as required by the trace. These are complex128 diagnostics,
not a verified odd-block sampler. C56 already has visible odd-block examples;
this probe supplies a coherent-route fixture beyond C65's scalar tail, not
an independent claim of novelty or universal odd-block visibility.

Main read and reran the experiment, removed a vacuous nonnegative-error clause
from the binary predicate, separated the odd-block hypothesis into its actual
counterexample check, corrected the distinction between production tables and
allocated diagnostic full-r matrices, and used the standard harness report.
The original agent report
`out/odd_block_tail_20260911T050055095625Z.json` is retained. TODO 25 owns the
remaining exact-input extension and same-accuracy comparison.

## Matched timing and its limits

The comparison fixes supplied r=2^61-2, b=2, t=63 and target TV 1e-6. It
rotates the ordering of three certified rejection implementations: the C63
PREFIX-COMPONENT baseline, the verified finite-work proposal, and the scalar
proposal. Median whole-construction-and-sampling times are respectively about
0.1440, 0.06286 and 0.05117 seconds. The scalar implementation is about 2.81x
faster than that particular prefix-component rejection baseline, not than all
classical simulators or the whole-coherent C61 sampler.

Only three seeded random traces are used, each repeated twice; there are not
six independent rejection traces. The baseline's attempt counts are 3,2,4;
the two new methods' counts are 13,8,2. Different kernels consume different
random words, so equal seeds do not imply equal sampled attempts. An observed
attempt count exceeding an expected-count bound does not falsify that bound.
All rejected attempts and refinements are charged. An earlier single finite-
work smoke timing was already known when this timing hypothesis was declared;
the experiment header discloses it, and no timing order was predicted.

This is a bounded implementation comparison on the now-simplified fixture.
Counts of retained matrices/scalars are not peak process memory, nor are
linear-depth arithmetic counts backend bit-runtime theorems. No orbit,
output-probability or random-word table is allocated by these samplers. The
optional backend's uniform primitive-error constant remains unproved.

## Preserved verifier corrections

1. The first tree control reported large relative error in one tiny mass,
   which did not itself establish large CONDITIONAL-law error. Main replaced
   it with the rare deterministic node described above, and added an actual
   zero internal prefix. The earlier report
   `out/prefix_mass_budget_20260911T044138447156Z.json` is retained.
2. The tree preflight initially charged an exponentially enumerated random-
   word space as though it were a retained array. Main replaced that fictitious
   allocation with a conservative allowance for the actually retained entries;
   it is still not a measured RSS bound.
3. The initial late-background report labeled squared scaled coordinates as
   actual joint fine-work probabilities. `direct_prefix` returns a sqrt(M)-
   scaled vector, so the probability difference requires division by M. The
   corrected report labels both quantities. The earlier
   `out/late_backgrounds_20260911T044808775803Z.json` remains available.
4. A lower-only normalization check was replaced with an absolute deviation
   check for both distributions. Main also converted final successful reports
   to the standard harness schema. Earlier component reports
   `out/verified_finite_work_20260911T044538130343Z.json` and
   `out/verified_finite_work_20260911T045010973947Z.json`, and the corrected
   agent late-background report
   `out/late_backgrounds_20260911T045033866463Z.json`, are not overwritten.
5. The read-only auditor initially suspected extracting `arb.rad()` could
   round an uncertain radius inward. That was NOT a certification bug:
   the backend documents `rad()` as exact, and `binary_fraction` rejects an
   uncertain Arb value. The auditor corrected the claim after checking the
   contract and examples. Explicit `.upper()` is retained for consistency,
   not described as repair of a real containment failure.

## Validation and reproduction

Fresh core passed before scientific changes: `out/finite_work_core.log`.
Affected suites pass in `out/finite_work_test_lab.log` and
`out/finite_work_test_claims.log`. The lab suite also passes without the
optional backend, explicitly skipping verified arithmetic rather than claiming
it ran: `out/finite_work_no_flint_lab.log`. The other six science suites were
not rerun for these helper changes.

Existing full-law regressions also pass:
`out/finite_work_verified_sampling.log`,
`out/finite_work_coherent_sampling.log`, and
`out/finite_work_rejection_regression.log`. Shared branch construction and
exact-input trigonometric helpers preserve their prior mathematics. Run the
new experiments from research/ with, for example:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_finite_work_comparison
```

The other experiment names in the table use the same invocation. The backend
is an optional environment dependency, not a new mandatory project dependency.
The documentation gate is recorded in `out/finite_work_docs.log`. Neither
manuscript nor either abstract workshop was edited in this follow-up.
TODO 24 owns the deferred backend-bit-complexity boundary; TODO 25 owns the
next discriminating research question. Larger binary size records would not
resolve the benchmark limitation identified here.
