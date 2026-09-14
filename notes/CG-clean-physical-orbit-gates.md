---
code: CG
date: 2026-09-11
title: "Clean physical orbit gates close the coordinate gap, without a sampling speedup"
outcome: confirmed
claims: [C56, C59, C74, C75]
todo: [14, 33, 34, 35]
---
# CG — Test coherent cleanup, not only basis labels

C75 owns the clean-conjugation construction, full-domain reversible extension,
arithmetic existence bound and static-baseline caveat. Main derived that
construction; lower-cost agents audited it and supplied initial gate tests
and an independent output reference. This investigation used the existing
Circuit/statevec engine and finite-work sequential_path, not a new propagator.

The physical work register holds a modular residue, not its discrete-log
index. The finite compiler deliberately enumerates tiny truth tables; it
does not implement the polynomial reversible arithmetic whose existence
follows from C75. The user-supplied exact order, orbit promise and small W
implementation remain part of the contract. Neither papers nor abstract
workshops nor production samplers were changed.

## Main gate validation

`out/clean_orbit_gates_20260911T090928946014Z.json` passes 29/29 checks.
The frozen (N,a,b,r) fixtures are (7,3,3,6), (13,2,3,12), (15,2,2,4), with
angles 0, pi/7 and -pi/3. Each compute/strip permutation is checked on its
ENTIRE small binary domain using complex amplitudes. Every promised orbit
column of the clean mixer and reflection rotation is compared to a separately
assembled indexed matrix, including zero scratch rows and relative phases.
Checking the full columns also tests their action on arbitrary coherent
superpositions, not merely basis transition probabilities.

Maximum positive complex-amplitude error is 1.069e-13; maximum scratch
leakage probability is 3.55e-30. The recorded cumulative propagation count is
46,758,400 Pauli-gate-times-complex-entry updates under a 100,000,000 cap.
Maximum simultaneous numerical payload estimate is 794,624 bytes under
16 MiB; largest tested individual circuit has 3,933 Pauli rotations under
30,000. These are not FLOP counts, bit runtime, total Python-object memory
or RSS. Gate-list construction, table scans and modular-power setup counters
are reported separately; the compiler's exponential dependence is explicit.

Controls are operational failures: omitting coordinate erasure leaves scratch
probability at least .75 on a witness; plain inversion differs on orbit
labels; changing controlled V to controlled (-V) gives the wrong coherent
rotation. For the nontrivial N=13 reflection the latter has a .937
phase-invariant matrix discrepancy while agreeing with the independently
predicted opposite-sign rotation. M=2 makes R identity, so those fixtures
are NOT counted as meaningful wrong-sign controls: there the difference is
only a global phase. This restriction was written before the initial run.

## Full physical output and stronger baseline

`out/clean_orbit_output_20260911T090820589417Z.json` passes 3/3 checks.
It keeps N=13,a=2,b=3,r=12,t=3 fixed and varies only the second reflection
angle over 0, pi/7, pi/2. There is an initial fine mixer, a fine mixer after
each ascending arithmetic control, and reflections after the first and
second controls. The ten-qubit physical circuit ends with the exponent
inverse QFT. Its complete eight-output law agrees with both an independent
12-dimensional instrument contraction and the exact static-sector grouping.

Maximum compiled probability error is 1.086e-13; maximum grouped-reference
error is 1.111e-16. Groups are {0}, {1,3}, {2}, so the strongest baseline
uses coherent dimension at most six. Cross-group branch leakage is below
5.66e-16. Physical scratch probability is below 7.11e-30 and total-mass
error below 4.59e-13. An additional compiled omitted-reflection row matches
its own independent reference; at pi/7 omission changes the compiled law
by TV .0232048832787. This guards against validating only an invisible gate.

The three main physical circuits have 16,086 rotations; the omitted row
has 15,647. All four together charge 65,438,720 gate-entry updates. Each
physical run's numerical payload bound is 327,680 bytes; the reference
bound is separately 185,856 bytes. The reference actually makes 136
sequential_path calls, with sum(t*d^3)=351,216 as an explicitly named work
unit, NOT an exact operation count. Setup enumeration and finite projections
are diagnostic costs, not evidence of an orbit-free compiler or a timing win.
The grouping experiment enumerates all tiny groups to obtain the FULL law;
C59's large-instance baseline samples a group without that enumeration.

## Preserved failures and audit corrections

The initial unfactored primitive test passed in
`out/clean_orbit_gates_20260911T090247197366Z.json`, with log
`out/clean_orbit_gates_initial.log`. The factored lower-cost repeat passed in
`out/clean_orbit_gates_20260911T090734089131Z.json`. Its separate read-only
all-pair primitive diagnostic also checked the corrected phase convention;
that diagnostic is not a standalone harness artifact.

Two complete-output attempts failed the unchanged 30,000-gate construction
cap: `out/clean_orbit_output_20260911T090345546010Z.json` (lower-cost initial)
and `out/clean_orbit_output_20260911T090611719348Z.json` (main). Their first
finite reference calculations HAD run; no physical state propagation or
complete output verdict had occurred. They were budget failures, not native
crashes. Main retained both reports and their logs, then made two exact
compiler simplifications without changing inputs, tolerances or caps:

- Adjacent U^-1/U pairs between a mixer and reflection cancel, so that post-
  arithmetic block enters/exits coordinates once. No arithmetic is reordered.
- The two endpoint projectors giving a true transposition have exactly
  cancelling diagonal Pauli coefficients. Their sum is emitted directly.

The first simplification alone still exceeded the cap; the second fit it.
New primitive checks preceded the successful full-output test. The old
two-level helper's named-CNOT global phase was corrected locally, and
controlled Rx(pi) was explicitly converted to a true endpoint swap. Existing
shared Circuit/propagator behavior was not altered.

Before/following the first failed output run, main audits also caught an
unsupported Fixture keyword, a shape predicate comparing a tuple to an
integer (which would reject every probability law), missing internal
reference payload, and an unsupported scalar-work accounting statement.
These were verifier defects. Final reports retain raw probability vectors,
normalization/leakage checks and instrumented calls. A proof audit's initial
suggestion to combine all q=0 rotations was corrected: intervening arithmetic
prevents that simplification, despite commutation with the fine W blocks.

## Scoped gates and handoff

Core passed first in `out/clean_orbit_core.log`. The affected claim suite,
including the new independent C75 algebra/flag regression, passes in
`out/clean_orbit_claims.log`. The two final main experiments pass as above.
No production helper was changed; the other seven science suites were not
rerun. System Python 3.12.3, NumPy 2.4.6 and single-threaded BLAS were used;
the claim gate additionally used python-flint 0.9.0. These passing small runs
do not settle TODO 34's host/runtime problem. No firmware/settings changes or
new long timing benchmark were attempted.

The documentation gate passes all 10 checks in `out/clean_orbit_docs.log`;
indexes were regenerated with the repository tool. No commits were made.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_clean_orbit_gates
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_clean_orbit_output
```

The useful outcome is closure of the clean physical input-model gap, not a
breakthrough in sampling. TODO 35 selects a different cheaply implemented
physical phase, motivated by its additive/multiplicative Gauss-sum structure.
The next question and its pilot limits live there rather than extending this
completed experiment or rewriting the manuscripts.
