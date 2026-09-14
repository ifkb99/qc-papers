---
code: SR
date: 2026-09-10
title: "Physical phases expose a broader shortcut: route a sector rather than conserve it"
outcome: mixed
claims: [C58]
todo: [14, 22, 23]
---
# SR — The pi revival was predicted, not fitted to the sweep

TODO 22 replaced the localized kick by Rz on one actual work qubit, keeping
the previous background and insertion fixed. Main derived the exact orbit
sign sequence before the measurements and noticed that the pi rotation
permutes coarse sectors. C58 owns the proof and its scope: final orthogonality
of control-independent bijective routes suffices even without individual-sector
conservation. It does not make generic sector mixing cheap or prove hardness
when the shortcut fails.

Three existing lower-cost `gpt-5.6-luna` agents tested physical circuits,
integer sign/period structure, and abstract routed-sector formulas. Main
implemented `RoutedOrbitCircuit`, audited the references and controls, reran
all three experiments, and added lab/claims regressions. The qsim-research
skill's strongest-baseline and must-fail-control rules motivated the routing,
regrouping and exponent-transfer comparisons; they also prevented the false
interpretation of an invisible but incorrect phase operator as an identity.

## Physical comparisons and the useful negative controls

The N=7,a=3,t=5 background is W01@s1,W12@s3,W01@s4, all Rx(pi/2) repeated
on b=3 orbit blocks. Physical work-bit-1 Rz(theta) is inserted after s=2.
The angle sweep includes zero and both signs of pi/8, pi/4, pi/2 and pi.
Existing compiled Circuit/statevec and full-r `sequential_path` agree to
below 2.3e-14 in every output probability.

Wrong initial coarse dephasing evolves the SAME branches from the density
(|orbit 0><orbit 0|+|orbit 3><orbit 3|)/2. Its TV discrepancy at positive
pi/8, pi/4 and pi/2 is about 0.0369, 0.0681 and 0.0963, respectively. It
vanishes within reference roundoff at zero and plus/minus pi. End-only Rz
changes no output law. A genuinely repeated diagonal phase agrees with both
the periodic helper and the coarse-dephased full-r reference at all angles.

The naive C53-style transfer substitutes Rz(theta) on early exponent bit 1
while leaving the earlier mixers unchanged. It fails: output TV reaches
about 0.327 and the full-state L2 difference reaches 1. This is one tested
replacement, not a no-go theorem for all exponent-side descriptions.

The next modulus is N=13,a=2,r=12,t=5, with the same orbit-block schedule and
physical bit-1 Rz(pi/4). Existing Circuit rotations and CNOT conjugations
implement the two-level background gates on the actual four work qubits.
Every one of the 16 work inputs is checked against an explicit matrix;
background gates fix invalid labels, and Rz has zero label leakage. The
compiled output agrees with full-r to 1.23e-14, and regrouped b'=6 sectors
agree with full-r to 2.78e-17. C58 explains the exact period-six sign sequence.

An important control initially did NOT fail. Repeating only the first three
N=13 bit-1 signs changes the phase operator (max entry error about 0.765),
but not this early-insertion output (TV below 5e-17): the reached labels
still see the same signs. Retain this as an input-specific positive control.
The separately identified bit-0 phase does expose the wrong period-three
substitution, with TV about 0.0428. Changing the tested bit is an explicit
control choice, not part of the fixed-bit angle sweep.

## Routed helper audit and resources

The independent routing experiment projects D_q between explicit sector
bases for b=2,3 and M=2,3,4, including negative charges and multiples of M.
Maximum projection error is below 1.7e-15. Full-r original-order contraction,
an independently assembled routed contraction and production agree below
8.4e-17 on the multi-route rows. An r=4,b=2 existing Circuit/statevec fixture
gives a separate physical reference. A wrong routing target fails at matrix
level; freezing the sector despite routing fails at probability level.

The width-63 execution uses supplied r=3,000,000,021,b=3 and several integer
routes. It draws an actual sample with no orbit/output table. The finite-work
routine reports a conservative 47,664-byte matrix payload estimate, not peak
RSS or total process memory. This is not a precision guarantee or order
discovery. Integer phase reduction is essential even for this diagnostic:
direct trigonometry of an exponentially large unreduced angle is invalid.

Dense references stay tiny: r<=18 in the physical/routing experiment and
r<=24 in the sign audit. Individual dense allocations are capped at 16 MiB;
these are per-allocation checks, not a 16 MiB process-peak claim. The actual
N=7/t=5 state has 13 qubits and 131,072 complex-vector bytes; N=13/t=5 has
15 qubits and 524,288 bytes. N=13 compiles to 22,496 gates, below the revised
30,000-gate fixture cap. Widths 4–7 are preflighted, NOT a measured physical
width sweep. The sign audit enumerates the two tiny orbits once and charges
18 modular orbit steps; no generic efficient recognition algorithm is claimed.

## Failures and main-audit corrections retained

- The sign agent's earliest log (`out/phase_structure_test.log`) includes
  overstrong half-complement checks on N=13 higher bits. The exact sign
  sequences refute that extension; C58 limits the general parity argument.
  Two early failures (`...022614096276Z`, `...022644924875Z`) confused modulus
  with orbit length; `...022719087901Z` had a wrong grouped block shape.
  All are retained as `out/phase_structure_failure_20260911T*.json`.
- The vacuous bit-1 control is visible in
  `out/phase_structure_test_20260911T0225Z.log`. The associated
  `out/phase_structure_failure_20260911T022733777240Z.json` actually records a
  subsequent complex-JSON serialization exception, NOT the failed-control
  detail. An agent summary conflated those; the raw log is authoritative.
  Main also replaced tautological period checks and fixed a cap off-by-one.
- The routing experiment's first fractional-phase schedule canceled in the
  output probabilities. Its failed report is
  `out/sector_routing_20260911T022843809698Z.json`. The corrected, explicitly
  changed insertion exposes the wrong deterministic substitution. A separate
  wrong-sign probability check remains vacuous; only its matrix check is
  counted. Both facts are retained rather than hidden behind the pass total.
- Main caught a character-matrix allocation guarded as a vector after an
  allocation, and unreduced large phases in the independent wide reference.
  Guards now precede the actual matrix allocation; integer reduction is used
  before trigonometry. Earlier wide-reference reports are historical.
- N=13 first exceeded the proposed gate budget before state propagation
  (`out/physical_phase_failure_20260911T023301695151Z.json`). The bounded gate
  cap was raised to 30,000. The next independent physical comparison caught
  incorrect local-to-physical work-qubit embedding, with output error 0.0414:
  `out/physical_phase_run_4.log` and
  `out/physical_phase_20260911T023325579579Z.json`. Corrected embedding passes.
- Main fixed an ignored theta argument in the new N=13 reference, added the
  missing explicit -pi assertion and checked the N=13 Rz matrix as well as
  leakage. An initial lab regression expected bool indices to be rejected,
  unlike the inherited integer API; the failed log
  `out/physical_phase_test_lab.log` is retained. The corrected validation test
  uses a noninteger float. No production integer semantics were changed.

## Reproduction and status

From research/, run each module with Python 3.12 and bounded NumPy:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_physical_phase
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_phase_structure
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_sector_routing
```

Main-reviewed reports:

| experiment | checks | report |
|---|---:|---|
| physical phase and N=13 circuit | 38/38 | `out/physical_phase_20260911T023618836612Z.json` |
| exact phase structure | 14/14 | `out/phase_structure_20260911T023609004174Z.json` |
| deterministic sector routing | 99/99 | `out/sector_routing_20260911T023350927347Z.json` |

Logs are `out/physical_phase_main.log`, `out/phase_structure_main_audited.log`
and `out/sector_routing_main.log`. Core, lab and claims pass in
`out/physical_phase_core.log`, `out/physical_phase_test_lab_audited.log` and
`out/physical_phase_test_claims.log`. The other six science suites were not
rerun for this helper addition. C58 records the primary-source positioning:
standard block-permutation mathematics, not a new normalizer/factoring theorem.
The documentation gate passes all ten checks in `out/physical_phase_docs.log`;
indexes were regenerated and `git diff --check` is clean.
No manuscript, abstract workshop, commit or publication change in this follow-up.
The broader goal stays active; TODO 23 owns the next discriminating question.
