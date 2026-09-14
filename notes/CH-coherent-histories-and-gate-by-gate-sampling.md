---
code: CH
date: 2026-09-11
title: "Coherent routes: a prefix-amplitude connection removes rejection, but not all its costs"
outcome: mixed
claims: [C59]
todo: [14, 23, 24]
---
# CH — Static sector closure is not finite-time sampling cost

TODO 23 began with coherent-component rejection. Main's primary-literature
check found a stronger alternative: Bravyi–Gosset–Liu gate-by-gate sampling
uses prefix amplitudes rather than output marginals. The route-history
contraction can provide those amplitudes without orbit or output tables.
C59 owns the proof, input promises, costs and citations. The result is a
useful specialized implementation, not a new sampling principle or a general
quantum-simulation breakthrough.

Three existing lower-cost `gpt-5.6-luna` agents independently tested formulas,
static regrouping, actual circuits and the induced sampling law. Main
implemented the production helper, read and corrected their verifiers,
reran all four final experiments, and added lab/claims regressions. The
qsim-research skill's prior-art and strongest-baseline requirements changed
the plan materially: both regrouping and the known gate-by-gate algorithm
were evaluated before interpreting the proposed rejection scheme.

## Finite formula and regrouping evidence

The formula fixture has r=10,b=2,t=5 and repeated Rx(pi/2) at insertions
1,3,4. Reflections q=0@s2 and q=1@s3 share one swept strength theta in
{-pi/4,0,pi/8,pi/4,pi/2}. A separate k=0..3 test takes prefixes of the fixed
schedule ((2,0),(3,1),(4,0)) at theta=pi/4. These are separate controlled
series, not a claimed scaling law. Initial/end insertions are also checked.

The verifier explicitly enumerates the tiny initial sectors and exponent
strings, using full-r matrices and existing `sequential_path`; it is NOT the
compressed production algorithm. Production target/proposal errors are below
6.9e-16, and its unnormalized accepted-law identity is within 7e-18.
The independent full-r comparisons and normalization checks also pass.
Deleting cross terms changes normalized output TV by up to 0.198133.
The wrong-inverse-route control has raw mass 0.2, so its raw half-L1
discrepancy 0.452169 is NOT TV; after explicitly normalizing, TV reaches
0.810436. Uniform final gamma differs from its target marginal by TV 0.148445.

The regrouping verifier checks one- and two-reflection groups for r=10,14
against full-r contraction. Maximum two-reflection discrepancy is below
5.9e-16. Dropping coherence inside a group changes TV by 0.099915; the other
negative control is a uniform OUTPUT law, not a uniform final-sector law.
The gcd orbit-size formula agrees with finite BFS in 468 cases, including
composite M, and is evaluated on supplied huge M using integers only.
One reflection label closes groups of at most two sectors; q=0,1 joins all
sectors. This only defeats that static-grouping baseline, not every possible
classical representation.

## Independent prefixes, sampling transitions and actual gates

For r=10,14,b=2,t=4, the sampling audit applies full-r coherent gates to
|orbit 0> and then projects the sector. It does NOT reuse the production
history expansion for this reference. Prefixes cover every arithmetic/work
boundary, all partial Fourier-output prefixes, all sectors and selected
fixed exponent strings (0,1,2^t-1). Maximum prefix error is below 1.7e-16.
Complete joint laws enumerate every exponent string and output.

A separate deterministic enumeration of every classical sampling transition
uses the production prefix oracle. Its final joint law matches the independent
reference below 1e-16, and its output marginal matches `sequential_path`.
This verifies the update rule's full tiny law, not just a sample histogram.
k=0, W0, K0, end-only and self-loop fixtures receive the same prefix/joint/
transition/mass checks. All reference/transition masses are within 5e-16 of
one. Normalized negative controls change TV by 0.757395 (frozen sector) and
0.0968264 (coarse dephasing or deleting history interference in this fixture).

The physical-circuit audit uses the indexed r=8,b=2 orbit with work bits
(p,m_LSB,m_MSB), a clean increment ancilla and four exponent bits. Existing
X/CNOT/Toffoli implement ascending controlled addition mod 8; every clean
work/control basis input is checked for powers 1,2,4,8,16. Basis-probability
and normalization errors are below 4.3e-15. This does not test arbitrary
dirty-ancilla semantics.

Two distinct coherent reflections are compiled using existing Pauli rotations:
q=0@s1 and q=1@s2, with alternating fixed Rx/Rz p mixers. Although J_1's
individual Pauli terms need not all commute, they split into two orthogonal
projector groups. Both the product-zero/group-commutator identities and the
complete compiled unitary are checked. Holding theta0=pi/5 fixed, vary theta1
over {-pi/3,-pi/7,0,pi/8,pi/3}. Circuit/reference output error is below
8.1e-16; production/reference is below 1.7e-16. Omitting or moving the second
reflection changes a probability by 0.0541266. The r=8 construction is not
a constant-size physical realization for general M.

## Supplied-wide draws and the cost tradeoff

Both samplers draw two outputs for each k=2,3, using supplied
r=2,000,000,014,b=2,t=63. All reflection angles are pi/4; the background has
three fixed repeated blocks. No orbit, sector or output table is constructed.
The gate-by-gate draws use respectively 133 and 135 prefix-vector queries,
each with up to four/eight histories. Rejection uses 4,4 and 5,5 proposals
in the fixed seeded diagnostic. These counts are not estimates of its mean.

Main's two-draw totals were approximately 0.242/0.476 seconds gate-by-gate
versus 0.031/0.050 seconds rejection. Setup is excluded and the sample is
deliberately tiny: this is a smoke test and a counterexample to an assumed
universal practical speedup, not a performance benchmark or crossover study.
The 256-byte reported owned matrix payload is ONLY stored background and
identity arrays; temporary matrices, instrument effects, Python objects and
process memory are additional. C59 gives the broader scalar storage bound.

Finite dense references are capped before allocation: r<=14,t<=4 for the
sampling-law audit, r<=14 for regrouping, fixed r=10,t=5 for the formula
fixture, and at most nine qubits for the indexed physical circuit (t<=5).
Individual array guards are 16 MiB, not a process-memory guarantee. No
amplitude cutoff is used. Integer phase reduction avoids unreduced huge
angles but does not certify rounding, cancellation or sampling accuracy.

## Failures and audit corrections retained

- Early formula reports `out/coherent_routes_formula_20260911T024614583753Z.json`
  and `...024640189545Z.json` fail routing and/or Fourier normalization checks.
  Their corrected successors retain the same derived target identity.
- `...025551835418Z.json` fails an overstrong normalization requirement on
  the intentionally wrong-sector control and a comparison of acceptance
  RATIOS at numerical null events. Final reports retain those raw diagnostics
  but gate normalized control TV and accepted probability mass instead.
  Tiny absolute law error does not make a roundoff-null ratio meaningful.
- The first physical fixture used q=0,2 at later insertions and its visibility
  controls were vacuous: `out/coherent_route_circuit_20260911T024715274549Z.json`.
  The final fixture explicitly changes to q=0,1 and earlier insertions.
  Main derived the orthogonal-projector compilation after an initial concern
  about noncommuting Pauli terms. This was a construction issue, not a no-go.
- The first sampling report `out/coherent_route_sampling_20260911T025539026613Z.json`
  fails its independent joint-law reference by about 0.0934. Later early PASS
  reports, including `...025712217303Z.json` and `...025755381441Z.json`, are
  ALSO superseded: the verifier discarded weights below 1e-15 and called an
  unnormalized |c_h|-weighted control TV. Final code uses exact-positive
  branches, |c_h|^2 for the dephased-history control, and checks all masses.
- Main further replaced a history-based prefix reference by direct coherent
  full-r gates, fixed a supposed k=0 fixture that still had a rotation,
  required full edge transition-law checks, and added preallocation caps and
  harness predictions. Reports preceding these audits are historical, even
  where labeled PASS.
- Regrouping's early large-M checks were too weak and its character matrix
  was guarded as a vector. Main fixed the actual allocation guards and added
  finite composite-M gcd/BFS checks. The uniform-output control was relabeled
  correctly. The historical 15-check report is not the final audit.
- Main's initial lab regression failed because existing Circuit H/CNOT
  conventions carry a global phase (`out/coherent_routes_test_lab.log`).
  The corrected test aligns ONE phase over the entire prefix state, never
  separately by history or sector. Main also corrected a physical norm check
  that double-counted ancilla leakage already included in the marginal.

## Reproduction and handoff

Run from research/ with bounded Python/NumPy:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_coherent_routes
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_route_regrouping
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_coherent_route_circuit
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_coherent_route_sampling
```

Main-reviewed and rerun reports:

| experiment | checks | report |
|---|---:|---|
| history formula/rejection | 69/69 | `out/coherent_routes_formula_20260911T025959598748Z.json` |
| exact regrouping | 17/17 | `out/route_regrouping_20260911T025952835224Z.json` |
| compiled indexed circuit | 21/21 | `out/coherent_route_circuit_20260911T030447193969Z.json` |
| prefixes and sampling law | 6/6 aggregate | `out/coherent_route_sampling_20260911T030634804485Z.json` |

Logs: `out/coherent_routes_formula_main.log`,
`out/route_regrouping_main_audited.log`,
`out/coherent_route_circuit_main_audited.log`,
`out/coherent_route_sampling_main_audited.log`.
Core, lab and claims pass in `out/coherent_routes_core.log`,
`out/coherent_routes_test_lab_audited.log` and
`out/coherent_routes_test_claims.log`. The other six science suites were not
rerun for these helper additions. The documentation gate passes all ten checks
in `out/coherent_routes_docs.log`; indexes were regenerated and
`git diff --check` is clean. No manuscript, abstract workshop, commit or
publication changes. TODO 23 is closed at bounded scope; TODO 24 owns numerical
certification. The broader user goal remains active.
