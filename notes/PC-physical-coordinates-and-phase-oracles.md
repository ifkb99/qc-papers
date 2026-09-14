---
code: PC
date: 2026-09-11
title: "Cell reflection is easy to evaluate; a phase evaluator can hide discrete logs"
outcome: confirmed
claims: [C52, C58, C59, C73, C74]
todo: [14, 32, 33, 34]
---
# PC — Audit what the gate description supplies

C74 owns the elementary proofs and their oracle/physical distinction. While
auditing the wider word-envelope comparison, main asked whether its supplied
indexed gates could be evaluated directly on actual modular orbit encodings.
The useful split is between fine-coordinate/cell-reflection evaluation and
high-order phase evaluation. This is a scope clarification, not a new
discrete-log algorithm or a lower bound on quantum gate synthesis.

A lower-cost agent independently checked the algebra and prepared an exact
integer/Fraction probe. Main corrected the proposed inverse-doubling error
threshold, strengthened the verifier and reran it. The decoder receives only
noisy oracle values, not the underlying coordinate. Those values are generated
from known test coordinates: the experiment DOES NOT implement the hypothetical
physical x-to-phase evaluator required by the reduction.

## Exact finite evidence

Main's `out/physical_orbit_coordinates_20260911T084518735050Z.json` passes
7/7 checks using system Python 3.12.3, with NumPy 2.4.6 imported by the harness.
The scientific calculations are integers/Fractions, not floating contractions.
Seven frozen prime/composite fixtures test every orbit label, with orders up
to 96, fine dimension up to six and coarse modulus up to 32. Exact-order and
simultaneous orbit/table preflights run before the relevant allocations.

The fine-coordinate table and physical cell reflection match every reference
label. Plain inversion fails on each fixture. All tested nontrivial q phases
recover their promised coordinate residues at rational error 3/(8s), including
wraparound; an error of one full grid step supplies the insufficient-precision
control. Subgroup-table recovery and independently generated reflection-group
closure agree with C74/C59. The extra N=97,a=5,b=3,M=32 fixture uses labels
(0,8,16): their phase orders already permit a small static grouping. This
guards against presenting an easy low-order family as a new simulation win.

Inverse doubling checks 4,636 residue/error-pattern cases: all signs of
1/16 errors for M=2,...,16, and four fixed patterns at M=17,31,64,127.
It checks both grid recovery AND the contracted circular-error bound.
The epsilon=.24 control has a concrete failure already at M=2,m=0:
observations (19/25,6/25) produce estimate 31/50 and the wrong grid residue.
This refutes the agent's initial epsilon<1/4 sufficient-condition suggestion;
the claim uses the proved epsilon<1/6 condition and tests epsilon=1/16.

Charged fixture/reference updates total 16,544; decoder/reference updates
total 85,526, below separate caps 250,000 and 150,000. They are explicitly
not native instruction, bit-runtime or RSS measurements. The small full-orbit
reference is test-only; it is not hidden inside the proposed b-entry evaluator.

## Corrections and limits

The first initial report, `out/physical_orbit_coordinates_20260911T084134955450Z.json`,
stopped at the unwrap preflight: its 206,348-operation plan exceeded 150,000.
Earlier fixture computations had run; no unwrap measurements were completed.
The requested four-pattern wider coverage had not been fully wired into the
caller and still expanded all 64 patterns for two wider moduli. Main corrected
that mismatch without raising the cap. The first main pass is preserved as
`out/physical_orbit_coordinates_20260911T084229232279Z.json`; the final pass
adds the larger small-phase-order fixture, exact error-bound predicates and
raw negative-control observations.

Before that first run, audits caught missing Fraction serialization, a static
comparison made by rebuilding the claimed cosets rather than generating the
group, incomplete allocation/cost guards and excess repeated setup powers.
These were verifier defects, not evidence against the reductions.

This test checks the coordinate/phase algebra and oracle decoder. It
does not independently compile the full physical fine-work gate or a clean-
ancilla reflection rotation, and does not establish that arbitrary supplied
high-order characters admit a cheap classical evaluator. TODO 33 retains
those distinctions; C75/CG subsequently complete the bounded compiled-circuit
question without changing this algebra/decoder experiment.

Core/law checks in the unchanged system-Python environment pass in
`out/word_envelope_core_system_python.log` and
`out/word_envelope_laws_system_python.log`. The new exact claim regression
is recorded in `out/physical_coordinates_claims.log`. WE records earlier
affected lab/backend-absent and complete-law gates. No additional six-suite
coverage is implied. The host/runtime reliability question remains open in
TODO 34; exact references and redundant checks do not prove hardware health.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_physical_orbit_coordinates
```

Main inspected the normalizer paper's coordinate-conversion/phase premises
and the classical discrete-log bit-security theorem statements. C74 records
the exact prior-art reading scope. No new cryptographic-hardness or simulation
breakthrough is claimed. No manuscript, abstract or production sampler was
changed by this follow-up.
