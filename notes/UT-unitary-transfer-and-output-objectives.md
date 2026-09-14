---
code: UT
date: 2026-09-11
title: "The best joint-state transfer loses to omission on the measured output"
outcome: confirmed
claims: [C53, C76]
todo: [14, 34, 36]
---
# UT — Optimize the observable, not automatically its purification

C53 owns the general early-unitary transfer criterion, root-fidelity optimum
and parity obstruction; C76 owns the additive gate's dense sector boundary.
This is a stronger-baseline follow-up to GS, not a new phase-estimation or
simulation algorithm. Main derived the matrix optimization, requested a
lower-cost initial experiment and independent proof audit, then corrected
the verifier and added the actual physical-output comparison.

## Same circuit, stronger candidate

The input and insertion are exactly GS's N13,a2,b3,r12,t3 schedule with
s=2 early controls and k=0,1,2. The four-column pre-kick boundary A is formed
from the existing branch matrices, normalized to Frobenius norm one. A work
G_k gives B=G_k A. An arbitrary early-exponent E gives A E^T, allowing more
than the failed scalar phase-per-history replacement. This operation is
at the insertion, not commuted through the earlier arithmetic.

The SVD of A^dag B supplies the best joint-state candidate. The mathematical
identities and exact-transfer criterion follow from standard purification
freedom/polar optimization, not a new theorem. Main read Uhlmann's definition,
explicit fidelity formula and amplitude/polar discussion, Sections 2–3,
equations (3), (12), (17)–(20), in
[Transition Probability (Fidelity) and its Relatives](https://arxiv.org/pdf/1106.0979).
The unsquared root-fidelity convention is explicit. No claim relies on the
paper's later parallel-transport discussion. An independent lower-cost proof
audit agreed, including the transpose and rank-deficient cases.

## Main finite findings

The final report is
`out/additive_phase_transfer_20260911T095307109554Z.json` (5/5 PASS).
At k1 and k2, the joint-optimal root fidelities are .9441123032 and
.9170944664. The best diagonal-only values are .8091199332 and .7230247906;
omission gives .3135672902 and .5970156977. Nevertheless, neither general
candidate transfers exactly: squared vector residuals are .1117753935 and
.1658110672, with work-density commutator Frobenius norms .2969898042 and
.3763019961. These are small floating diagnostics, not an asymptotic theorem.

Main embedded each target/general/diagonal/omitted boundary into the same
ten-qubit physical layout, with the unprocessed high control in |+>, then
used existing Circuit/statevec for the remaining controlled power, terminal
mixer and inverse QFT. The target agrees with the independent full-r law
within 2.074e-14 per outcome; every candidate law normalizes and cleans its
scratch. The joint-optimal candidate's output TVs are .1511231980 (k1) and
.1581159429 (k2). The diagonal-only TVs are .2051832075 and .1086860902.
Omission is substantially BETTER on this measured task: its baseline values
are owned by GS and reproduced in this report. No ordering of output errors
was predicted from the state-fidelity ordering; the observed reversal is
retained rather than selecting a different input.

The k0 state and output are nulls. A separately normalized four-column
Fourier purification supported on four work labels has diagonal work density.
It permits exact nontrivial general transfer (root fidelity one) while
diagonal-only transfer has squared residual 1.2007600702. This synthetic
control demonstrates that the stronger test is not just the old scalar test;
it is not substituted for the physical initial state.

## Independent output formula and an observable-specific restriction

After the physical comparison, main asked a lower-cost agent to check the
surprise by a finite closed formula. With X any candidate boundary and
v_y[l]=exp(-2*pi*i*y*l/8), the output law is

    p_y = ||(I+(-1)^y U^4) X v_y / 4||^2.

The last common W cancels under the norm; it need not commute with U^4.
The agent's independent read-only calculation matched the complete physical
laws and ordering within 2.2e-14. That pilot has no standalone JSON; main
then added the formula and normalization/error predicates to the released
experiment. This is a fixed boundary identity, not a second propagator.

Summing even y suggested a sharper test, derived before that marginal was
measured: output parity depends only on the work density at the split. An
early-register-only trace-preserving operation cannot change it. Main and an
independent lower-cost audit derived C53's no-signalling lower bound and
compared the direct work-overlap formula to the physical parity marginal.
Their agreement is within 5.5e-14. The target/background gaps evaluate to
.001399550303399 (k1) and .001246125230684 (k2). Thus the formula obstructs
the entire restricted early-register-only class on these numerical witnesses,
not merely the particular polar candidates. These decimals are NOT interval-
certified lower bounds, and the restriction excludes work-coupled rewrites,
high-register operations, postselection and general classical samplers.

## Budgets, precision and corrected failures

The final twelve suffix propagations each use 3,849 Pauli rotations and
jointly charge 47,296,512 gate-entry updates under the 100,000,000 cap.
The aggregate pre-allocation numerical estimate is 1,459,712 bytes under
16 MiB. Instrumented diagnostic matrix products total 157 calls and 54,800
scalar product terms; parity inner products charge 576 entries. Eight SVD
dimensions and four eigensystem dimensions are recorded separately. The
independent reference makes 24 sequential_path calls and charges 124,416
t*d^3 units. These are not native FLOP/bit-runtime counts, total process
memory, physical synthesis costs for E or a timing comparison. Setup of the
supplied tiny reference matrices is not claimed free or included in the
diagnostic matrix-loop count.

An independent square-root-density route checks root fidelity; its tolerance
was conservatively set to 2e-7 because of rank-deficient square-root roundoff.
Actual errors are below 3.4e-16. Residual identities use 3e-12 tolerance;
physical laws and the parity/formula comparisons use 3e-10. Near the k0
null, cancellation in sqrt(1-F^2) can produce a spurious approximately
2e-8 bound or zero, so these are float diagnostics, not numerical certificates.

The lower-cost initial report
`out/additive_phase_transfer_20260911T094223554480Z.json` failed P1 because
its synthetic Fourier matrix had Frobenius norm squared four rather than
one. Its log `out/additive_phase_transfer_initial.log` is preserved. Main
caught that normalization mistake before interpreting the synthetic result.
The initial actual-state measurements are retained, but its resource and
identity predicates were incomplete: an accepted classification string was
tautological, residual identities were recorded but not checked, some guards
followed allocation, and work counts were an uninstrumented expression.
Main corrected all of these, handled a zero diagonal overlap with a unit phase,
removed import-time experiment creation, and added the physical suffix check.

The first corrected physical report
`out/additive_phase_transfer_20260911T094836405466Z.json` passes 4/4; the
final report additionally includes the independently derived closed formula
and parity test. The new predicates were declared before their new
measurements. The original k values, mixer schedule and insertion did not
change. A separate read-only agent audited the corrected physical embedding,
normalization, optimization and counters and found no substantive error.

## Validation and handoff

GS's core-first pass also covers this uninterrupted scientific follow-up.
The affected claims suite passes in `out/additive_transfer_claims.log`,
including independent two-state general-transfer and state-optimal-but-
output-worse regressions and a nonzero parity witness. No shared production helper changed; the other
seven science suites were not rerun. System Python 3.12.3 and NumPy 2.4.6
with one BLAS thread were used; the claims gate also uses python-flint 0.9.0.
The documentation gate is recorded in `out/additive_phase_docs.log` after
regenerating indexes. No manuscript, abstract, default, host-setting or
firmware changes or commits were made. TODO 34 remains unresolved.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_additive_phase_transfer
```

The skill's same-output/strong-baseline requirement materially changed this
investigation: a state-only optimization initially looked favorable, but the
complete measured law favored omission. C53 now records the stronger criterion
and restricted output obstruction. TODO 36 retains observable-native
contraction as the next question, not further tuning of this state optimum.
No new efficient sampler, general hardness claim or breakthrough is established.
