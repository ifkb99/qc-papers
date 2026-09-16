---
code: DB
date: 2026-09-14
title: "The differential bridge: DDT counts for off-diagonal Pauli pullbacks"
outcome: identity, counting theorem and APN extremality derived before measurement; all predictions and must-fail controls pass on dense, gate-level and a 14-qubit GPU dense reference that settled the out-of-sample F10 prediction
claims: [C99]
todo: [51]
---
# DB — Differential cryptanalysis counts X/Y-type Pauli paths

C99 owns the statement. This note records how it was found and checked.

## How it was found

The session started by cancelling two superseded Claude board tasks
(Tff977bfd7b8b42b5, T15db99d348534ca4). The user then asked for an original
idea. The ledger's cryptanalysis dictionary covers only the linear half:
C12 and C25 connect Walsh sparsity and nonlinearity to diagonal observables.
Its stated boundary is C13/F10, where X/Y pullbacks "blow past the cap".
Differential cryptanalysis is the other half of that dictionary, so the
question was whether the off-diagonal cost is its exact counterpart.

A permutation conjugates any Pauli into a signed permutation matrix, so the
pullback can be written in closed form before any measurement. The pair
structure of the involution σ gave the odd-k and k = 2 class counts. The
uncertainty principle and Cauchy–Schwarz then gave the δ-bound, APN
extremality and the 4-uniform closed form. All exact values (16, 256, 52,
976) were written into the experiment header before the first run.

## Checks

`experiments/experiment_differential_bridge.py`, first run on CPU;
log `out/differential_bridge/run1.log`, exit 0, 12/12. The board-recorded
reproductions run2.log and run3.log are byte-identical to it apart from
run1's trailing exit line.

| check | reference | result |
|---|---|---|
| identity, all 4^n labels | dense PTM from explicit kron Pauli matrices; adder6 from `Circuit.to_unitary()` | 0 mismatches on 11 fixtures (17,024 labels) |
| identity, gate level | `pps.propagate` on 8-qubit Cuccaro adder, X_q, Y_q, X_qX_r, X_qZ_r | 0/100 mismatches; 1..52 terms |
| a = 0 | Walsh sparsity (C8) | agrees |
| APN extremality | Gold x³, inversion at n = 3, 5 | exactly 16 / 256 for every a ≠ 0, b |
| 4-uniform closed form | inversion at n = 4, 6 | exactly 52 / 976 for every a ≠ 0, b |
| bounds and class counts | all fixtures incl. random n = 4, 5, 6 (seed 20260914) | no violation |
| affine | random invertible map + constant, n = 5 | 1 term for every label |

Must-fail controls:
* C1 (Schrödinger-direction σ') disagrees at gold3_gf32 a=1, b=0.
* C2 (naive Walsh count) predicts 1 where the true count is 200.
* C3 (dropped sign) changes the support at a=1, b=1.
* C4 (the unsquared bound 4^n/δ = 512) is violated by the APN minimum 256.

Uniform passes call for scrutiny (METHOD). The dense reference shares with the formula the permutation table, the
U[perm, arange] = 1 construction, the qubit-0-least-significant convention
and the label encoding. adder6's to_unitary route, P2's gate-level
propagator and P10's recorded F10 counts anchor the absolute direction. C1 shows that the dense reference distinguishes the two conjugation
directions. The δ-bound was not vacuous on the S-box fixtures (16 to 256),
but it is vacuous on both adders (δ = 2^n).

## Out-of-sample check on the F10 fixture

F10 never measured this count. The formula predicted it first, as logged in
run1.log: X_x0 → 849,836 terms and Y_x0 → 849,442 terms on the Toffoli
modexp with N=5, a=2, n_exp=1 (14 qubits). The difference between the two
can only come from classes with an even number k ≥ 4 of pairs. Script
`out/differential_bridge/f10_pps_check.py` reruns the unchanged propagator
with a 6M cap. Its first launch failed on an import path before any
computation (`f10_pps_check_importfail.log`).

The dictionary-propagator run (`out/differential_bridge/f10_pps_check.log`)
was **inconclusive**. X_x0 hit the 6,000,000-term cap after 2,639 s with
6,019,010 live intermediate terms. The 3,000 s timeout (exit 124) then
killed the run after Y_x0 had started and before it printed. The Clifford+T decomposition therefore passes
through more than seven times the final support: C99 counts the final
pullback, not the gate-level peak (compare C17). That route had no GPU path,
and choosing it was the mistake.

`experiments/experiment_differential_bridge_scale.py` then settled the
count on GPU (log `out/differential_bridge/scale_run1.log`, 4/4, about
4 minutes on one A4500). `lab/differential.py` was extended for C100 after
that run, so scale_run1 used an unpreserved earlier helper. The run of
record is the board-recorded reproduction `scale_run3.log`
(R8ee3a5298ba44c23), which ran on the archived helper with identical
verdicts and numbers. scale_run2 (R9ff21679dd1a4787) is an earlier
reproduction, identical apart from timings, from an attempt reclaimed after
its lease expired. U and X^a Z^b are built as dense 16384 × 16384
float64 matrices and multiplied with cuBLAS. The product is projected onto
the Pauli basis by a batched Walsh transform over V[c,y] = M[y⊕c, y], so
the route never forms σ, D_c or a DDT; it shares the permutation table and the
conventions listed above with the formula.

* The reference first matched the formula on all 4096 labels of rand6 and of
  adder6 (P9).
* It reproduced F10's recorded diagonal counts, Z_x0 = 3086, Z_x1 = 2926 and
  Z_x0Z_x1 = 2848 (P10).
* It gave X_x0 = 849,836 and Y_x0 = 849,442, with full (c,d) support sets
  equal to the formula's (P11).
* The Schrödinger-direction mutant gives 847,292 terms over 244 X-parts
  against 247, failing as required (C5).

The shared helpers now live in `lab/differential.py`; the first experiment
keeps its inline copies as its record.

## Interpretation and limits

The Pauli basis is poorly suited to off-diagonal observables through
nonlinear permutations. Each DDT class is a set of c-pairs, and a pair's
indicator is dense over half of Walsh space. The same pullback is one signed
permutation, which needs 2^n entries as a table. This explains F10's cap
hit, with the exact final count now known; it does not make simulation faster. The δ-bound is compilation-
independent, like C25, but it was vacuous on both tested adders (δ = 2^n).
Whether that holds for the F10 modexp is unverified. The per-row bound and the exact identity remain informative.

Prior art: brief searches found PTM surveys and DDT/boomerang Walsh
characterizations, but no source combining Pauli pullback counts with the
DDT. The identity is elementary and may well be known, so no novelty is
claimed.

Candidate directions, not opened as TODOs:
1. The same derivation applies to any abelian group. Z/2^t differences in
   the exponent register match the post-QFT observable (TODO14) better than
   XOR differences, and controlled modexp has orbit-structured differences
   there.
2. A monomial-frame hybrid could keep signed-permutation operators
   symbolically between non-permutation gates, as C82 does for linear frames.
   Extracting the expectation value would still carry the cost.
3. Paper A §3.3 could cite C99 to state where the Walsh identity stops.

## Review corrections (V88c6aae3d06240d2)

The first referee review (changes_requested) confirmed the proof and every
logged number, and it required five corrections, all applied above and in
C99:

* **R1:** make a reproduction on the archived helper the run of record
  (now scale_run3, run under attempt Aa544a40df7804b52), because of the
  helper change.
* **R2:** state explicitly that the DDT is taken of π⁻¹. On the PRESENT
  S-box at a = 1 the orientations give 40 against 16. Archived reproductions:
  `out/agent-board/reviews/referee-c99-1-v2/present_check.log` (referee)
  and `out/agent-board/workers/Af1dcc75393f94218/present_orientation.log`
  (coordinator).
* **R3:** say "final count", not the gate-level peak.
* **R4:** scope the vacuous δ-bound to the tested adders.
* **R5:** only an even k ≥ 4 makes the count depend on b.

Its non-blocking notes: the Y_x0 timing wording and the independence
statement are corrected; the run1/run2 exit-line difference is stated. The
`lab/differential.py` docstring now names both claims. The `__main__`
experiment name in run1–run3 log headers comes from the script's own
`Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)` call under `python -m`, not from the
harness, and is left as is. The following remain as stated limitations:
* adder8 bounds were checked on 36 × 4 labels;
* the comparisons use support and |coefficient|, not signs;
* C2 is a weak control;
* `lab/differential.py` carries out-of-scope helpers.

