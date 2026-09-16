---
code: HD
date: 2026-09-15
title: "Heisenberg decision diagrams: the PPS object in a canonical Boolean store"
outcome: adder Walsh/ROBDD/affine counts proved and reproduced; u_a Walsh ~8^n against ROBDD ~4^n measured to n=8 (29 qubits) against a slowly growing random-circuit control; one registered prediction (P5a) failed at n=4 and is preserved; no bound for u_a, no escape from the order barrier
claims: [C102]
todo: [54]
---
# HD — Storing the pulled-back observable as a decision diagram

C102 owns the statements. This note records how they were found and checked.

## Why

The user invited Claude to pick a direction on "the exponential growth of
memory related to CNOT gates", and later authorized spending the remaining
usage on it, with the arb board kept current (task Ta294b7d88a4c4015). C82
already showed that CNOTs cost PPS nothing beyond a frame update. The
support growth comes from Toffoli branching, and C8 identifies the final
support with Walsh sparsity. So the question became whether the Pauli
(Walsh) basis is a good store for this Boolean function at all.

## Path, including dead ends (scratch pilots, not evidence)

The pilot scripts and logs are archived in
`out/agent-board/workers/A46391e1cf2954646/pilots/`. Board message
M741d61dd419947ef discloses every number seen before registration.

1. **Signed affine-subspace indicators.** A Toffoli splits a piece into a
   signed three-term sum. On the adder this gave 2^m − 1 pieces, only about
   1.5× below Walsh. Branching on one control instead (y_u = 0 kills the
   product, y_u = 1 turns it into y_w) gave 2m pieces or 2^m, depending on
   which control is branched.
2. **The same rule on `cc_add_mod`.** Peaks were 2–2.3× below Walsh, but
   final counts were above it (456 vs 86). The representation is not
   canonical, so compute–uncompute cancellations were invisible. Sibling
   merging (two cosets of the same subspace with equal coefficient) fixed
   the final counts: u_a N=7 and N=11 gave 732 and 6081 pieces against
   Walsh 3206 and 31022. This is a constant factor, and the merging is ad hoc.
3. **Canonical ROBDD via `dd`.** With a good order it beat the affine pieces
   on u_a. In register-block order it was exponential on the adder, the
   classic XOR/ordering blowup.
4. **A scratch scaling run.** It suggested the ratio was growing, but it
   used complement-edge counts and dd peaks. A scratch prediction for n=6
   (log2 peak in [16.2, 16.9]) was refuted at 15.64–15.89. A scratch n=7 run
   was killed before printing, so n = 7 stayed unobserved.
5. **A same-multiset random-circuit control.** It kept the ratio near its
   generic value (the ROBDD of a random function on q variables is about
   2^q/q). That separated "BDDs are compact" from "arithmetic is special".

Coordinator speculation in chat, that the ROBDD was heading towards table
scale (about N), was not registered and is not supported. The all-N n = 8
step appeared to rise again, but review V06b5088b79994f1c showed this was a
class-mixing artefact (see below). The asymptotics remain unresolved.

## Checks

The reducer `lab/bdd_count.py` counts ROBDD nodes from a truth table by
bottom-up hash-consing, and never touches gates. Before any registered run
it agreed with a brute-force distinct-subfunction count on 300 random tables
(n ≤ 5, random orders, both conventions). That session check was not
archived at the time and cannot now be recovered. A reconstruction of it,
written after review V06b5088b79994f1c, is archived as
`out/agent-board/workers/A46391e1cf2954646/bruteforce_reducer_v1.py`
(600 comparisons, 0 mismatches; the first attempt failed on an import path
and its log is kept beside it). It is the same procedure, not the same run:
identity with the unarchived session check is not verifiable. The `dd` cross-check covers only the
complement convention: `dd` counts complemented edges, terminal included. The truth tables come from the
existing verified permutation replay: `walsh.classical_permutation`, or
`accel.classical_permutation` for q ≥ 20, which `test_accel.py` gates.
Walsh counts use `walsh.pullback_coefficients` below q = 20 and the exact
integer GPU path above it. The `dd` composition is an independent canonical
route that goes gate by gate and never builds a truth table.

**`experiments/experiment_heisenberg_bdd.py`**
* Board run R2fe53ff3836846c2, log `out/heisenberg_bdd/run_v1.log`: exit 139,
  a native crash with no traceback during u_a n = 4, after the adder section
  had passed (CLAUDE.md trap 3, CPython 3.14). Retained.
* Board run Rdb87a87850f64ccc, log `out/heisenberg_bdd/run_v2.log`, the same
  source under Python 3.12.10: exit 1 by design, 13/14. Its report is
  `out/heisenberg_bdd/report_v1.json`; run_v1 crashed before writing a report.
  * P1–P4 pass (derived; P1–P3 had already been evaluated for m ≤ 9 while
    the reducer was being developed, as disclosed in M956c9a970d43415a).
  * P5a FAILED: Walsh step 3.375 at n = 4, outside [2.8, 3.2].
  * P5b, P5c, P6, P7 and P8 pass.
  * Controls C1 (random ratio growth R_6/R_3 = 2.02 < 3; u_a 8.05), C2
    (register-block order at m = 10: 3585 nodes), C3 (quasi-reduced mutant)
    and C4 (adder without the first MAJ/UMA Toffoli) all failed as required.

**`experiments/experiment_heisenberg_bdd_n8.py`**
* Board run Rb9caf0557f6e4b48, log `out/heisenberg_bdd/n8_run_v1.log`,
  report `out/heisenberg_bdd/n8_report_v1.json`: exit 0, 5/5. Predictions
  were registered in M2d005ce2efc84ed9 before any n = 8 value existed.
  * Q1: ROBDD step 1.804, in [1.0, 2.0].
  * Q2: Walsh step 3.005.
  * Q3: R_8/R_7 = 2.299.
  * Q4: all eight n = 7 rows reproduced run_v2 exactly.
  * C1 failed as required: random R_7/R_6 = 1.108.

Each 29-qubit fixture took 222–273 s of wall time (harness row times),
probably dominated by the host-side 2^29-entry reduction (an estimate; not
profiled). The core gate (`out/heisenberg_bdd/core_v1.log`) passed before
run_v2 and the n8 run. Its timestamp is 6.5 s after run_v1's run.start, so
"before execution" is not established for run_v1, which crashed and is not
evidence. Other science
suites were not rerun, because no existing module changed.

## Weak points a reviewer should press

* **C1's construction and margins.** The random control keeps u_a's numbers
  of X, CNOT and Toffoli gates but draws qubits uniformly, one seed per n; it
  is not u_a's exact gate multiset. Margins are 2.02 against a threshold of 3,
  and 1.108 against 2, the latter comparing random R_7/R_6 with u_a R_8/R_7. The
  control's own ratio grows slowly, as expected for random functions
  (their ROBDD is near 2^q/q).
* **Medians over changing N, and an N mod 4 confound.** Found by review
  V06b5088b79994f1c. B separates by N mod 4 (the N ≡ 3 class is larger at
  n = 5, 6, 8), all eight n = 7 fixtures are N ≡ 1 (mod 8), and n = 8 has two
  N per class. Within N ≡ 1 the ROBDD log2 steps are 2.220, 2.117, 1.979,
  1.762, 1.643 (`class_split_v1.log`); the all-N rise from 1.629 to 1.804 was
  class mixing. At n = 8 B spans 512,765–683,089 (max/min 1.33).
* **One order for u_a, chosen from pilots.** No reordering search was done.
* **Peaks.** Only n ≤ 5 (ROBDD, complement convention, nodes of the
  intermediate function only) and n ≤ 4 (PPS). Everything at n ≥ 6 is a
  final count.
* **The affine-piece counts depend on the branching rule.** They are not a
  function invariant.
* **No primary-source audit** of decision-diagram quantum simulators (Zulehner
  and Wille, DDSIM) was done beyond confirming that they exist. Bryant 1992
  was read for the adder statements (§1.3–1.4, Table 1, Fig. 4). Its
  bounded-cross-section bound (Berman; McMillan) is a possible route to a
  proved u_a upper bound, not attempted.

## Candidate directions, not opened

* A proved upper bound on the u_a ROBDD from the compiled macro structure,
  using cross sections or residue-state streaming, which would settle
  H-table against H-square.
* A dynamic-order or linearly transformed BDD (a C82-style frame that makes
  CNOTs free) on the full ToffoliModExp, compared with Paper B's β
  dichotomy (TODO13).
* Gate-level peaks at n ≥ 6 with a compiled BDD package.

## Review corrections (V06b5088b79994f1c)

The first referee review (changes_requested) confirmed the adder proofs line
by line, recomputed every u_a row with q ≤ 20 independently, and verified
the chronology and the Bryant citation. It required eight corrections, all
applied above, in C102 and on the task thread:
* **R1:** clear the lint ERROR on the smoke helper (v2 omits it).
* **R2:** disclose the N mod 4 split.
* **R3:** correct the B/4^n range to 9.7–14.3.
* **R4:** make the ranges exact.
* **R5:** archive the brute-force check and state P6's convention.
* **R6:** describe the random control accurately.
* **R7:** state the dd peak convention.
* **R8:** record the n = 8 budget and write-scope deviation.

Non-blocking notes applied:
* the timing statement is labelled an estimate;
* the core-gate chronology is corrected;
* the Bryant wording now notes the free carry-in;
* C102 names the complement merges and the highest non-propagating
  position, and says the Walsh side is near-dense.

Also recorded:
* run_v1's control_result says run_v2 used Python 3.13, while it used
  3.12.10 (thread note).
* No mutant exercises the complement branch or the dd cross-check.
* Rows with q ≥ 23 have no independent recomputation.
