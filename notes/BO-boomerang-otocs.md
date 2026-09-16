---
code: BO
date: 2026-09-14
title: "OTOCs of reversible arithmetic as boomerang counts"
outcome: X-X OTOC = BCT identity verified; exact Cuccaro carry butterfly; compiled tail-bit OTOC symmetry explained by block inverses (beta in {1,3} when t-2 >= alpha) after a must-fail control did not fail; ideal-model sweep found non-discriminating
claims: [C101]
todo: [53]
---
# BO — Scrambling in reversible arithmetic, through the boomerang table

C101 owns the statements. This note records the path to them.

## Why

After C99 and C100, the user asked for one more area of Claude's choice. A
repository search found no OTOC or scrambling work: every earlier "OTOC"
match was the word "protocol". The OTOC is the standard operator-scrambling
diagnostic and is measured on hardware. Tracing V, W(t), V†, W(t)† through
a permutation gives the boomerang condition of Cid et al., the third
cryptanalytic table after the linear (C25) and difference (C99) tables. It
also suggested a probe of the r = β·2^α invariant, in the spirit of TODO13.

Primary source, body read: Boura & Canteaut, *On the Boomerang Uniformity of
Cryptographic Sboxes*, ToSC 2018(3):
* Definition 3 (p.295) matches the derived condition, with π in both places;
* Prop. 2 (p.296): the BCT of π⁻¹ is the transpose, which makes the
  orientation control C1 meaningful;
* Prop. 6 (p.303): inverse-map values;
* Table 5 (p.302): PRESENT.

Searches for a BCT–OTOC link, and for exact OTOC counting in classical
reversible circuits, found none.

## Run 1: `experiment_boomerang_otoc` (10/11, exit 1 on purpose)

Everything derived beforehand passed: the identity, the gate-level unitary
check, the Boura–Canteaut values, the butterfly (0, ½, ¾, ⅞ on m = 4), and
the ideal-model F = 1 at tail bits. The C99/C100 helper `power_map` was
reused; the OTOC code in `lab/otoc.py` is new.

P6 predicted identical tail rows for compiled β = 1, citing the C23/C46
identical-involution mechanism. Its must-fail control C5 used β = 3
(N = 7, g = 3) and did not fail: the rows were identical there too. The
file is left failing.

## Reading the construction again

Flipping e_k turns the boomerang condition into the commutation of M_k with
the later-block conjugate of the XOR by q. Compiled blocks satisfy
M_{c⁻¹} = M_c⁻¹ as permutations. The gate lists are not literal reverses
(the controlled swaps appear in the opposite order), but those swaps commute.
For N = 7 the constants are 3, 2, 4, 2: 4 and 2 are inverses. So equality
follows whenever c_{t−2}³ = 1 or c_{t−2} = 1, i.e. β ∈ {1, 3} when t−2 ≥ α.
My control had picked the other β covered by this mechanism.

## Run 2: `experiment_boomerang_otoc_tail` (6 checks plus an always-passing record)

Predictions were registered after that derivation:
* β = 3 compiled (N = 9, 13): tail rows equal.
* β = 5 compiled (N = 11, g = 2 and 3): rows differ, which serves as the
  control.
* The block-inverse identity holds, and its mismatched pair fails.
* The ideal-model sweep had no violations.

That last result was a tell. The open count Q4 showed that 124 of 136 ideal
cases *without* c_{t−2}³ = 1 also have equal rows. A post-hoc diagnostic
(`out/boomerang_otoc/q4_zero_rows.py` and its log) showed that 24 of those
124 rows are all zero. The other 100 are identical across every g at fixed N
(for example N = 23: ¼, ¼, ¼, 1/16, 0). My first reading attributed them to
the fraction of v with v and v ⊕ 2^j both in the padded region w ≥ N.
Review V7cf842d53f784f58 showed that this misses N = 22 (j = 4) and N = 25
(j = 3, 4). There, residues fixed by the tail multipliers pair with padded
points: (11, 27) for N = 22 and (20, 28), (10, 26), (15, 31) for N = 25.
The first entry of each pair is a fixed residue and the second its padded
partner. Post hoc, the rows look determined by points the multipliers fix.
In the ideal model the tail equality is mostly such a fixed-point artifact,
and Q3 passes for that reason. The compiled circuits are what discriminate.

Provenance: by file time, tail_run1 and tail_run2 ran the tail source
revision cddcdc7c, which differs from the archived 5a7838f1 only in its
STATUS docstring (V16ef0281b61e45a1 N1).

## What it means

* **Scrambling reads directly off the table.** OTOCs of permutation
  circuits are exact combinatorial counts, and every published BCT becomes
  a scrambling table. PRESENT has a nonlinear direction pair that does not
  scramble at all.
* **Adders scramble only locally.** Scrambling is total on the diagonal
  and decays as 2^(i−j) along the carry chain.
* **The invariant appears only as symmetry.** In compiled modexp the
  r = β·2^α structure shows up in OTOCs through tail-bit symmetry, not as
  F = 1, and it cannot tell β = 1 from β = 3.

Candidate directions, not opened: signed OTOCs (d, e ≠ 0) as a
Walsh-weighted BCT, and whether post-QFT OTOCs have a clock-shift (C100)
boomerang analogue.

## Review corrections (V7cf842d53f784f58)

The first referee review (changes_requested) confirmed the identity, the
Boura–Canteaut citations, the sufficient direction of the tail symmetry and
the chronology. It required seven corrections:

* **R1:** literal predict/must_fail registration, a new recorded run and a
  clean lint.
* **R2:** a carry-butterfly proof that handles j = m−1 through the carry-out.
* **R3:** "as permutations", not "gate for gate".
* **R4:** "if", not "iff", and state t−2 ≥ α.
* **R5:** restate the padding explanation with fixed residues and zero rows.
* **R6:** the orientation anchor is PRESENT, not adder6.
* **R7:** minor wording and counts.

The coordinator reproduced R2 (128, 256, 512 and 1024 carry-out-only
failures at m = 4), R3 (the gate lists are not reversed) and R6 (adder6 X–X
OTOCs are orientation-blind) before editing. The experiment docstrings keep
their original pre-registered text, and their STATUS headers carry the
corrections.

