---
code: LM
title: "TX29 (Pauli-LIMDDs). The post-modexp state in exponent-first orders needs > r/(2 sqrt(2N)) nodes at one level (C110, derived under its hypothesis on a); in x-first orders the reduced diagram has exactly ceil(beta/2) nodes at the first exponent level (C111, derived under its hypotheses). On the measured samples LIM merging buys a constant factor, except when the odd part of r is small."
outcome: record
claims: [C110, C111, C112]
todo: [69]
---
# LM — Pauli-LIMDDs of the post-modexp state (2026-09-19)

**Origin.** An open session the user gave to free research ("use it as you see fit").
A single coordinator session did the derivation and the check. The user then approved
the next steps, and fresh referees reviewed the results; each claim's heading records
what its referees confirmed. TX29 was chosen because it names a real simulator family with no lower bound,
and because the ordered-linear bounds (C107, C108) do not reach it.

**Object choice.** TX29's dictionary speaks of subfunctions of the Boolean bit f. A LIMDD
simulator holds a state, though, so the object here is the state after modexp,
Σ_e |e⟩|a^e mod N⟩. That is a different object from C107's f, and it is what a
gate-by-gate LIMDD run of `build_shor` passes through before the inverse QFT.

**Route.** Prior art first: the knowledge-compilation map (Vinkhuijzen–Coopmans–Laarman,
arXiv:2401.01322, appendix Lemmas 4 and 23) proves LIMDD lower bounds by counting LIM
classes at one level. It says nothing about modexp, periodic or coset states. The
same counting applies here. At the last exponent qubit of an exponent-first order, the
sub-states are two-point states, and LIM equivalence reduces to the XOR invariant
φ(c) = c ⊕ Ac. The fibre bound came from writing c ⊕ d in two ways, as
c + d − 2(c∧d) and as 2(c∧¬d) + d − c, so that c is pinned by a submask of d or of ¬d.
C110 owns the statements.

**What the measurements say.** The derived bound gives up √N, but on random semiprimes
the class count is 0.79–1.00 of the prefix count D_j (⌈r/2⌉ ≤ D_j ≤ r; 255 cases that
meet C110's hypothesis, out of 320 drawn, at n ∈ {10, 12, 14, 16}, j ∈ {0, t − 1}). On
those samples a LIMDD is within a factor 1.27 of the plain QMDD width at that level.
Sweeping every base under the hypothesis takes the ratio down to 0.22, and outside the
hypothesis it collapses (C110, Evidence). A constant lower bound under the hypothesis is
C110's conjecture. Run 2
first reported "0.77–0.96". Its scan used {a^e : e < r} instead of the level's prefix
set, and the upper figure was a median. The referee found both (review Va10c9887fa3f4116, B1). TX29's register text says an escape would need
multiplication by an element of ⟨a⟩ to act as an XOR translation. Lemma 1 makes that
exact: classes = #φ(prefix set), and an escape needs φ to collapse on the level's
prefix set C_j ⊆ ⟨a⟩ (⟨a²⟩ at j = 0 for even r).

**x-first orders (claim C111).** After every x qubit is read, the sub-states are
indicators of {e < 2^t : e ≡ ℓ mod r}. Take r = β·2^α. Complementing e (XOR by 2^t − 1)
maps ℓ to (2^t − 1 − ℓ) mod r. For t ≥ α it induces q ↦ (2^(t−α) − 1 − q) mod β on
q = ⌊ℓ/2^α⌋. For t < α every nonempty set is a singleton.
* Upper bound: at most ⌈β/2⌉ classes (the referee re-derived it).
* Lower bound: in the first version it was measured only. C111 now derives it:
  carry-free pairs force any XOR shift to hold all or none of each window m + supp(β),
  and the Fine–Wilf periodicity lemma then leaves only 0 and all-ones when the set-bit
  positions of β have gcd 1 and t − α ≥ 3·bitlen(β) − 1.
* Measured, and extended by a persistence lemma: C111 owns the range and the statement.

This is the project's β = 1 / β > 1 split (Paper B) appearing again. At the level below
the last x qubit, a LIMDD collapses when the odd part of r is small (at most ⌈β/2⌉
classes). That it does not collapse when β is large is C111's derived statement, whose
hypotheses are given there.

**Reading for the goal.** Both extreme order classes cost about r (or β/2) nodes at one
level. For e-first this is derived up to a √N factor under C110's hypothesis on a: it
fails at every j when 3 | N, and for large j when N has a Fermat-prime factor. For
x-first it is derived (C111) under its hypotheses. So Pauli-LIM merging does not remove
the exponential from modexp states in these orders, within those hypotheses. A small-size
search over interleaved orders found none better than about 0.3·r (TODO 69; exploratory). Interleaved orders
and non-Pauli LIM groups remain open (C110 limits). Every result is a lower bound
(TX15).

**C112 (the interleaved-order count).** The last-exponent level of a *general* order was
first written into TODO 69 with ⟨a⟩ in place of the prefix set C_j, which a referee
refuted with the instance C112 records. The corrected count is
exact, was derived twice independently, and now lives in C112; TODO 69 keeps only the
question it feeds.

**Process.**
* The P1 docstring first said brute force ran to n ≤ 9; the code does n ≤ 7. This was
  corrected after run 1, with the correction noted in the docstring.
* P6 was added after scratch scans had been seen. It is labelled post-hoc and not graded.
* Claim status lives in the claim files, and live review state on the board.
* **What the round cost, and where.** Twelve submissions and twelve scientific reviews
  over three claims: five on C110, two on C111, and five across the three follow-up
  tasks (board history of T4d562bcf23374072, Td6d62de660cb4ba3, Td25d387835d241f1,
  T69a2353121af4979, T642367af3d7549a2). No review refuted a proof. Two referees
  re-derived C110's lemmas, two re-derived C111's, the second of those also re-deriving
  Lemma P, which is new in v2, and two re-derived C112's count. One review did refute a
  statement: the ⟨a⟩ form of the interleaved-order count, which C112 records.
* **Where the defects were.** Mostly in prose, with two exceptions, both in the first
  review: P6a measured over ⟨a⟩ instead of the level's prefix set, a defect in the
  experiment script that forced a fourth run and moved the reported range from 0.77–0.96
  to 0.79–1.00, and the whole sweep ran a single base. The recurring prose defects were a
  hypothesis dropped from a sentence next to the one that fixed it, an extremal number
  quoted without its sample or branch, a fact given a second home, and an exact node
  count asserted of every diagram where only the reduced one earns it. The last had sat
  in C110 since its first version and survived five referees, because the detector
  written for it was pointed at the two newer claims only.
* **The rule that would have saved most of it.** METHOD.md already says a review is a
  correction list, not a source. Version 2 asserted it had re-derived its wording and had
  not: three transcribed sentences each carried a defect, one of them an error the review
  did not contain. Version 3 reintroduced the dropped-hypothesis class in the very
  sentences written to close it. Version 4 asserted a clean detector sweep it had not
  run. Each was caught by the next referee reading the bytes. A detector named in a
  revision is worth nothing until it is run and its output pasted — and this note's own
  first draft claimed every blocking finding was in prose, which the record above refutes.
* **Referee evidence and the lint token.** A referee's own logs tripped
  `EXIT-CONTRADICTS-LOG`, because its must-fail controls print an upper-case FAILED token
  in an exit-0 log. Note TL records the same collision twice. The referee prepared
  reworded copies, its host refused to run them, and it left its attempt blocked with
  two options for the coordinator; the coordinator declined to re-run them and
  cancelled the task, citing the three recorded runs the review rests on. The fix belongs in
  `tools/evidence_lint.py` or in the harness's lower-case token, not in reworded evidence.
