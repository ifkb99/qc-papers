---
code: VR
title: "TODO 66. Clean-scratch virtual reset: exact on clean inputs, collapses the output, but the peak lives inside one controlled multiplication."
outcome: record
claims: [C45, C50, C105, C106]
todo: [66]
---
# VR — Clean-scratch virtual reset (TODO 66, 2026-09-18)

**Origin.** A 2026-09-18 discussion with the user on "the exponential memory of
entangled CNOT gates". The premise was corrected first: a CNOT is one term to one
term in atomic PPS (`perm_pps.py` rule table; C82 removes even the relabel), so
the growth comes from Toffoli branching (C8, C25). The exponent register is already
handled (C45, C104); what remains is the work width q_w = 3n + 4. Three TODOs came
out of it (66, 67, 68); this note is 66. In-session work, no board task and no
referee yet.

## Mechanism

C105's clear rule Z^z ↦ Z^(z & ~mask) is the adjoint of a reset. Where a unitary
circuit returns scratch to |0> on an input set, inserting that clear leaves the
expectation on that set unchanged. `propagate_perm(..., reset_before={k: mask})`
implements it (opt-in; `test_perm_pps.py` [E] checks it against classical replay
with a reset on 24 random circuits, 15 of which the reset changes).

For `ToffoliModExp`, reading `toffoli_arith.py`: t, c0, anc are clean at every
cc_add_mod boundary; b is clean at every u_a block boundary when x at block entry
is < N (the inverse cmult returns b = x_orig to 0 only then). Block 0 sees
x_in XOR 1 because `build()` applies X(x0) first, a slip caught before the run.

## Result (`experiments/experiment_scratch_reset.py`, 46/46, `out/scratch_reset/`)

N ∈ {5, 7}, a = 2, t = 1..3 (one parameter varied). All six predictions and both
must-fail controls behaved as derived:
* R0 (input clear only) equals the WHT of f restricted to clean scratch, exactly.
* R1 (block-boundary clears) equals the WHT of the replay-with-resets function g,
  exactly, and agrees with f on every x < N input. Disagreements with f|S occur only
  at x ≥ N with ≥ 2 exponent bits set (1, 3 at N = 7, t = 2, 3; 3, 10 at N = 5).
* With C45, every block-boundary clear leaves ≤ 2^n terms (8, or 5 at N = 5, t = 1).
* Controls: a clear in the middle of a `b -= N` ladder disagrees on 4–16 x < N
  inputs; R1 differs from f|S at t ≥ 2.

**Exploration (not derived).** Peak retained terms, N = 7, a = 2:

| t | baseline | C45 | R1 + C45 | R2 + C45 | final C45 | final R1 + C45 |
|---|---|---|---|---|---|---|
| 1 | 6522 | 6522 | 6522 | 3936 | 1424 | 8 |
| 2 | 16386 | 8194 | 11443 | 7852 | 1596 | 8 |
| 3 | 48784 | 12197 | 11515 | 8338 | 3876 | 8 |

The output collapses about 500×, but the peak moves only ~1.46× beyond C45. The
peak locator (`out/scratch_reset/peak_locate.py`) shows why: per-segment peaks within
one block, forward order, are [7412, 4397, 1527, 29, 8, 8, 8] at N = 7, t = 2
(cmult's three cc_add_mods, the cswap, the inverse cmult's three). The inverse cmult
is nearly free. The growth happens while reverse propagation climbs back through the
forward cmult, where b is **live** (a partial sum of a·x), so no clean promise
applies. At the peak every b, t, x, c0 and anc sub-key occurs.

## Reading

* Virtual reset removes what C50/C106 Theorem C attributed to dirty inputs from the
  **output**, not from the **peak**. The peak is the Walsh expansion of a live
  modular multiplication in mid-computation, which is the arithmetic nonlinearity
  itself (C25's regime), not scratch junk.
* That points to a block-level Schrödinger/Heisenberg split, not finer clears:
  on clean inputs one block is the permutation (x, ctrl) ↦ (a^ctrl·x mod N, ctrl)
  on n + 1 bits. Its clean-return property can be **certified** by replaying the
  2^(n+1) clean inputs of the block, and its action on an operator supported on x
  computed from that permutation, in memory ~2^(n+1) instead of ~2^(q_w). For
  diagonal permutation circuits this only restates (e, x) ↦ a^e x mod N; its use is
  TODO 64 (Rx insertions between blocks, where X on clean scratch is dropped by the
  same promise) and TODO 14. Standard reduced-subspace simulation; no novelty claim.
* No claim is promoted from this note. P1 is a standard identity and P2 is C105's
  rule on a new circuit; a claim would need referee review first.

## Follow-up (2026-09-21)

The proposed block transfer was tested on a domain-preserving inverse-QFT
sampling contract. C115 owns its validation and negative total-resource
comparison; note MC records discovery and evidence. The work-Rx suggestion
above requires a separate domain-closure argument and is not covered by the
clean-subspace promise. TODO 66 records that residual.
