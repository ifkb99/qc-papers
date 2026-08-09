---
id: 12
state: done
title: Windowed / table-lookup arithmetic
outcome: DONE, and it corrected C29
claims: [C8, C24, C29, C30, C32, C36, C37, C38, C39]
---
# Windowed / table-lookup arithmetic

**The criterion was incomplete, and this is the item that found it.** See
`NOTES.md` §WD; files `windowed_arith.py`, `test_windowed.py`,
`experiments/experiment_windowed.py`. Claims C36–C39; C29 regraded to NARROWED.

The literal question has a yes answer that turns out not to decide anything:
**both** windowed designs satisfy V² = id exhaustively, and one of them loses
the invariance completely. So V²=id does not discriminate.

- **C36, the repaired criterion.** C29's *other* hypothesis — each block
  controlled on **its own qubit** — is load-bearing. The general condition is
  that the identity-tail block's dependence on the exponent register be
  **affine**: identity (one fresh control) or constant (uncontrolled) qualify,
  OR does not. Step 9 could not see this because every block it tested had
  exactly one control, where the two hypotheses coincide.
- **C37, lookup (Gidney-style) is covered, and more strongly than expected.**
  The tail table is all-ones, so the QROM permutation is `s ^= 1` regardless of
  j and the block never reads its window. The tail is **dead**, not merely
  constant: support confined to exponent bits {0…α−1}, verified 3/3 with α =
  1,1,2. Directly distinguishable from C24, which says the standard
  construction's support *contains* the all-ones-tail vector — same modulus and
  base, opposite geometry.
- **C38, the honest caveat.** Cost is **bounded and 2-periodic** in the
  tail-window count, not the exact constant C24 gives, because the tail applies
  W unconditionally and W²=id. 28078/62680/28078/62680 at K=0..3; K=0 and K=2
  support **sets identical**. Control β=3 grows.
- **C39, the cause is not windowing.** At w=1 — no window at all — the lookup
  form still has a dead tail. The difference is **emitting the multiply-by-1
  branch instead of optimising it away**. Skipping it (`SelectModExp`) grows
  **exactly ×4.00 = 2^w per window**, the derived rate (OR of w bits has Walsh
  sparsity 2^w), with density pinned at C30's ½ cap — i.e. a β=1 instance is
  put onto the generic β>1 curve. Mirror image of C32: there the standard
  reduction was accidentally cheap, here standard practice would be
  accidentally expensive.

**A prediction was refuted and the refutation was informative.** P5 (restoring
the j=0 branch restores the dead tail) was a conjunction: the flatness held
(support ~32k at every K, no growth) but the dead-bit half failed — the live
window's other bit stays live. Diagnosed, not waved away: the **activation
ancilla** is itself a scratch qubit the pullback ranges over; the lookup design
cancels it (2^w CNOT contributions into s XOR an even number of times), the
select design does not (the ancilla gates a block). Confirmed by checking the
block on the act=0 half-space, where the e1 dependence vanishes.

**VALIDATED AT SOURCE** (arXiv:1905.07682 §3.5, arXiv:1905.09749). C39's
mechanism is Gidney's own stated rationale — *"this also removes the need for
the multiplications to be controlled, because the table lookup can evaluate to
the factor 1 in cases where none of the exponent qubits are set"* — so C37–C39
describe the construction people actually propose. The joint `table[ei, mi]`
indexing we simplified away is confirmed harmless (in a tail window every
ke = 1, so the outer index goes degenerate). His *relabelling* swap, which we
modelled as a physical swap, turns out to be what makes the block an involution
at all, and his `if a is not target: swap(a, b)` line is C38's 2-periodicity
appearing in the real compiler. **Remaining real gap:** the analysis assumes a
**unitary** unlookup; Gidney's measurement-based uncomputation is not unitary,
so C8 does not apply to it as written.
