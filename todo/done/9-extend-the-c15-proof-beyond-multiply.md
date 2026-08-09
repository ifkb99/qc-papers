---
id: 9
state: done
title: Extend the C15 proof beyond multiply–swap–unmultiply
outcome: DONE
claims: []
---
# Extend the C15 proof beyond multiply–swap–unmultiply

**Confirmed: V² = id is the entire condition.** See `NOTES.md` §G; file
`experiment_c15_general2.py`.

Steps (i)–(iv) of the proof never use V's internals beyond V²=id — "identity on
the valid subspace" was context, not an ingredient. Tested with synthetic blocks
unrelated to modular arithmetic:

- controlled swap `x0↔b0` (order 2): support **constant** at 15362 across
  k=1..4 (and 15248 at N=5);
- controlled 3-cycle `x0→b0→b1` (order 3): support **grows** 15362 → 30984 →
  62088 → 123838;
- vacuity check passes (k=0 gives 3116 vs k=1's 15362), so the blocks act.

Generalised statement: *if a circuit contains k blocks each the same permutation
V controlled on its own qubit, V not modifying those controls, and V²=id, then
the Walsh support of any computational-basis pullback is independent of k.*

Paper B's scope caveat is replaced by a **criterion**: the theorem covers any
modexp construction whose a=1 block is an involution — checkable per
construction. Open follow-up: does windowed / table-lookup arithmetic
(Gidney-style) qualify? Now a well-posed question rather than a survey.

**Protocol catch:** pass 1 was vacuous — blocks acted only on b-qubits while the
observable was Z_x0, so both block types came out constant for a trivial reason.
The tell was the must-fail control failing to fail. Without it, a vacuous test
would have "confirmed" the conjecture for the wrong reason.
