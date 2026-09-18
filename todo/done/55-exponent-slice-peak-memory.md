---
id: 55
state: done
outcome: "C104/ES: the exponent-slice lemma, per-slice runs B1 and the dense recursion K2 are proved exact; K2 measured at (N=7,a=2) t<=10, (N=11,a=2) t<=6 and (N=7,a=6) t<=6, exact against the engine at the 14 phase-A points, with peak-byte ratios against the engine only at F1 t<=5, F2 t<=2 and F3 (F1 t=6 and F2 t=3 are not supported); applied known mathematics; no memory crossover with dense exact streaming (multi-pass dense trades time for memory), so time at matched memory moves to TODO 57"
title: "Does slicing the full-space modexp Heisenberg propagation by exponent Walsh bits remove the peak above a streamed exact output?"
claims: [C45, C82, C101, C103, C104]
---
# Exponent-slice decomposition of the full-space peak

Chosen from the SL slate (2026-09-16) under output contract **O2**: the exact
final operator of `perm_pps.propagate_perm` on `ToffoliModExp.build()` may be
delivered as a stream (each term exactly once, keys and values exact) rather
than held as one dictionary. Under O1 (materialized dictionary) no method gains
more than about 1.5× in term count, because the final operator has about
2^(q−1) Walsh terms; see SL.

## The candidate

Exponent qubits only ever occupy Toffoli control slots of their own block, so
after the blocks above block k are processed, the stored terms split into
slices keyed by the processed exponent Walsh bits, and later gates never mix
slices. A depth-first walk over those bits performs the same per-term updates
as the existing loop and holds at most 2^(q_w+1) + (t−1)·2^(q_w) terms beyond
the output stream [term-count proxy, derived in SL's deriver slate, unreviewed].
Each slice is expected to equal a C45 reduced observable
(`propagate_perm(zmask=obs|u, trace_plus=exponent)`).

Checked so far (diagnostics, not evidence for a claim): slice sums equal the
engine's `n_terms` at every step on six circuits at t ≤ 3; term-count proxy
peak 15,217 against 48,784 at N=7, a=2, t=3.

## Next step (METHOD.md loop)

Derivation with its experiment plan, then a design review by a fresh referee,
before any run. The derivation must settle: the exact slice lemma including
X gates and the lazy frame; the peak bound in terms and in a stated byte model;
total work against the existing loop; the output-stream contract; the relation
to C45 and to `lab/semiclassical.py`; and which baseline (the existing engine,
the dense work-state recursion, per-slice C45 runs) it must beat.

## Outcome (2026-09-17)

Done; see C104 and note ES. The question is answered for the dictionary engine and
left open against dense routes, where it becomes TODO 57 (time at matched memory).
