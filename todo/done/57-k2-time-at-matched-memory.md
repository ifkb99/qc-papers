---
id: 57
state: done
outcome: "C104 Limits + ES: at matched memory the fastest admissible dense-family member is K2's own endpoint (from the proved family identity, not measured); the fastest NON-endpoint member is the table-free route at 7.35-187.55x K2 at the four binding points (measured); no NON-ENDPOINT hybrid was measured-admissible at those points, because of a constant leaf transient of the unwindowed transform, a missing term in the plan's band, and at two of the four points the buffer cap decides both admissibility and the winner (counterfactual)"
title: "At matched peak memory, does the exponent-slice recursion K2 stream the exact modexp Walsh spectrum faster than a multi-pass dense route?"
claims: [C45, C103, C104]
---
# Time at matched memory: K2 against multi-pass dense

Opened 2026-09-17 from TODO 55. For the exact streamed Walsh spectrum of x0 of
`ToffoliModExp(N, a, t).build()` (output contract O2), peak memory alone does
not separate methods: a dense route can split the spectrum by the Walsh bits
of its top k input indices and make 2^k passes over a 2^(q−k) table, so dense
memory can be traded down for time (design review `Vf1d3e39b30b84bb1`, A1; exact at
k = 1, 2, 3 on a 7-qubit random permutation, `out/integrate-todo55/toy_multipass_check.py`). K2's measured footprint
(TODO 55, C104) is therefore a point on a memory–time frontier, not a memory floor.

The question: at the same peak memory, how does K2's time compare with the
fastest k-pass dense route (and with per-slice C45 runs), at fixed N and
growing t? The survey behind TODO 55 derived about 75× less work for K2 than
single-pass dense at N=7, t=6; nothing at matched memory is derived or measured.

Next step (METHOD.md loop): derivation with an experiment plan and design
review before any run. The derivation fixes the cost model (replay element
operations, transforms, emission, which floors every route at Ω(2^q)), what
"matched memory" means, and which dense routes are admissible.

## Outcome (2026-09-17)

Done; the graded statements are in C104's Limits and the round's three mechanism
defects are recorded in note ES. What is NOT answered: anything about dense routes in
general, the family as a whole, the t = 1 side, a windowed K2 or a windowed hybrid
(unbuilt, derived only). A windowed-transform hybrid is the obvious next prototype if
this direction is revived; it needs its own design review, not a re-read of this run.
