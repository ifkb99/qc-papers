---
id: 13
state: open
title: Third simulation method on the r = β·2^α invariant
claims: [C103]
---
# Third simulation method on the r = β·2^α invariant

**Triage (2026-09-22, METHOD.md barrier check 4):** success measures a real DD simulator's memory at small sizes and computes nothing known to be hard. Passes.

Paper B open problem 4. Two structurally unrelated methods (PPS, MPS) keying
on the same arithmetic invariant is suggestive; a third would make "property
of the algorithm, not the simulator" hard to argue with. Decision diagrams are
the natural candidate (MQT DDSIM): rerun the controlled design (fix N, vary a
so only r changes) and see whether memory keys on β. Mostly integration work,
no new theory.

C102/HD (2026-09-15) supply a related same-object baseline: final ROBDD
sizes of the full-space u_a pullback against its Walsh support, for one
fixed order and single controlled multiplications only. They do not run
the β-controlled modexp design this item asks for.

C103/HB (2026-09-15) do run the β-controlled design, and this item stays
**open**. What C103 settles: the exponent-first ROBDD of C15's own full-space
object has widths bounded by 2^α D_{μ(β)}(k − α), so a third store does key on
r — but through α and μ(β), not through β alone, and only asymptotically in t,
because below a crossover the bound is the trivial 2^k (C103's Limits own the
exact threshold and the trivial-level count). What it does not settle:
this is a canonical decision diagram of the pullback, not DDSIM's state
decision diagram, so the memory-keys-on-β question as posed here — rerun the
controlled design in a real DD simulator and watch its memory — remains
untouched. C103's Limits say the same.

DN (2026-09-15/17) closes the propagation-native branch of this item negatively and
**leaves it open** for its original target. A decision diagram of the propagated
operator (frozen pilot dd-pilot-20260915) was built, validated and measured in two
accepted board rounds: about 2x fewer nodes than a matched null, but 11.7-40.7x more
bytes than the dictionary route peak to peak. Its one arithmetic-looking feature, a
turn in the peak-node ratio at n_exp = 4 for N=7 a=3, is NOT C103's crossover:
α + μ(β) + 2 predicts 3 / 4 / 3 / 4 for a = 2, 3, 4, 5 and the coordinator's own
measurements (unreviewed; see DN's table and its labelling) give 5 / 4 / 5 / 4,
reversing the order it predicts. The promotion round separately showed, under review,
that C103's exponent-first order does not transfer to this
store (0.920 against 0.518 at n_exp = 2; see DN). So this store does not key on r through
α and μ(β) the way C103's widths do. Do not re-run the ordering design; DN records the
bases, the inverse-base confound that blocks the α reading, and the degenerate N=5
family. What remains untouched is what the item actually asks for: a real DD simulator's
state decision diagram (MQT DDSIM), which neither C102/C103 nor DN measures.

Reproducibility gap: the pilot propagator, the proposed `lab` module and the proposed
experiment all remain under gitignored `out/`, so nothing in git reproduces DN's
numbers. TODO 59 owns that promotion.
