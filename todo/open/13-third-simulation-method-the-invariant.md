---
id: 13
state: open
title: Third simulation method on the r = β·2^α invariant
claims: [C103]
---
# Third simulation method on the r = β·2^α invariant

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
