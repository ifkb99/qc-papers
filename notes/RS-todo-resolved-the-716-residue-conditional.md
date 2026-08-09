---
code: RS
title: "TODO 11 RESOLVED: the 0.716 residue is a CONDITIONAL linear structure"
date: 2026-08-08
outcome: solved
claims: [C15, C29, C32, C34, C35, C40, C41]
todo: [11]
---
# RS — TODO 11 RESOLVED: the 0.716 residue is a CONDITIONAL linear structure

`experiment_resid1.py`, `experiment_resid2.py`; predictions in the file
headers, written before measuring.

### The finding

> Injecting a single nonlinear monomial msb∧t0 does not destroy the linear
> structure — it **demotes it to a conditional one**. The support of every
> broken variant avoids the quarter **{z_msb = 1, z_anc = 0}** exactly
> (0 violations, 3/3 instances), capping density at **¾**. The measured
> 0.716 → 0.742 trajectory is the approach to that cap: at n_exp=6 the
> variant sits at 98.9% of ¾ where the baseline sits at 98.7% of its ½ cap.

Physical reading: **a Z on the accumulator's sign bit never appears without a
Z on the comparison ancilla.**

### Mechanism, read from the construction

**(i) Inert half-space.** The t register is restored between reductions
(Cuccaro MAJ/UMA restores the operand; loads are unloaded), so every wrap
sees the INITIAL t. On {y : t0(y)=0} the v4/v5/v45 wraps are literally the
identity — the circuit IS the baseline there — and L2's structure
w = msb⊕anc survives conditionally: g(y⊕w) = g(y) for all y with t0=0.
Proved by construction; verified pointwise, 0 violations (7584 on the
must-fail t0=1 control; the t0=0 slice equals the v0 slice bit-for-bit).

**(ii) Linear firing branch.** On t0=1 the wrap toffoli(msb,t0,·) degenerates
to cnot(msb,·) — a LINEAR gate — and XOR couplings cannot break an XOR
symmetry (L2 pass 1); they *rotate* it. The t0=1 slice carries its own exact
linear structure w′: **v4: msb alone; v5: msb⊕t1⊕anc**. Verified pointwise,
exhaustive, at N=5 and N=7.

**(iii) Walsh consequence.** For <z,w>=1 the t0=0 half cancels pairwise, so
Ĝ(z) is carried by the t0=1 slice alone, and the slice structure kills its
odd-coset spectrum. Chained: support ⊆ {<z,w>=0} ∪ {<z,w>=1 ∧ <z,w′>=0} —
one quadrant exactly empty, density ≤ ¾.

```
 v4 support density by (z_msb, z_anc) quadrant:
              (0,0)    (0,1)    (1,0)    (1,1)
  N=5 a=2    0.9099   0.9805   0.0000   0.9739     0 support elements in (1,0)
  N=7 a=3    0.9159   0.9802   0.0000   0.9822     0
  N=7 a=6    0.9136   0.9771   0.0000   0.9766     0
```

### Why one monomial ≠ two wraps ≠ two monomials

v4, v5 and v45 all inject the **same** monomial msb∧t0 (v45's two wraps share
it): same inert half-space, same cap — all sit at 0.716–0.722. v3x adds a
second independent monomial msb∧t1: the wraps are then inert only on a
QUARTER space, a quarter-space invariance forces no exact zeros at all, and
density jumps to 0.980–0.983. The residue was never "intrinsic to the
arithmetic" — it tracks the number of independent nonlinear monomials
coupling the msb, not the number of wraps.

### Independence from C15 (derived from C29, then verified)

The wrapped a=1 block is A′⁻¹SA′ — still a conjugate of the involutive cswap
layer — so C29 applies verbatim and constancy survives the breaking: v4 at
N=7 a=6 has support exactly **23488 at n_exp = 2, 3, 4**. The ~51% penalty
does not forfeit the free exponent register. Claim C35.

### The ladder, and corrections to C32's framing

> **CORRECTED 2026-08-08 by §GF (C40/C41).** The ladder below is the wrong
> parametrisation. The cap is **1 − 2^−d with d = dim span{per-cell
> structures}**, plus a **parity condition** on dependencies among them — not
> a function of the conditioning depth k. Two cells *sharing* a structure give
> d = 1 and a ½ cap, not ¾; and an odd dependency removes the cap entirely at
> unchanged d. The k = 0, 1 rows below are right for the reason §GF gives
> (d = 1 and d = 2), and the "two independent monomials → no forced cap" row is
> most likely the parity case, not a larger d. Also: the linear→hyperplane step
> is Carlet's **Proposition 29**, classical — cite it.

baseline: cap ½ (proved, L2) → one monomial: cap ¾ (this section) → two
independent monomials: no forced cap (0.98 measured). The pattern
1 − 2^−(k+1) for codim-k inert subspaces is a natural conjecture but is NOT
claimed; only k = 0, 1 are established.

C32's verb was wrong: a single nonlinear coupling **halves the forced-zero
set** rather than removing it. The "~51%" is a finite-size snapshot of a
¾-cap function beside a ½-cap function (asymptotic ratio 3/2 exactly). And
"0.716" is not a constant: it is instance-dependent (23464 / 23488 / 23579
per 32768) and width-dependent (rising toward ¾). TODO 11.5 closed.

### Open residue, honestly

The three live quadrants sit at 0.91–0.98 — below the generic ~0.995 zero
rate — with the (0,0) quadrant lowest. A further, smaller deficit exists,
presumably deeper conditional levels inherited from the baseline slice. Not
pursued.

**Grades.** Conditional survival on the inert half: **PROVED + verified**.
Quadrant law and ¾ cap: **verified 3/3 instances; mechanism proved except one
link** — the existence of the rotated slice structure w′ is exhaustively
verified per instance and argued from wrap-linearity, not yet proved in
general. Claim C34.
