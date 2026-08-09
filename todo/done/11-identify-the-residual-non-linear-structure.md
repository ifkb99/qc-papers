---
id: 11
state: done
title: Identify the residual non-linear structure
outcome: SOLVED (C33, C34, C35)
claims: [C15, C29, C30, C33, C34, C35]
---
# Identify the residual non-linear structure

**The residue is a conditional linear structure with an exact ¾ density cap.**
See `NOTES.md` §AF and §RS; files `experiment_affstruct.py`,
`experiment_resid1.py`, `experiment_resid2.py`. Density story now complete:
½ = linear structure (C30), ¾ = conditional structure (C34), no forced cap =
two independent nonlinear monomials, 1.000 = structureless.

How the planned sequence resolved:

1. **Affine structures** — counting lemma: density > ½ forbids linear AND
   affine structures, so the circuit half was settled without measurement.
   At function level the finder explained §I's exact-0.500 rows outright:
   all-ones structures (r=6 AFFINE, r=3 LINEAR, r=10 negative control clean),
   PROVED with a three-ingredient criterion (C33).
2. **Support complement** — it IS a recognisable set: the quadrant
   {z_msb=1, z_anc=0} is exactly empty, 0 violations in 3/3 instances (C34).
   Mechanism: t is restored between reductions, so the wraps are inert on the
   t0=0 half-space (baseline structure survives conditionally, proved) and
   degenerate to CNOTs on the t0=1 half (the structure rotates, verified,
   rather than breaks).
3. **Weight profile** — mooted; the structure was identified without it.
4. **Stacking** — decisive and surprising: v4/v5/v45 inject the same monomial
   msb∧t0 and all sit at the same ¾ cap; a second independent monomial (v3x)
   removes the cap (density 0.98). The residue tracks independent nonlinear
   monomials, not wrap count. Bonus, derived from C29 then verified: the
   broken variant keeps C15 constancy exactly (C35) — the ~51% penalty does
   not forfeit the free exponent register.
5. **0.716 the constant** — closed: instance-dependent and width-dependent;
   it is a finite-size snapshot of the approach to ¾ (98.9% of cap at
   n_exp=6, beside the baseline at 98.7% of ½).

Left open (logged in §RS): the live quadrants sit at 0.91–0.98, below the
generic ~0.995 — a smaller, deeper deficit, presumably further conditional
levels. The general 1 − 2^−(k+1) ladder is a conjecture, not claimed.
