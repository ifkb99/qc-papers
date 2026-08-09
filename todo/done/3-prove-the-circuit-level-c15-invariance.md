---
id: 3
state: done
title: Prove the circuit-level C15 invariance
outcome: MECHANISM SOLVED
claims: [C15, C23]
---
# Prove the circuit-level C15 invariance

Paper B's soft centre, now explained. Full write-up in `NOTES.md`
"C15 MECHANISM SOLVED"; files `experiment_c15_proof{,2,3,4}.py`.

**The mechanism (C23).** For β=1 and i ≥ α, every block applies the *same*
permutation V = u_a(·,1), controlled on its own qubit, and **V is an
involution**. So the circuit depends on the entire identity tail through a
single parity bit ⊕_{i≥α} e_i — one effective variable however many qubits it
spans. That explains constancy and liveness *together* (they had looked in
tension), predicts the onset n_exp = α+1, predicts the exact branch relation
h(y) = g(y ⊕ e_prev), and fails for β>1 as required.

Evidence chain, each link verified:
- branch symmetry |ĥ(z)|=|ĝ(z)|: BOTH=0 for β=1, magnitudes match bit-for-bit
  (err 0.00e+00); BOTH=14934/1036508 for β>1;
- sign pattern is a linear character, exact on 100% of support, with v a single
  bit = the previous exponent qubit;
- V² = id exhaustively (2^15 and 2^21);
- flipping any two identity-block controls leaves the function pointwise
  unchanged.

**Two failed routes, do not retry:** affineness of the identity block (14336 of
32768 violations) and the (z, z⊕e) pairing (an n_exp=1 artifact, not
r-dependent).

**FORMALISED — C15 is now a theorem for this circuit family.** Both owed links
closed once the block structure was read properly: `u_a(·,1) = A⁻¹SA` with S a
product of *disjoint* transpositions, so S²=id and V is an involution by
conjugation; the blocks then commute and compose to `V^p`; and averaging the
Walsh character over the tail confines the support to `z_I ∈ {0, 1_I}`, two
values regardless of tail length, so the size carries no dependence on it.

Its sharpest consequence was **derived before being tested** and then confirmed:
0 violations in 7/7 β=1 instances with the two halves individually constant
(7770/7779 at |I| = 2, 3, 4), against 48189 and 1556046 for the β>1 controls.

Scope limit to state in the paper: the proof uses only the
multiply–swap–unmultiply form, so it covers both compilations here and any
Vedral/Beauregard-style construction, but does not automatically transfer to a
modexp built otherwise.
