---
code: GF
title: "THE RECURRING GF(2) PATTERN: one fact, two corollaries, one folklore gap"
outcome: record
claims: [C15, C24, C30, C32, C33, C34, C36, C40, C41]
todo: []
---
# GF — THE RECURRING GF(2) PATTERN: one fact, two corollaries, one folklore gap

`experiments/experiment_gf2law.py`. Prior art checked at source first; the
derivation was done before any measurement; 7/7 checks pass with two must-fail
controls failing as required.

### The question

Three results here have the same shape — affine ingredients are free, nonlinear
ones are expensive: §L2/C30–C32 (linear structure ⟹ ½ cap; a linear wrap is
harmless, a Toffoli destroys it), §RS/C34 (one nonlinear monomial demotes it to
a *conditional* structure, ¾ cap), §WD/C36 (affine block control keeps C15,
OR does not). Is that one theorem?

### The answer: one FACT, two corollaries — and one of them was already textbook

The fact is the defining property of the transform: **the Walsh characters of
GF(2)ⁿ are exactly the affine functions**, so `(−1)^affine` is a *single*
coefficient and anything nonlinear is spread. Two distinct corollaries:

**(a) Affine symmetry confines the support.** This is **classical and we must
cite it, not claim it.** Carlet, *Boolean Functions for Cryptography and Coding
Theory*, **Proposition 29**: *"The derivative D_e f equals the null function
(resp. function 1) if and only if supp(W_f) is included in {0ⁿ, e}^⊥ (resp. in
its complement)."* That is C30 (linear structure ⟹ hyperplane ⟹ density ≤ ½)
and C33 (affine structure ⟹ the opposite coset) verbatim. Functions whose Walsh
support is an affine subspace are Carlet's **partially bent** functions.

**(b) Affine gating costs one Walsh coefficient.** If a circuit's dependence on
a register runs through a gating function φ, that register contributes a factor
`1 + |supp(φ̂)|` to the support. φ affine ⟹ φ̂ is a single character ⟹ factor 2
(C24's `z_I ∈ {0, 1_I}`); φ = OR of w bits ⟹ 2^w coefficients ⟹ factor 2^w per
window. §WD measured **exactly ×4.00 at w = 2**. Same fact, different place.

### C40 — the conditional law (the piece with no citation found)

Derived before measuring. Let a coset partition split F₂ⁿ into cells H_u, and
suppose f restricted to H_u has linear structure w_u. Decomposing the transform
over the partition, `c_z = 2⁻ⁿ Σ_u (−1)^{u·z_C} A_u(z')` and Proposition 29
kills `A_u(z')` whenever `w_u·z' = 1`. Hence

> **supp(f̂) ∩ E = ∅ for E = {z : w_u·z = 1 for every u}**, so
> **density ≤ 1 − 2^−d with d = dim span{w_u}** — when that system is
> consistent.

```
  planted, n=14                cells   d  |E|    cap     density  violations
    no conditioning (C30 form)     1   1  8192  0.5000   0.4919      0
    2 cells SHARING one w          2   1  8192  0.5000   0.5000      0
    2 cells, independent w         2   2  4096  0.7500   0.7466      0
    4 cells, independent w         4   3  2048  0.8750   0.8649      0
```

**This corrects §RS.** The conjectured "1 − 2^−(k+1) ladder" is the wrong
parametrisation: the cap is set by the **span dimension of the per-cell
structures**, not by the conditioning depth. Two cells *sharing* a structure
give d = 1 and a ½ cap (measured 0.5000), not ¾. The ladder is only the special
case where each new level contributes one new independent vector.

### C41 — consistency is a PARITY condition, and it can destroy the cap

The system `{w_u·z = 1}` is solvable iff every linear dependency among the w_u
has **even** support. An odd dependency (`w_0 ⊕ w_1 ⊕ w_2 = 0` forces
`1⊕1⊕1 = 1 ≠ 0`) makes E empty and the cap vanishes entirely.

```
  4 cells, ODD dependency    d=2  consistent=False  |E|=0     density 0.9845
  4 cells, EVEN dependency   d=3  consistent=True   |E|=2048  density 0.8671  (cap 0.8750)
```

Same cell count, same construction, caps 1.0 vs 0.875 **purely from the parity
of a dependency**. This was the must-fail control and it failed as required.

### P4 — the project's own numbers are instances

Read off the real v4 modexp (N=5, a=2, wrap `msb_t0_anc`): the t0=0 slice has
structure msb⊕anc, the t0=1 slice has msb. They span d = 2, so
`E = {z : z_msb⊕z_anc = 1 and z_msb = 1} = {z_msb=1, z_anc=0}` — **exactly the
quadrant C34 found empty** (counts 7454 / 8032 / **0** / 7978), cap ¾, measured
density 0.7161. C30 is the m = 0, d = 1 case. So C34's mechanism is not special
to modexp; it is this law with d = 2.

### Scope, and the honest altitude

- **Affine invariance tested, not asserted** (P5): under a random GF(2) change
  of basis the cells become cosets of a generic subspace, the structures
  transform to L⁻¹w_u, the per-cell structures still hold pointwise, E is still
  exactly avoided, and |support| is unchanged. So the law covers arbitrary
  subspace partitions, not just coordinate-aligned ones.
- **This is a COROLLARY of textbook material**, not a deep theorem: Proposition
  29 plus the standard coset decomposition of the Walsh transform. We did not
  find the combination stated and it may be folklore. The contribution is
  recognising that this project's density caps are all instances of it, plus
  the parity condition.
- **PRIOR ART — checked as far as open access allows; residual risk low but
  not zero.** The paywalled item is Carlet–Tarannikov, *Covering sequences of
  Boolean functions and their cryptographic significance*, DCC 25:263–279
  (2002) = ref [326] of Carlet's book. Its **body was not read**. What was
  checked instead, all legitimately:
  - **Carlet's own book** — the comprehensive modern survey *by the same
    author*, which cites [326] on pp. 205, 206, 319 and reproduces its
    definitions and its Walsh characterisation (Def. 47, Prop. 60). §5.5 was
    read in full. A covering sequence is a **single global** sequence λ with
    `Σ_a λ_a D_a f(x)` constant; **partial** covering sequences (ref [231])
    relax that to two levels on a set and its complement. Ours has a
    **different structure vector per cell**, which is neither — no single λ
    reproduces it. The book does not state C40.
  - **The paper's abstract**, from Carlet's own publications page: its
    contributions are characterisations of balancedness / correlation immunity
    / resiliency, subclasses of resilient functions, and degree and
    nonlinearity bounds. None is a per-cell support-confinement law.
  - **No self-archived preprint** exists on the author's page; there is no
    legitimate free copy of the body that was found.
  - **Adjacent and worth citing anyway: Maiorana–McFarland** (book §5.1.1) —
    functions whose restrictions to each coset of a subspace are **affine**,
    whose support is confined by the image of φ. Same *flavour* (structure on
    cells constrains the support), different hypothesis (affine restriction,
    far stronger than a linear structure) and different conclusion.

  Net: C40/C41 may still be folklore, and remain an easy corollary of
  Proposition 29 plus the coset decomposition — that is the altitude to claim
  them at. But "covering sequences already contains this" now looks unlikely.

### Why this is worth having

It converts density measurement into density *prediction*. Given a construction,
find the per-cell structures, take the span, check the dependency parities, and
the cap follows without simulating anything — and it says which modifications
can possibly help: only ones that add an independent structure vector, or that
introduce an odd dependency.
