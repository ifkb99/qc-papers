---
code: BI
title: "BIASED INPUTS (Gangopadhyay–Kumar–Stănică–Gangopadhyay, JAMC 2023)"
outcome: record
claims: [C8, C15, C18, C21, C24, C40, C42]
todo: []
---
# BI — BIASED INPUTS (Gangopadhyay–Kumar–Stănică–Gangopadhyay, JAMC 2023)

Paper dropped in the repo root:
`Stability-of-the-Walsh-Hadamard-spectrum-of-cryptographic-Boolean-functions-with-biased-inputs.pdf`
(J. Appl. Math. Comput. 69:3337–3357, 2023). Read; two things in it are
directly ours, one is a warning, and one earlier negative result of ours now
has a principled explanation.

### The dictionary — their biased inputs ARE our biased input STATES

They study `W_f^{(1+ε)}(u) = Σ_x (1+ε)^{wt(x)} (−1)^{f(x)⊕u·x}`, the WHT when
inputs are i.i.d. Bernoulli(½+δ), and define the **stability transform**
`S_f(u) = Σ_x wt(x)(−1)^{f(x)⊕u·x}` as its first-order term in ε.

In our setting the same object appears for a completely different reason. For a
**product input state** `⊗_i (cos θ_i|0⟩ + sin θ_i|1⟩)` we have
`⟨Z^z⟩ = ∏_{i∈z} δ_i` with `δ_i = cos 2θ_i`, so

> **⟨O⟩ = Σ_z c_z ∏_{i∈z} δ_i**

— the PPS expectation is exactly a δ-biased evaluation of the Walsh spectrum,
graded by wt(z). **Verified numerically**: random product state, N=7 a=6
n_exp=3, 16 qubits — direct statevector 0.133418110567 vs the weighted Walsh
sum 0.133418110567, err 8.2e-14. So their whole framework is "PPS with a biased
product input", and their stability spectrum is the sensitivity of a PPS
estimate to input-state bias. Their δ = 0 is our maximally-mixed direction and
δ = 1 is the computational basis.

### C42 — for the REAL Shor input, most of the support is dead weight

Shor's initial state puts the exponent register in |+⟩, i.e. **δ = 0 on those
qubits**, so *every* z with any exponent-register support contributes **exactly
zero** to ⟨O⟩. Combined with C24 (`z_I ∈ {0, 1_I}` on the tail, exp bits below α
free) this predicts a useful fraction of exactly **2^−(α+1)**. Derived, then
measured:

```
  N=7 a=6  α=1   n_exp=2,3,4:  |S|=15549   dead 75.03%   useful 3883
  N=5 a=4  α=1   n_exp=2,3,4:  |S|=15509   dead 74.99%   useful 3879
  N=5 a=2  α=2   n_exp=2:      |S|=15493   dead 75.00%   useful 3873
  N=5 a=2  α=2   n_exp=3,4:    |S|=32143   dead 87.51%   useful 4014
```

75% at α=1 and 87.5% at α=2 = 1 − 2^−(α+1), and the α=2 row **jumps exactly at
n_exp = 3 = α+1**, reproducing C21's onset from an independent direction. The
useful count is itself constant in n_exp, so C15 holds for useful work as well
as for total cost.

**Honest limit — this is not a free speedup.** PPS propagates the observable
backwards and only meets the input state at the end, so a term is only known to
be dead once propagation is finished. It does **not** reduce peak memory (C18)
as stated. What it does is separate *cost* from *useful work*: at α=2 seven of
every eight Pauli terms carried are irrelevant to the answer. Whether the dead
set can be predicted early enough to prune is open and worth asking.

### Reference for the p-biased machinery

O'Donnell, *Analysis of Boolean Functions*, is free on arXiv
([arXiv:2105.10386](https://arxiv.org/abs/2105.10386)) — **Chapter 8** is the
p-biased Fourier expansion, which is the proper setting for the dictionary
above and for anything further in this direction. Use it rather than
re-deriving.

### Why §W's negative result was structural, not bad luck

§W concluded that weight truncation is the wrong knob while δ (magnitude)
truncation is exactly right. Their **Theorem 10** gives the reason: the
S-spectrum is covariant only under **weight-preserving orthogonal** A, and they
exhibit a counterexample showing it is **not invariant under extended affine
equivalence**. Hamming weight is not an affine invariant; the Walsh support and
|c_z| are. So a weight-graded truncation is basis-dependent by construction,
whereas everything else in this project (C8, C40, the PPS cost model) is affine
covariant. §W's failure was forced.

Their grading is on wt(x) and ours on wt(z), but the transform swaps the two
sides and the argument is the same either way. Ties directly to §GF: affine
structure is what the Walsh basis respects, and weight is not part of it.

### What it does NOT give us

No cost model, no permutation circuits, no connection to the 2-adic story.
Their results are about bent and symmetric functions (Theorem 4's bound
|S_f(u)| ≤ (n/2)2^{n/2}, Maiorana–McFarland stability, symmetric
classifications) — none of which our pullbacks are. Cite it for the biased-input
dictionary and for Theorem 10; do not lean on it for anything else.
