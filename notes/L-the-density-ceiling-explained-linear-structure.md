---
code: L
title: "THE DENSITY-½ CEILING EXPLAINED: a linear structure from the reduction ancilla"
outcome: solved
claims: [C7, C12]
todo: []
---
# L — THE DENSITY-½ CEILING EXPLAINED: a linear structure from the reduction ancilla

Two observations had gone unexplained across the whole project:

- **C7:** circuit density converges to 0.498 — approaching ½ **from below** and
  never crossing it — where random Boolean functions give 1.000.
- **Step 6:** the real modexp table at r=6 gives density exactly 0.500 on even t
  where a random table of the same period gives 1.000.

Approaching ½ from below without crossing is the signature of a **linear
constraint on the support**. Confirmed:

```
                                |supp|  density   GF(2) rank  defect
  modexp N=5 a=2 n_exp=1          3086  0.188354      13/14        1
  modexp N=7 a=6 n_exp=2         15549  0.474518      14/15        1
  modexp N=15 a=7 n_exp=1        31176  0.237854      16/17        1
  random f, n=14                 16384  1.000000      14/14        0
```

The support has GF(2) rank **n−1**, so it lies in a hyperplane and the density
is capped at exactly ½. Dually, that means a **linear structure**: a nonzero w
with `g(y ⊕ w) = g(y)` for all y — verified pointwise.

### What w is

**w = b_msb ⊕ anc** — the accumulator's sign bit XOR the modular-reduction
comparison ancilla. `{b3, anc12}` at N=5,7; `{b4, anc15}` at N=15.

The pairing is visible directly in the source. `cc_add_mod` does

```python
qc.cnot(msb, self.anc)                        # set the negative flag
...
qc.x(msb); qc.cnot(msb, self.anc); qc.x(msb)  # uncompute it
```

so flipping msb and anc *together* is a symmetry of the reduction step.

### Scope: circuit property, not observable, and shared by both compilations

```
  Toffoli, N=7:  Z_x0 Z_x1 Z_x2 Z_b0  all give defect 1, w = b3,anc12
  Fourier, N=5:  defect 1, w = b3,anc7      (its own sign bit + ancilla)
  Fourier, N=7:  defect 1, w = b3,anc7
```

**Observable-independent** and present in **both** compilations, which share the
add / subtract-N / conditional-restore reduction. So it follows from the
reduction *discipline*, not from one implementation — but it is not a property
of modular exponentiation as an algorithm. A reduction that computes its
comparison flag differently might not have it.

### Three consequences

1. **C7's ½ is explained exactly** — the support is filling a hyperplane whose
   density is ½ by construction.
2. **C7's Θ(2ⁿ) conclusion is unaffected** (½·2ⁿ is still Θ(2ⁿ)), but **the
   constant ½ must be reported as construction-dependent**, not algorithmic.
   This is a correction to how C7 is currently phrased.
3. **Feeds back into C12.** Bent functions have **no** linear structures at all,
   so possessing one bounds modexp away from the bent bound automatically. That
   partly explains why modexp sits at 0.74 of the bent bound where random
   functions reach 0.97–0.99 — it is not merely "less random", it carries a
   specific, standard cryptanalytic weakness.

The ripple adder shows the same phenomenon more strongly: kernel dimension 5,
so density ≤ 2⁻⁵ = 0.031 (measured 0.0098).
