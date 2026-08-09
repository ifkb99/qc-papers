---
code: F9
title: "THE RESULT: PPS term count = Walsh sparsity, exactly"
outcome: solved
claims: []
todo: []
---
# F9 — THE RESULT: PPS term count = Walsh sparsity, exactly

For a circuit implementing a basis permutation π, `π†Z_jπ` is the diagonal
operator `(−1)^{g(y)}` with `g(y) = bit j of π(y)`. Expanding a diagonal operator
in the Pauli basis is *precisely* the Walsh–Hadamard transform of `(−1)^g`.
Therefore:

> **the number of Pauli terms PPS carries = the Walsh sparsity of g**

Verified to machine precision, 6/6 instances, both compilations, supports
identical (not merely counts):

```
 N,a  compilation   q   walsh     pps  match    maxerr
 5,2      Fourier  10     451     451    YES  6.66e-16
 5,2      Toffoli  15   15493   15493    YES  3.75e-16
 7,3      Fourier  10     461     461    YES  3.61e-16
 7,3      Toffoli  15   15539   15539    YES  4.44e-16
 4-bit adder Z_b0  10       1       1    YES  1.11e-16
 4-bit adder Z_b2  10      10      10    YES  1.11e-16
```

Consequences:

1. **Exact predictive cost model.** Walsh runs in O(2ⁿ·n) and took 0.05s where
   PPS took 49s. No extrapolation, no power-law fitting, no truncation
   heuristics — for permutation circuits the answer is computable outright.
   This sidesteps the paper's Eq. 17 machinery entirely *in this regime*.
2. **Compilation-invariant.** Kills the "compilation, not algorithm" thesis.
3. **Explains F1 correctly at last.** Adder low bit = `a0⊕b0⊕c0`, affine ⟹ Walsh
   sparsity 1. Sparsity 1 ⟺ affine is a standard theorem. The permutation
   property was never the operative fact.
4. **Bridge to cryptanalysis.** Walsh sparsity / linearity is *the* central
   quantity in linear cryptanalysis. "PPS-hard reversible circuit" ≈ "Boolean
   function resistant to linear approximation". That literature is deep and
   directly importable — likely the most valuable thread here.

**Caveat:** Walsh is itself O(2ⁿ), so this predicts cost rather than beating it.
Its value is as an *exact* cost model and an explanation, not a faster simulator.
