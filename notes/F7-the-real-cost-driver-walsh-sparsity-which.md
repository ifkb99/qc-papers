---
code: F7
title: "the real cost driver is Walsh sparsity, which is compilation-invariant"
outcome: record
claims: []
todo: []
---
# F7 — the real cost driver is Walsh sparsity, which is compilation-invariant

Z-closure bounds support to 2^n Z-strings instead of 4^n Paulis — a genuine
quadratic ceiling. It does **not** make the problem cheap. The Toffoli modexp
still needs 15476 Z-strings.

Reason: for a permutation π and observable Z_j, `π† Z_j π` is the diagonal
operator `(-1)^{f(x)}` where f is the Boolean function giving bit j of the
output. Its Pauli expansion is exactly the **Walsh–Hadamard expansion of f**, so
the term count is the Walsh sparsity of f.

- ripple adder: output bit = `a0 XOR b0 XOR c0`, **linear** ⟹ 1 Walsh coefficient
  ⟹ N_final = 1. That, not "permutation", is why F1 collapsed.
- modexp: bit of `a^e mod N` is **highly nonlinear** ⟹ dense Walsh spectrum
  ⟹ 15476 terms.

**⟹ Compilation sets the ceiling (2^n vs 4^n); the algorithm's Boolean structure
sets where you sit under it.** The draft abstract's thesis ("compilation, not
algorithm") is therefore wrong as stated and needs rewriting — it is *both*, with
distinct roles. F1's original explanation was also wrong (right conclusion,
wrong mechanism).
