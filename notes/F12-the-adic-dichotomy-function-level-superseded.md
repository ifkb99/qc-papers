---
code: F12
title: "the 2-adic dichotomy at function level (superseded by C15 above)"
outcome: superseded
claims: [C15]
todo: []
---
# F12 — the 2-adic dichotomy at function level (superseded by C15 above)

Function-level probe (`/tmp/c7b.log`), g(e) = bit 0 of aᵉ mod N over e ∈ [0,2ᵗ),
which isolates the algorithm from ancilla layout:

```
    N   a    r   t        2^t   sparsity  density  r pow2
   15   7    4  12       4096          4  0.00098     yes
   15   7    4  24   16777216          4  0.00000     yes    <- CONSTANT in t
   21   2    6  24   16777216   16777216  1.00000      no    <- FULLY DENSE
   35   3   12  24   16777216   16777216  1.00000      no
  143   5   20  24   16777216   16777216  1.00000      no
  323   5  144  16      65536      64293  0.98103      no
  323   5  144  24   16777216   16777216  1.00000      no
```

The robust statement is one-sided: if r | 2ᵏ then aᵉ mod N depends only on the
low k bits of e, so g is a function of k variables and its Walsh support is ≤
2ᵏ — **constant in t**. If r has an odd factor, no universal density conclusion
follows for every selected bit, because the scalar bit function may have a
proper minimal period (for example N=13, a=4, r=6, LSB sparsity 1). N=323
(r=144=16·9) shows aliasing/period effects at small t, washed out by t=24.

Two consequences:

1. **This is the same phenomenon Dang, Hill & Hollenberg report for MPS** — that
   memory depends on *the factors of r* rather than r itself
   ([arXiv:1712.07311](https://arxiv.org/abs/1712.07311)). Two unrelated
   classical methods, same 2-adic dependence. Worth stating as a shared
   structural fact rather than a coincidence of either method.
2. **The textbook demo is the degenerate case.** N=15, a=7 has r=4, a power of
   two — sparsity 4, constant forever. Every "we simulated Shor" result on N=15
   sits in the trivially-simulable corner. For cryptographic N, r is generically
   not a power of two, so dense circuit-level behaviour is common in the measured
   family, but it is not forced for every selected output bit.
