---
code: c15-confirmed-circuit-level-and-the-sharpest
title: "C15 CONFIRMED AT CIRCUIT LEVEL — and it is the sharpest result so far"
outcome: solved
claims: [C15]
todo: []
---
# C15 CONFIRMED AT CIRCUIT LEVEL — and it is the sharpest result so far

The tension in F12 (function-level dichotomy vs circuit-level density stuck at
~0.5) is resolved: **the earlier runs all used n_exp=2, where a 2-bit exponent
register leaves no room for periodicity to show.** Sweep n_exp and the dichotomy
appears immediately.

Controlled design: fix N (identical circuit, ancillas, gate structure), vary `a`
so that only r changes. `experiment_c15.py`, `experiment_c15b.py`.

**N=7, a=6 (r=2, a power of two), observable Z_x0:**

```
 n_exp  qubits        2^q   sparsity   density  vs n_exp=2    time
     2      15      32768      15549  0.474518        SAME    0.0s
     4      17     131072      15549  0.118629        SAME    0.4s
     6      19     524288      15549  0.029657        SAME    2.9s
     8      21    2097152      15549  0.007414        SAME   27.3s
    10      23    8388608      15549  0.001854        SAME  176.2s
```

**CONTROL — N=7, a=3 (r=6, has an odd factor), same circuit size:**

```
 n_exp  qubits   sparsity   density    growth
     2      15      15539  0.474213         -
     4      17      64353  0.490974     4.14x
     6      19     258691  0.493414     4.02x
     8      21    1037728  0.494827     4.01x
    10      23    4155634  0.495390     4.00x
```

Sparsity is **exactly** 15549 across a 256x growth in Hilbert-space dimension
for r=2, against a clean 4.00x-per-step growth for r=6. Same modulus, same
circuit, only the base differs.

**Independent confirmation at N=21** (`experiment_c15.py` §2): a=8 (r=2) gives
1037174 at both n_exp=2 and n_exp=4 while the dimension quadruples; a=2 (r=6)
goes 1037322 → 4186980 and a=4 (r=3) goes 512784 → 4152181.

**Mechanism.** aᵉ mod N depends only on e mod r. If r | 2ᵏ then only the low k
bits of e matter, so every additional exponent qubit adds a variable the Walsh
support cannot touch — the support is pinned to a fixed subspace. Any odd factor
in r makes the period incommensurate with the GF(2) basis and the support
spreads over everything.

**Statement of the result.**

> For reversible modular exponentiation with computational-basis observables,
> Pauli-path simulation cost is **independent of the exponent-register size when
> the order r is a power of two**, and **Θ(2^q) as soon as r has an odd factor.**

Consequences:

1. **Period-finding precision is free, or fatal, depending on r.** The exponent
   register is what sets the accuracy of the continued-fractions step; here
   enlarging it costs PPS *nothing* when r is a power of two and quadruples cost
   per two qubits otherwise.
2. **The textbook demo is the degenerate case, quantitatively.** N=15 a=7 has
   r=4. Every "we simulated Shor on N=15" result sits in the corner where PPS
   cost does not grow at all. Cryptographic N has r with odd factors
   generically, i.e. the Θ(2ⁿ) branch.
3. **Same 2-adic dependence as MPS.** Dang, Hill & Hollenberg report memory
   depending on *the factors of r*
   ([arXiv:1712.07311](https://arxiv.org/abs/1712.07311)). Two structurally
   unrelated classical methods keying on the same arithmetic property is worth
   stating as a shared fact about the algorithm, not a coincidence of either.

**Caveat.** `live_variables` reports all q bits live even in the r=2 case, so the
constancy is *not* simply "the function ignores the extra qubits". The extra
exponent qubits do influence the function (through behaviour on invalid inputs)
without enlarging the Walsh support. The clean subspace argument above explains
the valid-input structure; the full-space statement is empirical over 5 sizes
and 2 moduli, not proved.
