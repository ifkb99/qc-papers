---
code: F6
title: "C3 RESOLVED: Z-closure is exact, but the benefit is ~8x, not 10^4"
outcome: solved
claims: [C3]
todo: []
---
# F6 — C3 RESOLVED: Z-closure is exact, but the benefit is ~8x, not 10^4

`toffoli_arith.py` built and verified (`test_toffoli_arith.py` A–G pass,
including a cross-check that both compilations give identical `a^e mod N`).
The matched A/B, N=5 a=2 n_exp=2, observable Z_x0, **δ=0 (exact)**:

```
compilation  qubits  rots  nonCliff   N_max  N_final  Z-type  non-Z  max|non-Z|
   Fourier       10  5359      2262  131064   131064     508 130556    4.71e-02
   Toffoli       15 21247      4074   40774    15476   15476      0    0
```

**Z-closure is real and exact.** The Toffoli compilation ends with *zero*
non-Z-type terms — not "small", exactly zero. The Fourier compilation ends with
130556 non-Z terms carrying real weight (largest 4.7e-2, so not fp noise).

**Why Fourier lacks closure even though it computes a permutation:** Beauregard's
circuit is a permutation only *on the valid subspace* (b=0, anc=0, x<N). As a
full 2^n unitary it is not a permutation matrix — the phase rotations only
conspire on that subspace. Toffoli/CNOT/X is a permutation matrix on the whole
Hilbert space, unconditionally.

**But the payoff is modest:** 15476 vs 131064 final terms (~8.5x), N_max 40774 vs
131064 (~3.2x). Not the "four orders of magnitude" the draft abstract claimed.
