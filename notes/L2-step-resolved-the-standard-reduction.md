---
code: L2
title: "STEP 10 RESOLVED: the standard reduction is *accidentally* PPS-friendly"
outcome: solved
claims: [C7, C15]
todo: [11]
---
# L2 — STEP 10 RESOLVED: the standard reduction is *accidentally* PPS-friendly

First hypothesis (that the msb↔anc CNOT pairing causes it) was **refuted**:
conjugating the ancilla with extra couplings — CNOT from t0, Toffoli from t0&t1,
CNOT from b0 — left the structure completely unchanged, defect 1 with the *same*
w in every case. Adding XOR couplings cannot break an XOR symmetry.

### Why the structure is robust — three ingredients

**(a) An msb flip commutes with the arithmetic.** Flipping the msb of an m-bit
register *is* adding 2^(m−1), and mod-2^m addition is commutative:
`(b ⊕ 2^(m−1)) + c == (b + c) ⊕ 2^(m−1)`. Verified exhaustively — **0 violations
over all 256 (b,c) pairs at m=4**. So the flip propagates through every
`_add_const` untouched.

**(b) anc is coupled to msb only by XOR** (`cnot(msb,anc)` and the x/cnot/x
uncomputation), so flipping both together restores anc exactly.

**(c) The msb never reaches the observed register.** msb is `b[m−1] = b[n]`, and
the cswaps in `u_a` use `zip(x, b[:n])` — they **exclude** it. So the x register
cannot see the flip.

Together these force `g(y ⊕ w) = g(y)` with w = msb ⊕ anc, hence density ≤ ½.

### Breaking it requires NONLINEARITY in the msb

```
   mode                                 density  defect          w
   v0  baseline                        0.472809       1   b3,anc12
   v6  cswap(t0, msb, anc)             0.472809       1   b3,anc12   <- linear, survives
   v4  toffoli(msb, t0, anc)           0.716064       0   FULL RANK  <- BROKEN
   v5  toffoli(msb, t0, t1)            0.720764       0   FULL RANK  <- BROKEN
```

All four compute `a^e mod N` correctly — the wraps are controlled on scratch
qubits that are 0 on the valid subspace. v6 fails to break the symmetry because
a swap is *linear* over GF(2); only the Toffoli variants, which use the msb
nonlinearly, destroy it.

**Density jumps 0.473 → 0.716, a ~51% cost increase**, from a modification that
changes nothing about what the circuit computes.

### The finding

> **The standard Beauregard/Vedral modular reduction is accidentally friendly to
> Pauli-path simulation.** Its factor-of-two saving is not designed in — it is a
> byproduct of the msb being excluded from the swap network while commuting with
> the adder. A construction that touches the msb nonlinearly forfeits it and
> costs ~51% more.

This is the **first genuine compilation-dependent cost effect** in the project.
Everything else — the Walsh identity, C15, the 2-adic dichotomy — turned out
compilation-invariant, and the one earlier claim to the contrary (F4/F6) was a
bug. Note it is a *constant-factor* effect: 0.716·2ⁿ is still Θ(2ⁿ), so C7's
asymptotic conclusion is untouched.

The framing is also the reverse of the usual one. Not "compile to X to make PPS
cheaper", but "the standard compilation is already cheaper than it needs to be,
for a reason nobody designed."

**Still open:** the broken variants sit at 0.716–0.721, not 1.000, while random
functions reach 1.000. So residual non-linear structure remains after the linear
structure is destroyed. Unidentified.

**→ RESOLVED (TODO 11): see §RS.** It is a conditional linear structure with an
exact ¾ density cap; the quadrant {z_msb=1, z_anc=0} is exactly empty.
