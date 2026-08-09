---
code: conventions-get-these-wrong-and-everything
title: "Conventions (get these wrong and everything silently breaks)"
outcome: record
claims: []
todo: []
---
# Conventions (get these wrong and everything silently breaks)

**Pauli representation** (`pauli.py`): string = `(x, z)` bitmask pair,

```
P(x,z) = i^{popcount(x & z)} * X^x Z^z
```

The `i^{|x&z|}` is what makes P Hermitian (turns XZ into Y). Consequence:
**all coefficients in the Heisenberg expansion stay real.** Qubit 0 is the
least significant bit; `to_matrix` builds with `np.kron(m, M)` in that order.

`pauli_mult(p,q) -> (r, k)` with `P(p)P(q) == i^k P(r)`, where
`k = |a&b| + |x&z| - |rx&rz| + 2*|b&x|  (mod 4)`. Verified exhaustively n=3.

**Everything is a Pauli rotation.** `U = exp(-i θ σ / 2)`, one conjugation rule:

```
U† P U = P                            if [P,σ]=0
       = cos(θ) P + sin(θ) (i σ P)     if {P,σ}=0
```

Clifford = θ=±π/2 ⟹ cos=0 ⟹ maps to a *single* Pauli, no branching.
All branching comes from non-Clifford angles. This is the whole cost mechanism.

**Gate decompositions** (`circuits.py`, all verified up to global phase):
- `H = Rz(π/2) Rx(π/2) Rz(π/2)`
- `X = Rx(π)`, `T = Rz(π/4)`, `S = Rz(π/2)`
- `CNOT(c,t) = Rzx(-π/2) · Rz_c(π/2) · Rx_t(π/2)`, σ for the ZX term is `(1<<t, 1<<c)`
- `Toffoli` = standard 15-gate Clifford+T, **exactly 7 T gates** (asserted in tests)
- `CP(θ) = Rz_a(θ/2) Rz_b(θ/2) Rzz(-θ/2)`

**Ordering.** `Circuit.gates[0]` is applied to the state *first*. PPS therefore
iterates `reversed(circuit.gates)` (Heisenberg: innermost conjugation is the
last gate). Getting this backwards gives plausible-looking wrong answers.

**Expectation readout.** `<0|P(x,z)|0>` = 1 if `x == 0`, else 0. So
`<O> = sum of coeffs over Z-type terms`.

---
