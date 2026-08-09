---
code: verified-infrastructure
title: "Verified infrastructure"
outcome: record
claims: []
todo: []
---
# Verified infrastructure

`test_core.py` — all pass. Do not trust any result if these regress.

```
[1] pauli_mult / i·σ·P / commutation vs dense, n=3, all 4096 pairs
[2] H, X, T, S, CNOT(0,1), CNOT(1,0), Toffoli decompositions; Toffoli T-count == 7
[3] Cuccaro ripple adder correct on all inputs, nbits=2,3; T-count == 14·nbits
[4] PPS(δ=0) == dense expectation, random circuits (max err ~1e-9)
[5] PPS(δ=0) == dense on a real adder circuit
```
`experiment2.py` §0 additionally verifies `CP`, `QFT` vs DFT, `QFT∘QFT⁻¹ = I`.

---
