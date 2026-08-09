---
code: F10
title: "scope of the Walsh identity (exactly where it stops)"
outcome: record
claims: []
todo: []
---
# F10 — scope of the Walsh identity (exactly where it stops)

`experiment_scope.py` Q1, Toffoli modexp N=5, 14 qubits:

```
         observable   type   walsh     pps   nonZ  match
               Z_x0      Z    3086    3086      0    YES
               Z_x1      Z    2926    2926      0    YES
          Z_x0 Z_x1      Z    2848    2848      0    YES
     Z_x0 Z_x1 Z_x2      Z    2514    2514      0    YES
               X_x0      X       -   60358  60358    n/a (hit cap)
               Y_x0      Y       -   60358  60358    n/a (hit cap)
```

Holds for **every Z-type observable**, single- or multi-qubit (generalise via
`(-1)^{popcount(π(y) & zmask)}`). Fails completely for X/Y-type: the pullback
leaves the diagonal, every surviving term is non-Z, and it blows past the cap.

Honest scope statement: **computational-basis observables on permutation
circuits.** That is exactly what one measures in Shor (bits of the output
register), so the restriction is natural rather than convenient — but it must be
stated, not glossed.
