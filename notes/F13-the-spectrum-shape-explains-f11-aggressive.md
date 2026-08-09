---
code: F13
title: "the spectrum shape explains F11 (aggressive truncation wins)"
outcome: record
claims: [C12]
todo: []
---
# F13 — the spectrum shape explains F11 (aggressive truncation wins)

F11 was empirical; C12's spectrum data explains it. modexp N=5, Z_x0:

```
   keep |c| >    #kept     sum kept      |err|
        0e+00    15493    -1.000000   0.00e+00
        1e-03    13005    -0.987061   1.29e-02
        1e-02     1579    -0.955322   4.47e-02
        3e-02       75    -0.771973   2.28e-01   <- worst
        1e-01        4    -1.000000   0.00e+00   <- exact, 4 terms
```

**4 of 15493 coefficients reproduce ⟨O⟩ exactly; the other 15489 sum to zero.**
The spectrum is heavy-tailed by a factor of 33 (max|c| = 0.2676 against the
flat-spectrum value 1/√S = 0.0080). Keep the spikes → exact. Keep the spikes and
*some* of the cancelling sea → residue. That is the non-monotonicity, now
predicted from spectrum shape rather than observed.

Practical rule, counterintuitive and actionable: **for these circuits truncate
aggressively, not mildly.** Mild δ is the worst regime.
