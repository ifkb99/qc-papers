---
code: F11
title: "truncation error is NON-MONOTONIC in δ (aggressive beats mild)"
outcome: record
claims: [C5]
todo: []
---
# F11 — truncation error is NON-MONOTONIC in δ (aggressive beats mild)

`experiment_scope.py` Q2, Toffoli modexp N=5, truth = −1:

```
   delta    N_max        <O>     |err|  admissible
   0e+00    40770   -1.00000  6.66e-16         yes
   1e-04    39992   -1.00348  3.48e-03    VIOLATED
   1e-03    23482   -0.71532  2.85e-01         yes
   1e-02     2052   -1.00000  0.00e+00         yes
   3e-02      514   -1.00000  0.00e+00         yes
   1e-01       34   -1.00000  0.00e+00         yes
```

**δ=1e-1 gives the exact answer with 34 terms; δ=1e-3 is off by 0.285 with
23482 terms.** Three orders of magnitude more work, far worse answer.

Mechanism: the coefficient spectrum is a handful of large terms plus ~10⁴ small
ones (each ≈2⁻ⁿ) that **cancel among themselves**. Discard all of them and the
cancellation is preserved exactly. Discard *some* and the residue survives.
Mild truncation is the worst regime — it breaks cancellations without removing
the terms that would have completed them.

This independently reproduces the paper's own counterintuitive observation that
"reducing δ does not always improve accuracy" (Gharibyan et al., abstract) on a
completely different circuit family, and gives a concrete mechanism for it.

**C5 survives, quantified: 4/18 truncated runs were inadmissible (|⟨O⟩| > 1).**
All four were at mild δ (1e-4, 1e-3) — never at aggressive δ. The check is free
and one-sided: it can prove a run invalid, never valid.

---
