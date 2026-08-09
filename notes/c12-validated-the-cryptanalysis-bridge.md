---
code: c12-validated-the-cryptanalysis-bridge
title: "C12 VALIDATED — the cryptanalysis bridge is quantitative"
outcome: record
claims: [C12]
todo: []
---
# C12 VALIDATED — the cryptanalysis bridge is quantitative

With c_z normalised, L = 2ⁿ·max|c|, NL = 2ⁿ⁻¹(1 − max|c|), bent bound
2ⁿ⁻¹ − 2^(n/2−1). Parseval gives Σc² = 1, so a flat spectrum over S terms has
each |c| = 1/√S and NL ≈ 2ⁿ⁻¹(1 − 1/√S) — tying PPS cost S to nonlinearity.

```
                  function   n  sparsity   dens   max|c| 1/sqrt(S)  NL/NLbent
         adder b0 (affine)  12         1  0.000  1.00000   1.00000    0.00000
                  adder b1  12         4  0.001  0.50000   0.50000    0.50794
                  adder b4  12        46  0.011  0.50000   0.14744    0.50794
            modexp N=5 a=2  15     15493  0.473  0.26758   0.00803    0.73649
           modexp N=15 a=7  18    127936  0.488  0.25763   0.00280    0.74382
           modexp N=21 a=2  21   1037322  0.495  0.25588   0.00098    0.74463
            random f, n=14  14     16384  1.000  0.03357   0.00781    0.97404
            random f, n=18  18    262144  1.000  0.00965   0.00195    0.99229
```

Both endpoints land exactly: affine ⟹ sparsity 1, NL 0. Random ⟹ density 1.000,
NL/NLbent 0.97–0.99 (near-bent, as expected). **modexp sits at a stable ~0.74 of
the bent bound** — strongly nonlinear, but measurably *not* random, and the
figure is flat across 15→21 qubits.

**Caveat that must be stated:** the Parseval prediction over-estimates NL by
~35% for modexp (measured 12000 vs predicted 16252) because the spectrum is
heavy-tailed, while for random functions it is accurate to 1–3%. So S ↔ NL is a
**bound that is tight only for flat spectra**, exact at the affine and bent
endpoints, not an identity. Do not overclaim it.
