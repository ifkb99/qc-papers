---
code: X
title: "CRYPTANALYSIS IMPORT (TODO step 5). A usable bound, tight at the extremes."
outcome: record
claims: [C12]
todo: [5]
---
# X — CRYPTANALYSIS IMPORT (TODO step 5). A usable bound, tight at the extremes.

C12 gave a qualitative bridge. This makes it an inequality and tests it against
**published** nonlinearity values rather than self-measured ones.

### X1 — the bound

With normalised coefficients (Σc²=1 by Parseval) and NL = 2ⁿ⁻¹(1 − max|c|):

```
    S · max|c|² ≥ Σ c_z² = 1      ⟹      S ≥ (1 − NL/2ⁿ⁻¹)⁻²
```

Since S is exactly the Pauli support (Paper A), **any published nonlinearity
gives a lower bound on PPS cost, for every circuit computing that function,
compilation-independent, without running anything.**

```
                  function   n  S (actual)      bound  holds  tight       NL  NL/NLbent
         affine (x0^x1^x3)  10           1        1.0   True   True        0     0.0000
  inner product (bent) m=4   8         256      256.0   True   True      120     1.0000
  inner product (bent) m=6  12        4096     4096.0   True   True     2016     1.0000
           AES S-box bit 0   8         239       64.0   True  False      112     0.9333
       AES min over combos   8           -          -      -      -      112   (published: 112)
         modexp N=5 (Z_x0)  14        3086        4.0   True  False     4096     0.5039
                  random f  12        4000      264.2   True  False     1922     0.9534
```

**Independent validation:** the AES S-box nonlinearity comes out at exactly
**112**, the published value, over all 255 nonzero linear combinations. That is
a check of the whole Walsh pipeline against an external constant.

### X2 — closing the loop end-to-end

`Toffoli(x_i, y_i, out)` over i computes `out ^= <x,y>`, the inner product,
which is **bent**. Bent ⟹ flat spectrum ⟹ PPS must carry every Z-string:

```
    m  qubits  perm_pps |supp|  predicted   match
    3       7               64         64    True
    4       9              256        256    True
    5      11             1024       1024    True
```

Support is exactly 2^(2m) — the output qubit forced into every term, the bent
structure filling the rest. **No truncation can help**: every coefficient has
identical magnitude, so there is no "large" subset to keep. This is the clean
worst-case statement, and it is derived from a published property.

### X3 — HONEST LIMITATION: the bound is weak in the middle

Tight only at the two extremes, where Parseval is saturated by a single
magnitude:
- affine: bound 1, actual 1 ✓
- bent: bound 2ⁿ, actual 2ⁿ ✓
- AES: bound 64, actual **239** — off by 3.7×
- modexp: bound 4, actual **3086** — off by 770×

So as a *quantitative* predictor it is poor except at the endpoints, because NL
depends only on `max|c|` and discards the rest of the spectrum. The valuable
statements are the qualitative direction and the exact bent worst case, not the
numeric bound in between.

Sharper bounds would need more than NL — e.g. the full Walsh value/multiplicity
distribution, which *is* published for several crypto families (the AES inverse
has a known spectrum). Worth a follow-up: for functions with a fully published
spectrum, S is known **exactly**, not bounded.

### X4 — correction to a C12 claim

Earlier: "modexp sits at a stable 0.74 of the bent bound across 15→21 qubits."
That stability is across **N at fixed n_exp=2**, not across widths. At n_exp=1
the same instance gives 0.5039, because max|c| is 0.5 there versus 0.268 at
n_exp=2. **Do not state it as width-independent.**
