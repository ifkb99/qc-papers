---
code: W
title: "WEIGHT TRUNCATION (TODO step 2). Naive identity fails; knob is unusable."
outcome: record
claims: [C14, C15]
todo: [2]
---
# W — WEIGHT TRUNCATION (TODO step 2). Naive identity fails; knob is unusable.

PPS has two truncation knobs: coefficient threshold δ (all prior work here) and
Pauli **weight** (drop terms acting on > k qubits). This closes the gap.

Setup: for a permutation circuit the surviving Paulis are Z-strings, and the
Pauli weight of `Z^z` is `popcount(z)` — which *is* the Fourier **degree** of the
Walsh coefficient. So weight-k truncation is literally low-degree Fourier
truncation. The hope was that its error would be the Fourier tail, extending the
exact model to the second knob.

### W1 — the naive identity is FALSE (and the reason generalises)

**PPS truncates incrementally, at every gate.** A term discarded early never
branches, so the final result is not the truncated final operator. The Fourier
tail describes **terminal** truncation (truncate the finished operator) only.

```
modexp N=5, 14 qubits, exact <O> = -1.000000
  k   terminal <O>   term err   incremental <O>   incr err   tail mass
  2      -1.015625   1.56e-02         -0.500000   5.00e-01    0.499298
  3      -0.898438   1.02e-01         -1.000230   2.30e-04    0.494888
  4      -0.820312   1.80e-01         -1.000648   6.48e-04    0.470779
  8      -1.242188   2.42e-01         -0.644562   3.55e-01    0.188095
 14      -1.000000   0.00e+00         -1.000000   0.00e+00    0.000000
```

Both are **non-monotonic in k**, and incremental lands on either side of
terminal — at k=3 it is 400x *more* accurate, at k=8 far worse. This same
incremental-vs-terminal distinction applies to δ, and retroactively explains why
δ behaves non-monotonically too (C14): it is not just cancellation, it is that
truncation changes what subsequently branches.

### W2 — weight truncation is the WRONG KNOB for reversible arithmetic

Signed sum of coefficients per weight level (modexp N=5):

```
 weight   #terms      mass   signed sum
      1        3  0.250031    -0.500000
      2       12  0.250671    -0.515625
      3       54  0.004410    +0.117188
      4      163  0.024109    +0.078125
      5      316  0.052185    -0.140625
      7      601  0.071869    +0.312500
      8      574  0.084717    -0.515625
      9      451  0.089630    +0.296875
     13        8  0.002197    -0.031250
```

The levels carry **large, alternating** signed sums that cancel only when summed
over *all* of them. Any weight cutoff slices through that cancellation. There is
no k that both reduces cost and preserves ⟨O⟩.

### W3 — coefficient truncation is the right knob, now *exactly* explained

For the same circuit, the `|c| > 0.1` set (**4 terms** of 3086):

```
  they sum to             : -1.00000000    (exact <O> = -1.000000)
  everything else sums to : +0.00e+00      (3082 terms)
```

An **exact** split, not an approximation. δ-truncation works because the large
coefficients form a set that carries ⟨O⟩ exactly while the remainder sums to
exactly zero. Weight truncation fails because the same split is not aligned with
degree.

### W4 — the dominant coefficients are NOT all low-degree

Weights of the four: **1, 2, 8, 9**. Two of them sit high. So low-degree
truncation cannot substitute for δ — it would discard half the terms that carry
the answer. This is the concrete reason W2 holds.

Corrects an earlier misreading: the apparent "292x low-weight enrichment"
(mass 0.50 at weight ≤ 2 vs 0.026 for a random function) is **just those two
0.5-magnitude coefficients** at weights 1 and 2 (0.5² + 0.5² = 0.5). It is not a
broad low-degree structure — `k(99% mass) = 11` confirms the bulk is high-weight.

### W5 — the (z, z⊕e) pairing lead: REFUTED

Top-12 listing showed every coefficient appearing twice at equal magnitude,
differing by an exponent qubit. That looked like a hard group-theoretic
constraint and the best route to proving C15 (step 3). It is not:

```
   a   r  n_exp  |support|   all exp bits paired?
   6   2      1       3116                   True
   6   2      2      15549                  False
   3   6      1       3206                   True
   3   6      2      15539                  False
```

Exact at n_exp=1 for **both** r=2 and r=6, broken at n_exp≥2 for both
(≈0.4% of partners missing). So it is an n_exp=1 artifact, not r-dependent, and
**not the C15 mechanism. Step 3 remains open.**

Worth keeping: at n_exp=1 the pairs are almost all *opposite*-sign (3084
opposite vs 2 same), so nearly everything cancels pairwise and the two
same-sign pairs carry ⟨O⟩ — consistent with F13/W3.

**Bottom line for step 2:** the exact model does *not* extend to weight
truncation, and weight truncation should not be used on reversible arithmetic.
That is a negative result, but a sharp and actionable one, and W1 gives a
mechanism that also improves the understanding of the δ results.
