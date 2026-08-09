---
code: T
title: "DOES C15 SURVIVE TRUNCATION? (TODO step 7). Peak cost: yes. Accuracy: to a point."
outcome: record
claims: [C5, C15, C24]
todo: [7]
---
# T — DOES C15 SURVIVE TRUNCATION? (TODO step 7). Peak cost: yes. Accuracy: to a point.

Every C15 figure was δ=0, and nobody runs PPS at δ=0 — so the practical claim
("period-finding precision is free when r is a power of two") was unsupported in
the regime where it would actually be used.

**Predictions were written before measuring** (`experiment_c15_trunc.py` header):

- **P1** Terminal truncation *must* preserve constancy at every δ, since C24/M2
  give identical magnitude multisets across n_exp and a magnitude threshold on
  identical multisets keeps identical counts. Zero exceptions expected; failure
  would mean C24 or M2 is wrong.
- **P2** Incremental truncation (what PPS does) *may* break it — more n_exp means
  more gates, hence more truncation events, and W1 already showed incremental and
  terminal are different operations.
- **P3** Control: β>1 must fail at every δ.

### T1 — P1 confirmed exactly

```
 N=7 a=6 (r=2, beta=1)          delta:  0    1e-6   1e-4   1e-3   1e-2  3e-2  1e-1
   n_exp=3                            15549  15549  15549  13194  1712    75     4
   n_exp=4                            15549  15549  15549  13194  1712    75     4
   n_exp=5                            15549  15549  15549  13194  1712    75     4
 N=5 a=2 (r=4, beta=1)
   n_exp=3,4,5                        32143  32143  32143  26257  1704    36     8   (identical)
 N=7 a=3 (r=6, beta=3)  CONTROL
   n_exp=3                            30712  30712  30712  23849  1442    40     4
   n_exp=4                            64353  64353  62982  44480  1225    51     4
   n_exp=5                           129012 129012 123137  69849   994    43     4
```

Identical across every width at every δ for β=1; the control diverges. **The
invariance is not an artifact of working exactly — it holds at every truncation
level**, because the underlying magnitudes are identical.

### T2 — P2: the practically relevant cost survives; accuracy has a limit

Incremental (real PPS), N=7 a=6:

```
  n_exp   delta     N_max   N_final         <O>
      3   0e+00     24369     15549   -1.000000
      3   1e-02      1028         8   -1.000000
      3   1e-01        34         8   -1.000000
      4   1e-02      1028        16   -1.000000
      4   1e-01        34         0   +0.000000     <- accuracy lost
      5   1e-02      1028        32   -1.000000
      5   1e-01        34         0   +0.000000

  constancy across n_exp:  delta=0     N_max SAME  N_final SAME  <O> SAME
                           delta=1e-4  N_max SAME  N_final SAME  <O> SAME
                           delta=1e-2  N_max SAME  N_final differs  <O> SAME
                           delta=1e-1  N_max SAME  N_final differs  <O> differs
```

**N_max — peak memory, the quantity that actually bounds cost — is SAME at every
δ tested** (24369, 24369, 1028, 34). That is the practically relevant statement
and it survives intact.

What does *not* survive: `N_final` drifts at aggressive δ (8 → 16 → 32), and at
δ=1e-1 the estimate collapses to ⟨O⟩ = 0 for the larger widths. Mechanism is
W1's: more exponent qubits means more gates, hence more incremental truncation
events, so terms that a terminal threshold would have kept are destroyed en
route. The larger circuit is *more* fragile at the same δ despite having an
identical exact spectrum.

### T3 — control confirmed under incremental truncation

The β>1 control was cut short by the CPython 3.14 crash (trap 7) but got far
enough to settle the question:

```
  N=7 a=3 (r=6, beta=3) CONTROL
   n_exp   delta     N_max   N_final         <O>
       3   0e+00     48855     30712   -1.000000
       4   0e+00     98018     64353   -1.000000
       5   0e+00    196060    129012   -1.000000     <- doubling per step
       3   1e-04     48036     30779   -1.001772     <- norm violation
       4   1e-04     94827     62819   -0.997549
```

N_max grows ~2× per added qubit, against SAME at every δ for β=1. And δ=1e-4
produces ⟨O⟩ = −1.001772, i.e. |⟨O⟩| > 1 — another instance of C5, now on the
control. N=5 a=2 (β=1) likewise gave N_max 48972 at both n_exp=3 and 4.

### Verdict for Paper B

The practical claim stands, with a stated boundary:

> Peak Pauli-path cost is independent of exponent-register width for β=1, at
> every truncation level tested. Accuracy is likewise independent up to
> moderate δ, but at aggressive δ the wider circuit degrades first — not because
> its exact spectrum differs (it is identical) but because incremental
> truncation has more gates to act on.

So "period-finding precision is free" is true for *memory* and true for
*accuracy up to moderate δ*, and must not be stated unqualified.
