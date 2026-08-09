---
code: status-retracted-narrowed-after-building-real
title: "STATUS 1: F2 RETRACTED, F1 NARROWED (after building real modexp)"
outcome: retracted
claims: []
todo: []
---
# STATUS 1: F2 RETRACTED, F1 NARROWED (after building real modexp)

Building the real Beauregard modexp (`modexp.py`) invalidated two of the three
findings below. **Read this block before trusting F1–F3.**

- **F2 is RETRACTED.** The "QFT is the bottleneck" numbers came from a QFT with
  a bit-reversal bug (trap 4) *and* a toy where the adder acted on the b-register
  while the QFT acted on the a-register — so the observable barely met the
  arithmetic. With the QFT fixed, the toy sandwich collapses to a constant
  `N_max = 8` at every size. The old scaling table (138 → 872306, 1.57 bits/qubit)
  is an artifact. Do not cite it.
- **F1 is NARROWED.** It holds for *gate-level* permutation circuits only.
- **F3 survives in corrected form**, see F5.

Real-circuit numbers (N=5, a=2, n_exp=3, 11 qubits; modexp = 8038 rotations /
3339 non-Clifford; δ=1e-3):

```
   observable   circuit     N_max         <O>
       Z_exp0    modexp         1   +1.000000     <- trivial: Z on a control qubit
       Z_exp0      shor     46756   +0.072884
         Z_x0    modexp     13313   +1.022676     <- |<O>| > 1 : INVALID estimate
         Z_x0      shor     13313   +0.222552
         Z_x1    modexp     17555   +1.021963     <- also invalid
         Z_x1      shor     17555   +0.162320
         Z_b0    modexp     17126   +0.455090
         Z_b0      shor     17126   +0.025288
```

### F4 — PPS cost depends on the *gate-level* arithmetic, not the logical map

Toffoli-based adder (`circuits.ripple_adder`) is a permutation **gate by gate**,
so Z-strings stay Z-type and everything cancels (F1). Beauregard modexp computes
the *same kind of logical map* but in Fourier space, where the individual gates
are CP/Rz rotations and are **not** permutations. Result: Z_x0 through the modexp
gives `N_max = 13313`, not 1.

**⟹ Same logical arithmetic, opposite PPS profile, depending only on whether you
compile to Toffoli or to Fourier rotations.** This was thread 3 (a curiosity);
it is now the main result and the main live question.

Also note: adding the inverse QFT changes `N_max` **not at all** for x/b-register
observables (13313 both, 17555 both, 17126 both) — it only shifts `<O>`. The
modexp dominates. This is the direct refutation of F2.

### F5 — Truncation produces provably invalid estimates (corrected F3)

`<Z_x0> = +1.022676` at δ=1e-3. For any Pauli observable |<O>| ≤ 1, so this is
not a slightly-wrong answer, it is an *impossible* one. Truncation error here is
not small and not bounded. This is much harder evidence than the old F3 slope
argument, and it does not depend on the toy circuits.

**Cheap validity check to keep using: assert |<O>| ≤ 1.** Free, and it caught
this immediately.

---
