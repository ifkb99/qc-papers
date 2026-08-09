---
code: findings-narrowed-retracted-superseded-see
title: "Findings (F1 narrowed, F2 retracted, F3 superseded by F5 — see above)"
outcome: superseded
claims: []
todo: []
---
# Findings (F1 narrowed, F2 retracted, F3 superseded by F5 — see above)

### F1 — Gate-level permutation circuits are FREE for PPS with Z-type observables

**Scope: Toffoli/CNOT/X circuits only. Does NOT extend to Fourier arithmetic.**

3-bit Cuccaro adder, 42 T gates, observable `Z_b0`:

```
peak Pauli terms 128  →  final terms 1  →  <O> exact at every δ tested (incl 1e-2)
```

**Mechanism:** a classical reversible circuit is a permutation matrix.
Conjugating a *diagonal* operator by a permutation stays diagonal. Z-type Pauli
strings span the diagonals ⟹ Z-strings map to Z-strings, and every intermediate
branch into an X/Y string **cancels exactly**. Branching is real but transient.

**⟹ Toffoli-density does NOT imply PPS-hardness.** This is elementary once seen
and is almost certainly known. Treat as corrected foundation, not a result.

### F2 — The QFT is the bottleneck, not the arithmetic

Adding an inverse QFT on the a-register drives the observable off-diagonal, and
*then* the T gates branch for real:

```
                 Z-type weight fraction:  start   min    end    mean
arithmetic only                           1.000  0.000  1.000  0.484
+ inverse QFT                             1.000  0.000  0.000  0.004
```

Scaling (δ=0, exact PPS; "sandwich" = H layer → adder → inverse QFT):

```
 nbits qubits Tgates | arith N_max fin | sand N_max    fin  ratio
     2      6     31 |          32   1 |        138     52    4.3
     3      8     51 |         128   1 |       1158    384    9.0
     4     10     74 |         512   1 |      10306   2916   20.1
     5     12    100 |        2048   1 |      94006  22592   45.9
     6     14    129 |        8192   1 |     872306      -  106.5 (CAP)

sandwich: log2(N_max) ~ 1.57 bits per qubit   (1.0 = 2^n, 2.0 = 4^n worst case)
```

arith column *always* ends at 1 term (F1). Sandwich never collapses, and the
QFT's multiplier is itself growing exponentially.

**⟹ For Shor-like circuits the PPS bottleneck is the QFT, not the modular
exponentiation — the opposite of where the gate count sits.** This is the most
interesting thing found so far.

### F3 — δ is not a working dial for this family

```
nbits=5 (12 qubits):  δ=1e-2  N_max= 2974   slope 2.09
                      δ=1e-3  N_max=41398   slope 1.14
                      δ=1e-4  N_max=91446   slope 0.34
                      δ=1e-6  N_max=94006   slope 0.01
                      δ=0     N_max=94006   (exact)
```

Below ~1e-4 δ does nothing: coefficient spectrum has a **floor**, not a
power-law tail, so truncation has nothing left to discard. Accuracy fails as a
**cliff**, not gracefully — on a non-degenerate observable: exact through
δ=3e-2, then error jumps to 1.0 at δ=1e-1.

Matters because the paper's practical contribution (extrapolate N_max from cheap
test runs, their Eq. 17 `N_max ~ δ^-m`) needs a power law to extrapolate along.

---
