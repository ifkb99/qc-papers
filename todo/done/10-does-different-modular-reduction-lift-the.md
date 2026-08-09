---
id: 10
state: done
title: Does a different modular reduction lift the ceiling?
outcome: YES, ~51%
claims: [C31, C32]
---
# Does a different modular reduction lift the ceiling?

See `NOTES.md` §L2; files `experiment_reduction.py`, `experiment_reduction2.py`.

**First hypothesis refuted.** The msb↔anc CNOT pairing is *not* the cause —
conjugating the ancilla with extra CNOT/Toffoli couplings left defect 1 and the
same w in every variant. Adding XOR couplings cannot break an XOR symmetry.

**Real mechanism (C31), three ingredients:** (a) flipping the msb *is* adding
2^(m−1), which commutes with mod-2^m addition — verified exhaustively, 0/256
violations; (b) anc is coupled to msb only by XOR, so flipping both restores it;
(c) the msb is `b[n]` and the cswaps use `b[:n]`, so it never reaches x.

**Breaking it needs nonlinearity in the msb (C32):**

```
   v0 baseline                 density 0.473  defect 1
   v6 cswap(t0, msb, anc)      density 0.473  defect 1   <- linear, survives
   v4 toffoli(msb, t0, anc)    density 0.716  defect 0   <- BROKEN
   v5 toffoli(msb, t0, t1)     density 0.721  defect 0   <- BROKEN
```

All compute `a^e mod N` correctly. **~51% cost increase from a compilation
choice that changes nothing about the computed function** — the first and only
such effect in the project, everything else having been compilation-invariant.
Constant-factor only; Θ(2ⁿ) is unaffected.

Framing worth keeping: not "compile to X to make PPS cheaper", but "the standard
construction is already cheaper than it needs to be, for a reason nobody
designed".

**Left open:** broken variants sit at 0.716–0.721, not 1.000 (random reaches
1.000), so residual *non-linear* structure remains after the linear structure is
destroyed. Unidentified — a natural next thread.
