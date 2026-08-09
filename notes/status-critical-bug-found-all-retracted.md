---
code: status-critical-bug-found-all-retracted
title: "STATUS 2 — CRITICAL BUG FOUND; F4/F6/F7/F8 ALL RETRACTED"
outcome: retracted
claims: [C5]
todo: []
---
# STATUS 2 — CRITICAL BUG FOUND; F4/F6/F7/F8 ALL RETRACTED

**A bug in `pps.py` invalidated every result involving an X, Y or Z gate.**

```python
if abs(s) < 1e-12:      # WRONG: "identity up to phase"
    new = terms
```

`U†PU = cos(θ)P + sin(θ)(iσP)` for anticommuting P. At **θ = π**, `sin = 0` but
`cos = −1`, so anticommuting terms must be **negated**. The shortcut treated
every θ=π gate — i.e. every X, Y, Z — as a no-op. Fixed: the fast path now
requires `cos > 0`.

**Why the test suite missed it:** the random circuits in `test_core.py` [4] drew
from `{h, t, cnot, rx, rz}` and never emitted an `x()`. The ripple adder is built
only from CNOT/Toffoli — also no X. The bug could only fire in `ToffoliModExp`
(`_load` uncontrolled, `cc_add_mod`'s `qc.x(msb)`) and in `build()`'s `x(x[0])`.
Regression test added as `test_core.py` [4b], which fails on the old code.

**How it was caught:** Walsh (exact integer arithmetic) disagreed with PPS. The
discriminating test was running PPS in `longdouble` — the error was *identical*
to `float64` (1.000e+00), ruling out precision and proving a logic bug.
**Keep that trick: if extra precision doesn't move the error, it isn't precision.**

### Corrected results (all re-run post-fix)

```
                       exact PPS (delta=0)
 N,a  compilation   q  gates    N_max   N_fin  nonZ      <O>  true
 5,2      Fourier  10   5359     2530     451     0  -1.0000    -1
 5,2      Toffoli  15  21247    40770   15493     0  -1.0000    -1
 7,3      Fourier  10   5359     2494     461     0  -1.0000    -1
 7,3      Toffoli  15  21715    40796   15539     0  -1.0000    -1
```

Every prior claim that flipped:

- **F4 RETRACTED.** "Fourier lacks Z-closure" was pure artifact. **Both**
  compilations end with **zero** non-Z terms. Z-closure is a property of the
  *unitary being a basis permutation*, not of the gate set.
- **F6 RETRACTED.** The direction reverses: Fourier is *cheaper* (451 vs 15493),
  not 8.5x more expensive.
- **F7 SUPERSEDED by F9** (below) — right instinct, wrong framing.
- **F8 RETRACTED.** ⟨O⟩ was wrong (+1 where truth is −1), so the δ-sweep
  conclusions were meaningless.
- **C5 partially survives**: norm violation still occurs under truncation
  (Fourier N=7 δ=1e-3 gives ⟨O⟩ = −1.02086), but it is now the *only*
  truncation claim standing.

### CONFOUND that kills the A/B as designed

The two compilations use **different qubit counts** (10 vs 15) with different
ancilla layouts. They agree on the *valid subspace* but implement **different
permutations of their respective full Hilbert spaces**, so their pullbacks are
different operators over different-sized domains. 451 vs 15493 is therefore
**not** a compilation effect — it is mostly an ancilla-count effect.
**Any future A/B must match total qubit count.**

---
