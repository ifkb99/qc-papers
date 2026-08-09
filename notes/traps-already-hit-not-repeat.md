---
code: traps-already-hit-not-repeat
title: "Traps already hit (do not repeat)"
outcome: superseded
claims: []
todo: []
---
# Traps already hit (do not repeat)

1. **Floating-point residue counted as real Pauli terms.** Exact cancellations
   leave ~1e-16 junk. `pps.py` now floors at `abs(v) > 1e-13` even when δ=0.
   Before the fix, `experiment.py` §3 was measuring pure fp noise (coefficients
   reported at 2^-53, 2^-105 — garbage). **`experiment.py` §3 output is stale/
   invalid**; §1 and §2 are still fine.
2. **QFT bit-reversal acts on the *input* index (columns), not rows.**
   `U ≈ phase * F[:, bitrev]`. Cost an hour. `QFT∘QFT⁻¹=I` passing does *not*
   validate the convention.
3. **Degenerate observables.** In the H→adder→QFT toy, every a-register
   observable has `<O> = 0` by symmetry (an adder has no period structure). The
   only non-zero ones found were ancillas that return to |0>. **Cannot do a real
   accuracy-vs-δ study without genuine modular exponentiation.** This is the
   main reason the modexp generator was built — and building it immediately
   killed F2, so the instinct was right.
4. **QFT swap network goes BEFORE the body, not after.** The H/CP body puts a
   bit reversal on the *input* index (`U = phase * F[:, bitrev]`). Forward QFT is
   therefore `Body . Swaps` — swaps applied to the state first. Getting this
   backwards silently produces a unitary that still satisfies `QFT∘QFT⁻¹ = I`
   **and still passes a casual eyeball test**, but breaks Fourier arithmetic
   (bit reversal does not commute with addition) and *silently changed all of
   F2/F3*. Costliest bug of the session. `circuits.qft(..., swaps=True)` is now
   verified against the plain DFT with sign +1.
5. **(SUPERSEDED — environment changed.)** Run everything with `uv run python`
   from `research/` (Python 3.14 + `.venv`, numpy 2.4.6, scipy, stim, qiskit,
   quimb, numba, juliacall). Only `source .env` when Julia/PauliPropagation.jl is
   needed — its `LD_PRELOAD` of Julia's libstdc++ **segfaults numpy longdouble**
   (exit 139), which cost a debugging cycle. The old advice to use
   `/usr/bin/python3` predates the venv and no longer applies.
6. **`pkill` before a heredoc in the same command kills the write.** Two scripts
   vanished this way. Write the file first, kill second — or use the Write tool.
7. **CPython 3.14 crashes on long perm_pps runs.**
   `Fatal Python error: _TAIL_CALL_CACHE: Executing a cache.` — a bug in 3.14's
   new tail-calling interpreter, not in this code. It killed two background jobs
   mid-run and looked like a hang or a stray `pkill`. If a long run dies with no
   traceback, suspect this first. Workaround: rerun (it is intermittent), or use
   a 3.12/3.13 interpreter for long jobs.

---
