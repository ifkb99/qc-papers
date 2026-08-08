# Pauli-path simulation of reversible quantum arithmetic

Research code for two related results about Pauli Path Simulation (PPS, also
called sparse Pauli dynamics) applied to reversible arithmetic and Shor's
algorithm.

**Paper A — `ABSTRACT.md`.** For any circuit implementing a permutation of the
computational basis, the Pauli support carried by PPS with a computational-basis
observable is *exactly* the Walsh–Hadamard spectrum of the corresponding output
bit function. Cost is therefore computable in closed form rather than
extrapolated. Verified to machine precision, supports identical, on 6/6
instances across two independent compilations.

**Paper B — `ABSTRACT_SHOR_2ADIC.md`.** Applying that identity to modular
exponentiation: PPS cost is governed by the 2-adic structure of the
multiplicative order r. Writing r = β·2^α with β odd, cost is independent of
exponent-register width when β = 1 and Θ(2^n) otherwise. Verified in a
controlled design (fix modulus, vary base) at three moduli, for both final
support and peak memory, and **proved** via a parity-reduction argument.

## Start here

**`HANDOFF.md`** — if you are picking this up fresh, read that first. It states
where things stand, the single open thread, how to get running, and what not to
redo.

`METHOD.md` records how this work was actually conducted — the research loop,
the failure modes hit, and the techniques that earned their keep. Worth reading
before extending anything, since several of the lessons cost real time.

`CLAIMS.md` is the single source of truth for claim statuses (C-numbers,
F-numbers, and the retracted/dead ones). Both abstracts point at it rather than
carrying ledgers of their own.

`NOTES.md` is the working record — conventions, verified results, retracted
claims, traps, and an honesty log of everything believed and then killed. Read
it before trusting any number. Several early findings (F1–F8) are retracted;
the status blocks say which.

## Environment

Python 3.14 via `uv`. Run everything from this directory:

```bash
uv run python test_core.py          # correctness gate -- run first, always
```

`source .env` only when Julia / PauliPropagation.jl is needed — its `LD_PRELOAD`
of Julia's libstdc++ segfaults numpy `longdouble` otherwise.

## Test suites (all must pass)

```bash
uv run python test_core.py           # Pauli algebra, gate decompositions, PPS vs dense
uv run python test_modexp.py         # Beauregard modexp, layers A-G
uv run python test_toffoli_arith.py  # Toffoli modexp, layers A-G + cross-check
uv run python test_walsh.py          # logical trace, Walsh machinery
uv run python test_perm_pps.py       # permutation-native PPS vs Walsh and vs PPS
uv run python test_lab.py            # lab/ engine, pinned to logged numbers
uv run python test_claims.py         # CLAIMS.md headline rows, re-verified
```

Sanity anchor: `ModExp(15, 7, n_exp=4).build_shor()` must give exponent-register
peaks at y = 0, 4, 8, 12 with p = 0.25 each (r = 4).

## Layout

| file | role |
|---|---|
| `pauli.py` | symplectic Pauli algebra, the single rotation/conjugation rule |
| `circuits.py` | `Circuit` as (σ,θ) list + logical trace; gate decompositions; Cuccaro adder; QFT |
| `pps.py` | rotation-level Pauli propagation with δ-truncation |
| `perm_pps.py` | permutation-native PPS — stays Z-type throughout, 2× lower peak, ~10× faster |
| `walsh.py` | permutation extraction, FWHT, pullback coefficients, sparsity, affineness |
| `statevec.py` | O(2ⁿ)-per-gate state-vector simulation (batched) |
| `modexp.py` | Beauregard (Fourier-arithmetic) modular exponentiation |
| `toffoli_arith.py` | Toffoli-compiled modular exponentiation |
| `lab/` | experiment engine: protocol harness, cached measurement, GF(2)/structure analysis, modexp variants, null models |
| `experiments/` | the experiment scripts — lab-notebook records, filenames unchanged; run as `uv run python -m experiments.<name>`; start new ones from `TEMPLATE.py` |
| `archive/` | fully retracted scripts, kept for the record |
| `out/` | gitignored: logs and the measurement cache |

Notable experiments: `experiment_c8` (the Walsh identity, Paper A core),
`experiment_c7` (scaling to 24 qubits, density → ½), `experiment_c15{,b}` (the
2-adic result, Paper B core), `experiment_resid{1,2}` + `experiment_affstruct`
(the conditional-structure resolution, C33–C35). `experiment.py` and
`experiment2.py` are partially retracted — see their headers.

## Scope and limits

- The identity covers **computational-basis (Z-type) observables on permutation
  circuits**. It fails for X/Y-type observables, where the pullback is not
  diagonal.
- It gives an *exact cost model*, not a faster simulator: the Walsh transform is
  itself O(2ⁿ).
- Nothing here bears on the classical hardness of factoring. Efficient
  simulation of Shor's algorithm on general inputs would be a classical
  factoring algorithm; these are diagnostic results about where PPS breaks.
