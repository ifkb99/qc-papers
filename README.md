# Pauli-path simulation of reversible quantum arithmetic

Research code for two related results about Pauli Path Simulation (PPS, also
called sparse Pauli dynamics) applied to reversible arithmetic and Shor's
algorithm.

The research goal behind them is broader: search existing mathematics for
anything that removes an exponential from these simulations. `METHOD.md` ("The
goal and the barrier check") gives the rules, and `transfers/INDEX.md` lists
what has been imported or ruled out.

**Paper A — `PAPER_A.md`** (full draft; `ABSTRACT.md` is now the abstract workshop and prior-art dossier). For any circuit implementing a permutation of the
computational basis, the Pauli support carried by PPS with a computational-basis
observable is *exactly* the Walsh–Hadamard spectrum of the corresponding output
bit function. Final full-operator support is therefore characterized exactly
by a Walsh transform, which remains exponential to evaluate in general.
The identity is verified to machine precision on the reported
instances; the supplied Fourier and Toffoli examples agree on the valid
subspace but do not implement the same full-space permutation.

**Paper B — `PAPER_B.md`** (full draft; `ABSTRACT_SHOR_2ADIC.md` is now the abstract workshop). Applying that identity to modular
exponentiation: PPS cost is governed by the 2-adic structure of the
multiplicative order r. Writing r = β·2^α with β odd, the studied circuit has
width-independent cost once the β=1 identity tail is present; Θ(2^n)-type
growth for β>1 is measured in the reported circuit families, not proved for
every selected bit. Verified in a controlled design at three moduli, for both
final support and peak retained terms in uncontracted PPS, with the β=1 mechanism proved by parity
reduction.

**New state-aware extension (2026-09-09).** For the pre-inverse-QFT work
observable, finished exponent controls can be contracted exactly. The reduced
support is bounded independently of exponent width in both order branches;
the β=1 tail has an additional idempotent shortcut. See
`notes/CT-finished-control-contraction.md` and claims C45/C46 for proofs,
measurements and limitations. This does not simulate Shor output sampling.

```python
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm

me = ToffoliModExp(N=7, a=3, n_exp=6)
result = propagate_perm(me.build(), 1 << me.x[0], trace_plus=me.exp)
```

Omit `trace_plus` to retain the original full-operator behavior.

**Follow-up investigations (2026-09-09/10).** The four directions were pursued
in order; each has a saved experiment with predictions and negative controls:

| Direction | Record | Experiment module |
|---|---|---|
| Balanced controls and exact tail error | `notes/BC-balanced-control-channels.md` | `experiments.experiment_balanced_controls` |
| Tensor rank versus Walsh support | `notes/TR-walsh-tensor-memory.md` | `experiments.experiment_tensor_memory` |
| Conditional order-finding output | `notes/CF-conditional-order-finding.md` | `experiments.experiment_conditional_order_finding` |
| Scratch-extension equivalence | `notes/SE-scratch-extension-equivalence.md` | `experiments.experiment_scratch_equivalence` |

The conditional pilot now samples actual order-finding outputs, but retains
exponential work-space setup/storage; it is distinct from the earlier
pre-QFT observable contraction. Claims C47–C50 separate the proved identities
from finite-instance measurements and remaining limitations.

**Consolidated papers and reachable-state follow-up (2026-09-10).** The current
synthesis is in the existing `PAPER_A.md` and `PAPER_B.md`, including the
distinction between full support, clean-code predictions, tensor rank and
conditional sampling. TODO 16 / C51–C52 / `notes/RO-reachable-orbits-and-spectral-sampling.md`
record on-demand compiled replay without scratch-space tables and comparisons
against charged classical baselines. It removes scratch overhead, but does
not outperform those baselines. An explicitly order-informed latent-eigenphase
sampler needs no work-state vector; obtaining the order remains a separate cost.
TODO 14 contains the next controlled-perturbation question.

## Start here

**`CLAUDE.md`** — conventions, traps and the discipline this project runs on.
Auto-loaded for agents; read it first if you are a person.

**`HANDOFF.md`** — if you are picking this up fresh, read that second. It states
where things stand, the open threads and their ranking, how to get running, and
what not to redo.

`claims/INDEX.md` is what is known · `todo/INDEX.md` is what is next ·
`notes/INDEX.md` is how things were found. **All three, and `CLAIMS.md`, are
generated** — edit the individual files and run `tools/reindex.py`.
`RESTRUCTURE.md` explains why the repo is laid out this way.

`METHOD.md` records how this work was actually conducted — the research loop,
the failure modes hit, and the techniques that earned their keep. Worth reading
before extending anything, since several of the lessons cost real time.

`CLAIMS.md` is the single source of truth for claim statuses (C-numbers,
F-numbers, and the retracted/dead ones). Both papers point at it rather than
carrying ledgers of their own; each keeps identifiers out of its prose and
resolves them through its Appendix A claim map.

`NOTES.md` is the working record — conventions, verified results, retracted
claims, traps, and an honesty log of everything believed and then killed. Read
it before trusting any number. Several early findings (F1–F8) are retracted;
the status blocks say which.

## Environment

Python >=3.14 is declared by `pyproject.toml`; run everything from this
directory via `uv`:

```bash
uv run python test_core.py          # correctness gate -- run first, always
```

**Use 3.12 or 3.13 for long `perm_pps` jobs.** CPython 3.14 dies intermittently
on them with `Fatal Python error: _TAIL_CALL_CACHE` — an interpreter bug, not a
bug in this code. It presents as a hang or a death with no traceback.

The project environment declares NumPy, Qiskit, Quimb, SciPy, Numba, Matplotlib,
Stim, Click, and `cupy-cuda13x`; `uv sync` installs that declared set even when runs
use the CPU. Set `LAB_GPU=1` to opt into CUDA when a compatible GPU is present
(7–15× on the large sweeps; see `accel.py`). For long `perm_pps` jobs, use
Python 3.12 or 3.13 as documented in `CLAUDE.md` because CPython 3.14 has an
intermittent tail-call-cache failure. For NumPy-only experiments, an isolated
supported invocation is:

```bash
uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_control_trace --extended-width 16
```

## Test suites (all must pass)

```bash
uv run python test_core.py           # Pauli algebra, gate decompositions, PPS vs dense
uv run python test_modexp.py         # Beauregard modexp, layers A-G
uv run python test_toffoli_arith.py  # Toffoli modexp, layers A-G + cross-check
uv run python test_walsh.py          # logical trace, Walsh machinery
uv run python test_perm_pps.py       # permutation-native PPS vs Walsh and vs PPS
uv run python test_windowed.py       # windowed modexp arithmetic + tail-block structure
uv run python test_accel.py          # CUDA backend vs CPU reference (skips without a card)
uv run python test_lab.py            # lab/ engine, pinned to logged numbers
uv run python test_claims.py         # CLAIMS.md headline rows, re-verified
```

Documentation has its own gate, run before committing doc changes:

```bash
uv run python tools/reindex.py       # regenerate the generated files
uv run python tools/check.py         # 10 consistency checks; must pass
```

Sanity anchor: `ModExp(15, 7, n_exp=4).build_shor()` must give exponent-register
peaks at y = 0, 4, 8, 12 with p = 0.25 each (r = 4).

## Layout

| file | role |
|---|---|
| `pauli.py` | symplectic Pauli algebra, the single rotation/conjugation rule |
| `circuits.py` | `Circuit` as (σ,θ) list + logical trace; gate decompositions; Cuccaro adder; QFT |
| `pps.py` | rotation-level Pauli propagation with δ-truncation |
| `perm_pps.py` | permutation-native PPS — stays Z-type throughout, ~2× lower peak in the measured arithmetic rows, ~10× faster |
| `walsh.py` | permutation extraction, FWHT, pullback coefficients, sparsity, affineness |
| `statevec.py` | O(2ⁿ)-per-gate state-vector simulation (batched) |
| `modexp.py` | Beauregard (Fourier-arithmetic) modular exponentiation |
| `PAPER_A.md` | **the Paper A draft** — full paper, supersedes `ABSTRACT.md` |
| `PAPER_B.md` | **the Paper B draft** — full paper, supersedes `ABSTRACT_SHOR_2ADIC.md` |
| `toffoli_arith.py` | Toffoli-compiled modular exponentiation |
| `accel.py` | optional CUDA backend for permutation replay and FWHT; opt-in via `LAB_GPU=1`, gated by `test_accel.py` |
| `windowed_arith.py` | windowed modexp: `WindowedModExp` (table lookup) and `SelectModExp` (select-multiply, with the `skip_zero` knob); `replay` for circuits too wide to hold a permutation array |
| `CLAUDE.md` | conventions, traps, research discipline — auto-loaded |
| `todo/` | one item per file; `open/` vs `done/` is the state. `INDEX.md` generated |
| `claims/` | one claim per file; `retracted/` may not be cited. `INDEX.md` and `CLAIMS.md` generated |
| `notes/` | one investigation per file, keyed by its §code. `INDEX.md` generated |
| `tools/` | `reindex.py` regenerates the generated files; `check.py` is the doc gate |
| `lab/` | experiment engine: protocol harness, cached measurement, GF(2)/structure analysis, modexp variants, null models |
| `experiments/` | the experiment scripts — lab-notebook records, filenames unchanged; run as `uv run python -m experiments.<name>`; start new ones from `TEMPLATE.py` |
| `archive/` | fully retracted scripts, kept for the record |
| `out/` | gitignored: logs and the measurement cache |

Notable experiments: `experiment_c8` (the Walsh identity, Paper A core),
`experiment_c7` + `experiment_c7_scale` (scaling to 30 qubits, density → ½),
`experiment_c15{,b}` (the 2-adic result, Paper B core), `experiment_c21_onset`
(the invariance onset measured at α = 3 and 4), `experiment_c17_deficit`
(the measured peak-ratio regularity, C44), `experiment_resid{1,2}` +
`experiment_affstruct` (the conditional-structure resolution, C33–C35).
`experiment.py` and `experiment2.py` are partially retracted — see their
headers.

## Scope and limits

- The identity covers **computational-basis (Z-type) observables on permutation
  circuits**. It fails for X/Y-type observables, where the pullback is not
  diagonal.
- It gives an *exact cost model*, not a faster simulator: the Walsh transform is
  itself O(2ⁿ).
- Nothing here bears on the classical hardness of factoring. Efficient
  simulation of Shor's algorithm on general inputs would be a classical
  factoring algorithm; these are diagnostic results about where PPS breaks.
