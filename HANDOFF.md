# Handoff — pick up here

Written for a session with **no prior context**. Read this file, then
`PAPER_A.md` / `PAPER_B.md` if writing, or `TODO.md` if experimenting.
Everything is committed; working tree clean.

---

## Where things stand

**Both papers are drafted in full.** `PAPER_A.md` (the exact cost model) and
`PAPER_B.md` (the 2-adic result for modular exponentiation) supersede
`ABSTRACT.md` and `ABSTRACT_SHOR_2ADIC.md`, which are retained as abstract
workshops and prior-art dossiers. `CLAIMS.md` remains the single source of truth
for claim *status*; where a paper and the ledger disagree, the ledger wins.

Paper A: for any circuit implementing a permutation of the computational basis,
the Pauli support carried by Pauli-path simulation with a computational-basis
observable is *exactly* the Walsh–Hadamard spectrum of the corresponding
output-bit Boolean function — an exact cost model where the field uses empirical
extrapolation. Paper B: applied to modular exponentiation, cost is governed by
the 2-adic structure of the multiplicative order r = β·2^α — independent of
exponent-register width when β = 1, Θ(2ⁿ) otherwise, with the mechanism proved.

**No live research thread is blocked.** TODO items 1–12b are closed; 12c–12e,
13 and 14 are open and each states its own rationale. The highest-value open
item is **12e**, and specifically **C21's onset at α = 3 and 4**: a
pre-registered prediction that went unmeasured only because it was too slow,
and is now cheap (see GPU, below).

---

## GPU — read this before running anything slow

The two hot paths (permutation replay, FWHT) are memory-bound array passes, and
a CUDA backend gives **7–15×, growing with n**. This is already built, gated and
wired in. It will save hours.

```bash
LAB_GPU=1 uv run python -m experiments.<name>     # that is the whole interface
```

- **`accel.py`** is the backend; **`test_accel.py`** (suite 9) gates it against
  the CPU reference and skips cleanly with no card.
- **Opt-in by design.** `walsh.py` stays the reference and is never silently
  substituted. `LAB_GPU=1` routes `lab.measure.support` through the GPU; nothing
  else changes behaviour.
- **Measured** (RTX A4500): permutation replay 7.1× / 6.9× / **11.5×** at
  q = 17 / 19 / 21 (29.3 s → 2.6 s at q = 21); FWHT 11.8× / 12.8× / **14.7×** at
  2²⁴ / 2²⁶ / 2²⁸ (18.1 s → 1.2 s at 2²⁸).
- **Capacity: n ≤ 30** on one 20 GiB card (`accel.MAX_QUBITS`); it raises rather
  than thrashing past that. Measured: float64 FWHT at n = 30 needs 12 GiB and
  fits; n = 31 needs 24 GiB and does not.
- **Does the second card scale it further? Only by +1 qubit, and there is a
  cheaper way.** Memory doubles per qubit, so 40 GiB buys exactly one more than
  20 GiB. A split would be genuinely easy — with the array halved by its top
  bit, every FWHT level except the last is local to a half, and only the final
  butterfly crosses cards — but +1 qubit is a poor return, so it is **not
  implemented** (logged as TODO 12f).
  **Use `wht_exact` / `pullback_support_exact` instead:** the FWHT of ±1 data is
  *exactly integer-valued* (verified bit-for-bit), so int32 gives the same
  answer in half the memory — the same +1 qubit, no complexity — **and it makes
  the support test exact (≠ 0) rather than a magnitude threshold.** That second
  property matters: it removes the thresholding artifact that otherwise makes
  measured densities drift below their true value as n grows (see `NOTES.md`
  §GF and `experiment_gf2law_scale.py`). Safe to n = 30, since intermediate
  magnitudes are bounded by 2ⁿ and 2³⁰ < 2³¹.
- The best use of the second card is **throughput**: two independent sweeps at
  once, one per device, via `CUDA_VISIBLE_DEVICES=0` / `=1`. Zero new code, and
  it is what TODO 12e actually wants.
- **Install is machine-specific.** `pyproject.toml` pins `cupy-cuda13x` to match
  this machine's CUDA 13.3. On a CUDA 12 host swap to `cupy-cuda12x`; both were
  tested and perform identically. Without cupy everything still works on CPU.
- `accel.device_info()` prints what was detected; `accel.enabled()` tells you
  whether `LAB_GPU` actually took effect.

---

## Getting running (2 minutes)

```bash
cd /home/djneko/Workspace/qsim-test/QuantumSimTest-master/research
uv run python test_core.py          # correctness gate — run first, always
uv run python test_claims.py        # headline results, pinned to logged numbers
```

**Traps that will bite immediately:**

1. **Run with `uv run python` from `research/`.** Not bare `python3`.
2. ~~`source .env` only for Julia work~~ — **TRAP RETIRED 2026-08-08.** Julia
   was unused (no `.py` imported `juliacall`); `juliacall`, the 1.1 GB depot and
   `.env` are all deleted, so the `LD_PRELOAD`/`longdouble` segfault cannot
   recur. The venv went 1.6 GB → 534 MB.
3. **CPython 3.14 crashes on long `perm_pps` runs** —
   `Fatal Python error: _TAIL_CALL_CACHE`. A 3.14 interpreter bug, not this
   code. Presents as a hang or a mysterious death with no traceback. Rerun (it
   is intermittent) or use 3.12/3.13 for long jobs.
4. **Do not `pkill` in the same command as a heredoc write** — it kills the
   write. Two scripts vanished this way.
5. **A GPU run needs `LAB_GPU=1`** — `accel.enabled()` returns False without
   it and everything silently runs on CPU at 1/10th the speed.

---

## HISTORICAL — the former live thread (RESOLVED; kept for context)

> **RESOLVED 2026-08-08 (post-handoff).** The residue is a **conditional
> linear structure**: the support exactly avoids the quadrant
> {z_msb=1, z_anc=0}, capping density at ¾. See `NOTES.md` §AF/§RS, claims
> C33–C35 in `CLAIMS.md`, and TODO item 11. The section below is kept as the
> state of knowledge at handoff time.

### What is known

Modexp pullbacks carry a **linear structure** w = b_msb ⊕ anc: the GF(2) rank of
the Walsh support is n−1, so the support lies in a hyperplane and density is
capped at exactly ½. This explains why C7 saw density converge to 0.498 *from
below* and never cross. Verified pointwise: `g(y ⊕ w) = g(y)`.

It is forced by three things acting together (§L, §L2 of `NOTES.md`):
(a) flipping the msb *is* adding 2^(m−1), which commutes with mod-2^m addition;
(b) anc is coupled to msb only by XOR; (c) the msb is `b[n]` while the cswaps use
`b[:n]`, so it never reaches the observed register.

Breaking it needs **nonlinearity in the msb**. Conjugating the reduction with
`toffoli(msb, t0, anc)` destroys it — full rank, density 0.473 → 0.716 — while
still computing `a^e mod N` correctly.

### The open question

**The broken variants sit at 0.716–0.721, not 1.000. Random Boolean functions
reach 1.000. So something non-linear survives after the linear structure is
destroyed, and it is unidentified.**

### Concrete things to try

1. **Is it another linear structure at higher order?** Check for *affine*
   structures (w with `g(y⊕w) = g(y) ⊕ const`, not just `= g(y)`) — these also
   constrain the support but are not caught by the rank test currently used.
2. **Is the missing ~28% structured?** Look at the *complement* of the support.
   If the absent z form a recognisable set (a coset, a weight band, a subspace
   union), that names the constraint.
3. **Weight profile.** Compare mass-by-weight of the broken variant against the
   binomial baseline (machinery exists in `experiment_weight.py`). Deviation
   would localise the structure by Fourier degree.
4. **Is 0.716 a recognisable constant?** It was 23464/32768 at N=5 and
   23488/32768 at N=7 — close but not equal, so probably not an exact rational.
   Worth ruling in or out early.
5. **Does more nonlinearity push it to 1.000?** Stack several independent
   Toffoli conjugations. If density saturates below 1, the residue is intrinsic
   to the arithmetic rather than to the reduction.

### Where the tools are

> **Restructured 2026-08-08.** The patterns below now live in the `lab/`
> package — `lab.gf2` (rank/kernel plus the affine-aware structure finder),
> `lab.variants` (the wrap registry that replaced the copy-pasted `Variant`
> subclasses, with `verify_correctness`), `lab.measure` (cached
> support/density/peak), `lab.modarith`, `lab.nulls`, and `lab.harness`
> (`Experiment`: predictions-before-measurement and must-fail controls,
> enforced). Start new experiments from `experiments/TEMPLATE.py`. Finished
> scripts moved to `experiments/` unchanged; run them as
> `uv run python -m experiments.<name>`.

- `walsh.py` — `pullback_coefficients`, `classical_permutation`, FWHT,
  `walsh_sparsity`, `is_affine`.
- `experiments/experiment_linstruct.py` — GF(2) rank + kernel extraction
  (now `lab.gf2.rank_kernel`).
- `experiments/experiment_reduction2.py` — the `Variant` subclass that injects
  modified reductions while preserving correctness (now `lab.variants`).
- `perm_pps.py` — permutation-native propagation, ~10× faster than
  rotation-level, with `delta` and `max_weight` truncation.
- `windowed_arith.py` — the two windowed modexp constructions (`WindowedModExp`
  table-lookup, `SelectModExp` select-multiply with the `skip_zero` knob), plus
  `replay` (single-basis-state image, for circuits too wide to hold a
  permutation array) and `verify_modexp`. Gated by `test_windowed.py`.

---

## The `computational-research` skill

There is a skill at `~/.claude/skills/computational-research/SKILL.md`
generalising the working protocol used here. It was **created mid-session, so it
was not loadable at the time** — a fresh session should be able to invoke it
with `Skill(computational-research)`. Worth loading before doing experimental
work on this project.

It encodes the one thing that most distinguishes this domain: **in empirical
science a surprising result may be a discovery, but in research on your own code
it is almost always your own bug** — and the corollary that a result *agreeing*
with your hypothesis is more dangerous than one contradicting it, because you
won't look. Also: a bug-check procedure ordered by cost, derive-then-test,
always including a control that must fail, promoting regularities to proofs by
reading the construction rather than measuring more, and explicit claim grading.

It is not decoration. In this session its rules caught three of my own errors
that had already produced confident-looking numbers: a vacuous test whose
must-fail control failed to fail, a sweep that varied two parameters at once,
and degenerate random inputs at small sizes.

## Reading order

1. **This file.**
2. The `computational-research` skill (above), if doing experimental work.
3. `METHOD.md` — how the work is conducted, the failure modes hit, and the
   techniques that earned their keep. Short. Read before extending anything;
   several lessons cost real time. The skill is its generalisation; this file is
   the project-specific version with the concrete instances.
4. `NOTES.md` — the working record. **Read the STATUS blocks and the honesty log
   before trusting any number.** Roughly half of what was believed at various
   points is now marked dead, with the reason.
5. `TODO.md` — ranked, with rationale, so it can be re-ranked rather than
   followed blindly.
6. `CLAIMS.md` — the canonical ledger of every claim ID and its status, split
   into Paper A, Paper B and retracted/dead. Read it before citing any claim:
   where `NOTES.md` prose disagrees, `CLAIMS.md` wins.
7. `PAPER_A.md` / `PAPER_B.md` — the full drafts, when writing. `ABSTRACT.md`
   and `ABSTRACT_SHOR_2ADIC.md` are demoted to abstract workshops and prior-art
   dossiers. None carries a ledger; all point at `CLAIMS.md`.
8. `accel.py` — before running anything large. See the GPU section above.

---

## Do NOT redo these

Logged in full in `NOTES.md`; listed here so a fresh session does not burn time.

- **The compilation thesis** ("Toffoli vs Fourier compilation determines cost").
  Dead twice over — a propagator bug plus an ancilla-count confound. Both
  compilations have exact Z-closure and the identity is compilation-invariant.
- **"The QFT is the PPS bottleneck."** Dead — bit-reversal bug plus a toy that
  was not a proxy for Shor.
- **The toy "Shor sandwich"** (H → adder → QFT). Not a proxy. Use `modexp.py`
  or `toffoli_arith.py`.
- **Affineness of the identity block** as a route to the C15 proof. Refuted
  (14336/32768 violations). The proof goes via the parity reduction instead.
- **The (z, z⊕e) pairing** as the C15 mechanism. An n_exp=1 artifact.
- **An intermediate 2-adic law in α.** There is none; the β=1/β>1 split is the
  real structure and the apparent intermediate case was sampling aliasing.
- **The msb↔anc CNOT pairing as the cause of the linear structure.** Refuted;
  see §L2 for the actual three-part mechanism.

---

## Remaining TODO items (see `TODO.md` for the ranked list)

- **Exact-spectrum crypto import.** The bound S ≥ (1 − NL/2ⁿ⁻¹)⁻² is tight only
  at the extremes (AES: bound 64, actual 239). For crypto families whose *full*
  Walsh value/multiplicity distribution is published, S is determined **exactly**
  rather than bounded — a much stronger import.
- ~~Does windowed / table-lookup arithmetic satisfy the involution criterion?~~
  **Closed 2026-08-08 — see the box at the top and `NOTES.md` §WD.** The
  answer forced C29 to be narrowed and produced C36–C39. What it leaves open:
  the cost is 2-periodic rather than constant, and nothing here says whether
  that survives Gidney's *measurement-based* unlookup, which is not unitary and
  so is outside every argument in this project so far.
- **Yao.jl source-level check.** Closed as far as possible without installing
  it; its docs list Toffoli under "Clifford Gates: Two-qubit gates", wrong on
  both counts. Low value.

---

## Publication state

**Paper A is drafted in full: `PAPER_A.md` (v1, 2026-08-08).** It supersedes
`ABSTRACT.md` as the live document; `ABSTRACT.md` is retained as the abstract
workshop and prior-art dossier. **Paper B is drafted in full too: `PAPER_B.md` (v1)**, superseding
`ABSTRACT_SHOR_2ADIC.md`, which is retained as the abstract workshop. A third paper is now viable and was not before —
the compilation-dependence results (C30–C32, C33–C35, C36–C39) acquired a
spine in C40/C41 and are currently scattered as caveats inside A and B.

Both abstracts carry honest status markers and scope caveats; the claims ledger
they all point at is `CLAIMS.md`. Before submitting anything:

- The **prior-art position is settled and narrowed**: the diagonal↔Walsh
  ingredient is standard ([arXiv:1306.3991](https://arxiv.org/pdf/1306.3991))
  and the Pauli-spectrum↔Fourier analogy is established
  ([arXiv:2311.09631](https://arxiv.org/pdf/2311.09631)). The contribution is
  their composition into an exact PPS cost model for permutation circuits, plus
  what follows. Quipu/stabilizer frames and Cîrstoiu must be cited and
  distinguished — see `ABSTRACT.md`.
- Keep the honesty notes in. They are load-bearing, not decoration.
