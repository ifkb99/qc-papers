# Handoff — pick up here

Written for a session with **no prior context**. Read this file, then
`PAPER_A.md` / `PAPER_B.md` if writing, or `TODO.md` if experimenting.
Everything is committed; working tree clean.

**Branch: `todo-12e-gpu-sweeps`**, three commits ahead of `main` (`10a9288`) and
not pushed. Nothing on it is speculative — all nine test suites pass — so it can
be fast-forwarded into `main` whenever wanted.

---

## → START HERE: REFINING PAPER A AND PAPER B (in progress)

The experimental backlog is done. **The current job is editorial, not
experimental.**

**The two open questions were answered on 2026-08-08 and no longer need asking:**

1. **Target venue / format.** Shape both papers **as if submitting to a
   journal** — focused, every paragraph earning its place — but the likely
   first destination is **arXiv**. The user has not published before.
2. **Order.** **A consistency pass across both first, then Paper A end to end**,
   then Paper B. Paper B uses Paper A's identity as Lemma 1, so changes to A
   propagate into B.
3. **The third paper** (compilation-dependence, C30–C41) stays **scattered as
   caveats inside A and B** for now; revisit after the first two are submitted.

**The consistency pass is DONE** (commits `184f49c` and its successor). It found
**sixteen** defects, mostly 2026-08-08 edits applied in one place but not their
mirror — see those commit messages for the itemised list. Table 1's arithmetic,
all section cross-references, every claim ID against the ledger, and every
reproducibility command were verified mechanically rather than by eye.

**Remaining: Paper A end to end.** The three things flagged as worth attention:

- **§9 "Demonstration of reach"** is eight lines that say "there is a companion
  paper". At journal length that is a sentence in §1 or §12, not a section.
- **§11.1 / §11.2** are strong and load-bearing, but are where a length-conscious
  editor points first. Decide deliberately whether they stay in-body.
- **§4.4 Table 1** carries the paper, but its rot-PPS column stops at 17 qubits
  while the text claims 30. That asymmetry is disclosed honestly in §11.1 —
  check it is also disclosed *at the table*, not only 300 lines later.

**What changed in the papers that day, so a fresh session does not re-derive it:**

- `PAPER_A.md` §5 was rewritten from "the factor of two is an observation, and
  we state it as one" to **Proposition 2** (`rot = 2·perm − |B|`), plus the
  identification of B and the observable-dependence paragraph. The **abstract**
  was updated to match — it had still been claiming "a measured factor of 2.000
  to 1.9997".
- `PAPER_A.md` §11.1's instance sizes went 24 → **30 qubits**, and §4.4's
  closing paragraph with it.
- `PAPER_B.md` §6 gained the measured onset table for α = 3, 4 and the
  "invisible at a single width" caution; §7.1 gained C43 and the |I| = 1
  vacuity caveat; §9's density series went to 30 qubits / slope 1.006; §12.1's
  scope caveat "the α = 3 and α = 4 onsets are consistent-with but not
  confirmed" was **removed because they are now confirmed**; §12.3 gained the
  two new reproducibility lines.

**Known rough edges to look at while refining** (none are errors, all are
things a referee would poke):

- Paper A §5's Proposition 2 is derived, but **|B| = 2 itself is measured, not
  proved** — the paper says so; keep it that way unless someone proves it.
- Paper B §6's table is six instances; α ≥ 5 is genuinely out of reach (needs
  r = 32, whose smallest instance is q = 32) and the text says so.
- The **third paper** flagged below (compilation-dependence, C30–C41) is still
  only a possibility, and its material is deliberately left scattered as caveats
  inside A and B (decided 2026-08-08). A consequence worth tracking: **C33, C34
  and C35 are now cited in neither paper** — three established claims with no
  home. That is the Paper C material, and it is where to start if it is ever
  extracted.
- **Both papers now keep claim IDs out of the prose**, resolving them through an
  Appendix A claim map. Paper B was converted to match Paper A on 2026-08-08; it
  previously carried IDs in its section headings. Keep new text to that
  convention.
- `ABSTRACT.md` / `ABSTRACT_SHOR_2ADIC.md` are abstract workshops and prior-art
  dossiers, **not** superseded copies to edit in parallel. Prior-art text lives
  there; do not duplicate it into the papers.

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

**No live research thread is blocked.** TODO items 1–12b **and 12e** are
closed; 12c, 12f, 13 and 14 are open and each states its own rationale.

> **12e closed 2026-08-08 — see `NOTES.md` §OS.** Every sweep that had been
> deferred for want of compute is now run. C21's onset is **measured** at
> α = 3 (two moduli) and α = 4, not merely predicted; C24 strengthened to set
> level (**C43**); C7's circuit series reaches **30 qubits** (|S| =
> 536,271,623, slope 1.006 bits/qubit); C38's 2-periodicity holds to K = 5.
> `PAPER_A.md` §11.1 and `PAPER_B.md` §6 and §12.1 were updated accordingly.
> One method note is worth carrying: the must-fail control caught a **vacuous**
> test of mine (the tail-confinement check has only one possible answer at
> |I| = 1), the second time that rule specifically has earned its place.

> **12d closed too, same day — see `NOTES.md` §PK.** The peak-ratio deficit is
> **not a constant**: `N_max^rot = 2·N_max^perm − |B|` with
> B = {z ∈ S : z_c = 0}, c the target of the gadget containing the peak. Derived
> from the gadget (the four T gates on c rotate inside the *closed* space
> span{X_c, Y_c}, so they branch once between them, not 2⁴ times) and verified
> 9/9. It is 2 for modexp because B = {Z_x0, Z_x0·Z_e0} — **exactly the two
> |c| = ½ Walsh coefficients**, the dominant Fourier modes §W3/W4 had already
> found from the final spectrum by an unrelated route. And the constant belongs
> to the **observable**, not to modexp: move it and you get 4004, 4014, or 0.
> Paper A §5 is now Proposition 2 rather than "an observation"; claim **C44**.
> Still unproved and marked so: |B| = 2 itself.

> **12g closed negatively — see `NOTES.md` §DF.** The one loose end 12e left
> (a suspected "cliff" in the hyperplane deficit between n_exp 2 and 3) does
> not exist: sweeping n_exp by 1 at fixed (N, a) over 10–12 widths gives a
> **monotone** decrease every time, and the cliff was two points straddling the
> single step at which N = 21 disagrees with every other instance. §OS4 is
> retracted at its source and in the honesty log. The §I period-ord₂(β)
> structure is function-level only and **does not transfer to circuit level** —
> that transfer failing is the reusable part. The experiment fails 3 of its 4
> predictions and exits nonzero by design.

The best remaining candidate is **14** (through the inverse QFT — highest risk,
highest reach, a Paper C candidate). **12c** is mostly dead (its naive form was
refuted by derivation) and **12f** is deliberately not done.

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
- **The binding constraint is the REPLAY, not the transform** (found 2026-08-08
  while doing TODO 12e). `idx ^= ((idx >> c) & 1) << t` keeps the array plus two
  temporaries live, so an int64 index array needs ~24 GiB at n = 30 and fails.
  Images are < 2ⁿ, so `accel._replay` builds it in **int32** for n ≤ 30 (~12
  GiB), which is what makes a q = 30 circuit-level run possible at all. Gated
  against the CPU int64 reference in `test_accel.py` [A]. A real q = 30 modexp
  now takes ~14 min end to end.
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
  once, one per device, via `CUDA_VISIBLE_DEVICES=0` / `=1`. Zero new code —
  and it is how 12e was actually run, two sweeps in parallel throughout. The
  measurement cache is content-addressed, so concurrent writers cannot collide
  and a later single-process run of the experiment replays everything for free.
  That is the pattern to reuse: **warm the cache in parallel, then run the
  experiment file once for the record.**
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
6. **Do not pipe a long experiment through `tail`** — the pipe buffers the whole
   run, so a job that streams progress for 40 minutes shows nothing until it
   exits, and looks hung. Redirect to a file and read that.
7. **Do not build large `frozenset`s of Python tuples over a support.** At 33.5M
   elements that is a multi-GB structure and an OOM; encode the pair as one
   int64 and compare sorted numpy arrays instead (see `experiment_c21_onset`'s
   signature construction). Same answer, ~40× less memory.

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

### Concrete things to try — **ALL FIVE ARE DONE. This is not a to-do list.**

Resolved as TODO item 11 (C33–C35, `NOTES.md` §AF/§RS): affine structures found
the mechanism, the support complement *is* a recognisable set (the empty
quadrant), the weight profile was mooted, 0.716 turned out instance- and
width-dependent rather than a constant, and stacking showed the residue tracks
independent nonlinear monomials rather than wrap count.

<details><summary>The original five, kept for the record</summary>

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

</details>

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

**Added 2026-08-08 — reach for these before writing new machinery:**

- `lab.measure.support(qc, q, exact=True)` — the exact integer path (support
  test is `!= 0`, no tolerance). Cached under its own key so exact and
  thresholded results can be compared rather than silently swapped.
- `lab.measure.stats(qc, q, masks=(w,))` — `count`, `density`, and per-mask
  GF(2) parity counts **without moving the support to the host**. At q = 30 a
  support is 4 GiB; a density sweep wants three scalars. `odd[w] == 0` is
  exactly "w is a linear structure" (C30). Cached as JSON.
- `accel.pullback_stats` / `accel._replay` — the backend for both. See the GPU
  section for why the replay is int32.
- **Peak-state extraction without a new propagator** (`experiment_c17_deficit`):
  Heisenberg propagation of the *last m gates* IS the state after m steps, so
  a peak Pauli set is obtained by running the existing verified `pps.propagate`
  on a gate suffix, and a perm-level state by running `propagate_perm` on a
  logical suffix. Writing a second propagator to inspect the first is exactly
  how this project has been bitten before; do not.

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

It is not decoration. Its rules have now caught **six** of my own errors that
had already produced confident-looking numbers:

- three **vacuous tests** — step 9's blocks acting only on unobserved qubits;
  TODO 12e's `z_I ∈ {0, 1_I}` check at |I| = 1; TODO 12g's sign test on
  monotone data. Two were caught by a must-fail control failing to fail; the
  third by noticing the pass was unanimous and unearned. **A unanimous pass is
  a tell, not a triumph.**
- one **control whose predicate was wrong** (12g's β = 1 control demanded
  monotone growth including a pre-lock step) — the control was right, the code
  asking it was not, and fixing it produced a sharper control;
- a sweep that varied two parameters at once;
- degenerate random inputs at small sizes.

The two rules that pay most, by count: **every experiment carries a must-fail
control**, and **derive before measuring** (which also means the prediction
must be written where it can be seen to have come first — `lab.harness`
enforces both).

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

The four experiments written on 2026-08-08, all runnable and all carrying their
own OUTCOME block in the header — read the header, not just the code:

```
LAB_GPU=1 uv run python -m experiments.experiment_c21_onset      # 7/7, §OS
LAB_GPU=1 uv run python -m experiments.experiment_c7_scale       # 6/6, §OS4
LAB_GPU=1 uv run python -m experiments.experiment_windowed_scale # 4/4, §OS6
         uv run python -m experiments.experiment_c17_deficit     # 6/6, §PK
LAB_GPU=1 uv run python -m experiments.experiment_c7_deficit     # 2/5, §DF
```

**`experiment_c7_deficit` exits NONZERO ON PURPOSE** — three of its predictions
are refuted and it is left failing rather than re-scoped, the same treatment
`experiment_windowed.py`'s P5 got. A red line in a log is not always a bug
here; check the header first. Everything is cached in `out/cache/` (gitignored,
~2 GB), so a re-run replays in seconds rather than hours.

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
- **A "cliff" in the hyperplane deficit between n_exp = 2 and 3.** Refuted the
  same day it was conjectured (§DF). D = 1 − 2·density is **monotone
  decreasing** in n_exp in every β > 1 instance over 10–12 widths; the apparent
  cliff was a two-point comparison landing on the one step where N = 21
  disagrees with every other instance.
- **The §I period-ord₂(β) oscillation at CIRCUIT level.** It is a
  **function-level** phenomenon only and does not transfer. Do not import it
  again — that import is what made 12g look worth a sweep.
- **The identity string Z⁰, or a parity constraint, as the peak deficit.** Both
  were the named candidates in TODO 12d and both are wrong; the answer is the
  two |c| = ½ Walsh modes (§PK).

---

## Remaining TODO items (see `TODO.md` for the ranked list)

**Open, in the order they are worth doing:**

- **14 — through the inverse QFT.** The only genuinely open frontier, and a
  Paper C candidate. Highest risk, highest reach. Note the standing warning in
  `TODO.md`: **Cîrstoiu was pulled for this and does not help** — their group
  indexes *circuit parameters* and their results are about *ensembles*, while
  Shor's circuit is fixed. Do not re-pull it for this.
- **12c — early pruning of the dead 2^−(α+1) fraction.** Its naive form is
  already dead by derivation (exponent support is not monotone under
  back-propagation, proved with an explicit `CCX` counterexample). What
  survives is a much harder question about certificates.
- **13 — a third simulation method (DDSIM) on the r = β·2^α invariant.** Mostly
  integration work, no new theory; strengthens "property of the algorithm, not
  the simulator".
- **12f — two-GPU split FWHT.** Deliberately **not** done: the second card buys
  exactly one qubit and `wht_exact` already buys the same one for free.
- **Housekeeping:** six unused dependencies (`click`, `matplotlib`, `numba`,
  `quimb`, `scipy`, plus `qiskit`/`stim` which were the source-reading tools for
  the prior-art sweeps). Only `numpy` and `cupy` are imported anywhere. Decide.

**Longer-standing, still open:**

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

**Paper A is drafted in full: `PAPER_A.md` (v2, revised 2026-08-08).** It
supersedes `ABSTRACT.md` as the live document; `ABSTRACT.md` is retained as the
abstract workshop and prior-art dossier. **Paper B is drafted in full too:
`PAPER_B.md` (v1, revised 2026-08-08)**, superseding `ABSTRACT_SHOR_2ADIC.md`,
which is retained as the abstract workshop. **Refining both is the agreed next
task — see the START HERE section at the top for what changed and what is still
open.** A third paper is now viable and was not before — the
compilation-dependence results (C30–C32, C33–C35, C36–C39) acquired a spine in
C40/C41 and are currently scattered as caveats inside A and B.

**Claims added 2026-08-08: C43** (C24 holds at set level, not merely in
cardinality) and **C44** (the peak ratio is exact: `rot = 2·perm − |B|`).
**Regraded:** C21 to α = 1..4, C7 to 30 qubits, C17 from "the factor is
empirical" to derived, C30 re-verified at q = 30, C37/C38 extended to K = 5.
**Removed:** §OS4's hyperplane-deficit "cliff". `CLAIMS.md` is canonical.

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
