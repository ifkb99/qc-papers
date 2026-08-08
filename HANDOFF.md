# Handoff — pick up here

Written for a session with **no prior context**. Read this first, then
`NOTES.md`. Everything is committed; working tree clean at 16 commits.

---

## Where things stand in one paragraph

Two papers' worth of results, both with their central claims **proved**, not
merely observed. Paper A (`ABSTRACT.md`): for any circuit implementing a
permutation of the computational basis, the Pauli support carried by Pauli-path
simulation with a computational-basis observable is *exactly* the
Walsh–Hadamard spectrum of the corresponding output-bit Boolean function — an
exact cost model where the field currently uses empirical extrapolation. Paper B
(`ABSTRACT_SHOR_2ADIC.md`): applied to modular exponentiation, cost is governed
by the 2-adic structure of the multiplicative order r; writing r = β·2^α with β
odd, cost is independent of exponent-register width when β=1 and Θ(2ⁿ)
otherwise. Twelve of fifteen TODO items are closed. **One live thread**, below.

---

## THE LIVE THREAD — residual non-linear structure

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

- `walsh.py` — `pullback_coefficients`, `classical_permutation`, FWHT,
  `walsh_sparsity`, `is_affine`.
- `experiment_linstruct.py` — GF(2) rank + kernel extraction (reuse
  `gf2_rank_and_kernel`).
- `experiment_reduction2.py` — the `Variant` subclass that injects modified
  reductions while preserving correctness. Add new modes to `_wrap`.
- `perm_pps.py` — permutation-native propagation, ~10× faster than
  rotation-level, with `delta` and `max_weight` truncation.

---

## Getting running (2 minutes)

```bash
cd /home/djneko/Workspace/qsim-test/QuantumSimTest-master/research
uv run python test_core.py          # correctness gate — run first, always
uv run python test_perm_pps.py
```

All five suites must pass before trusting anything:
`test_core`, `test_modexp`, `test_toffoli_arith`, `test_walsh`, `test_perm_pps`.

**Traps that will bite immediately:**

1. **Run with `uv run python` from `research/`.** Not bare `python3`.
2. **`source .env` ONLY for Julia work.** Its `LD_PRELOAD` of Julia's libstdc++
   **segfaults numpy `longdouble`** (exit 139). Cost a debugging cycle.
3. **CPython 3.14 crashes on long `perm_pps` runs** —
   `Fatal Python error: _TAIL_CALL_CACHE`. A 3.14 interpreter bug, not this
   code. Presents as a hang or a mysterious death with no traceback. Rerun (it
   is intermittent) or use 3.12/3.13 for long jobs.
4. **Do not `pkill` in the same command as a heredoc write** — it kills the
   write. Two scripts vanished this way.

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
7. The two abstracts, when writing. They carry no ledgers of their own; both
   point at `CLAIMS.md`.

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

## Remaining TODO items besides the live thread

- **Exact-spectrum crypto import.** The bound S ≥ (1 − NL/2ⁿ⁻¹)⁻² is tight only
  at the extremes (AES: bound 64, actual 239). For crypto families whose *full*
  Walsh value/multiplicity distribution is published, S is determined **exactly**
  rather than bounded — a much stronger import.
- **Does windowed / table-lookup arithmetic satisfy the involution criterion?**
  The C15 theorem now covers any construction whose a=1 block is an involution
  (V²=id is the entire condition). Gidney-style windowed arithmetic is the
  natural test; it is a well-posed question, not a survey.
- **Yao.jl source-level check.** Closed as far as possible without installing
  it; its docs list Toffoli under "Clifford Gates: Two-qubit gates", wrong on
  both counts. Low value.

---

## Publication state

Both abstracts are drafted with honest status markers and scope caveats; the
claims ledger they both point at is `CLAIMS.md`. Before submitting anything:

- The **prior-art position is settled and narrowed**: the diagonal↔Walsh
  ingredient is standard ([arXiv:1306.3991](https://arxiv.org/pdf/1306.3991))
  and the Pauli-spectrum↔Fourier analogy is established
  ([arXiv:2311.09631](https://arxiv.org/pdf/2311.09631)). The contribution is
  their composition into an exact PPS cost model for permutation circuits, plus
  what follows. Quipu/stabilizer frames and Cîrstoiu must be cited and
  distinguished — see `ABSTRACT.md`.
- Keep the honesty notes in. They are load-bearing, not decoration.
