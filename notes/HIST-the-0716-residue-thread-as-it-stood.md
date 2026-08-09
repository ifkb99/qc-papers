---
code: HIST
title: "The 0.716-residue thread, as it stood at handoff (RESOLVED)"
outcome: solved
claims: [C30, C33, C34, C35]
todo: [11]
---

# HISTORICAL — the former live thread (RESOLVED; kept for context)

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
