# Draft abstract — Paper B (Shor / 2-adic)

> **SUPERSEDED AS THE DRAFT 2026-08-08 → `PAPER_B.md`.** The full paper draft
> now lives in `PAPER_B.md`; this file is retained as the abstract workshop and
> the honesty-note dossier that `PAPER_B.md` §12 condenses. Where the two
> disagree on wording, `PAPER_B.md` is live; where either disagrees with
> `CLAIMS.md` on a *status*, `CLAIMS.md` wins.

> **Review update, 2026-09-09.** The β=1 theorem concerns the full operator;
> β>1 growth is empirical, scalar output periods can be smaller than r, and
> the proved onset is sufficient, not a universal first-lock lower bound.
> Useful fractions are approximate measurements (C42). C45/C46 now supply
> exact finished-control contraction with bounded reduced support in both
> order branches, so full-operator growth is not a lower bound on this task.

**Working title:** *The 2-Adic Structure of the Order Determines Pauli-Path
Simulation Cost for Modular Exponentiation*

Baker, I.

---

> **STATUS — read before circulating.** Split out of the Walsh cost-model paper
> (`ABSTRACT.md`, "Paper A") because the claim, the audience and the evidence
> are separable: Paper A is a methodology result about Pauli propagation on
> permutation circuits; this one is a result about Shor's algorithm specifically,
> and it *uses* Paper A's identity as an instrument. Publishable order is A then
> B, or B alone with A's identity cited as a lemma.
>
> Post-bugfix (θ=π propagator bug, see `NOTES.md` STATUS 2) and post external
> review. Confirmed for both final support and peak memory. The core dichotomy
> is verified in a controlled design at **three** moduli (N = 5, 7, 21); the
> onset rule n_exp = v₂(r)+1 is measured for α = 1, 2 and consistent for α = 3, 4.
> **The circuit-level mechanism is now PROVED** (parity reduction via an
> involution) — see the honesty note, which also records the scope limit.

---

## Abstract for the arXiv metadata field (revised 2026-09-09)

The PDF keeps the long abstract in `PAPER_B.md`; this compressed version is for
the submission form only. Plain ASCII because the field mangles math glyphs.

> Shor's algorithm is the canonical target for classical simulation studies of
> quantum advantage, its difficulty usually parameterised by qubit count and
> depth. For full-operator Pauli-path simulation, a dominant cost parameter is
> instead an arithmetic property of the instance: the 2-adic structure of the
> multiplicative order r of the chosen base. Writing r = beta 2^alpha with beta
> odd, we prove the Pauli support carried through the studied reversible
> modular-exponentiation construction is independent of the exponent-register
> width when beta = 1 once n_exp >= alpha + 1, and measure its first lock at
> that width in the reported rows; we measure its
> density bounded away from zero in sampled beta > 1 circuits -- consistent with
> Theta(2^n). In a controlled design varying only the base, the support at
> N = 7, a = 6 (r = 2) is exactly 15549 terms across a 256-fold growth in
> Hilbert-space dimension. The mechanism is proved: every identity-tail block
> applies the same involution, so the circuit depends on the whole tail through
> a single parity bit. The criterion this yields is checkable: Gidney-style
> windowed lookup arithmetic preserves the invariance in a strictly sharper
> form, while an obvious optimisation (skipping the trivial multiply) destroys
> it entirely. Since every order divides lambda(N), a modulus with lambda(N) a
> power of two lies in the free branch for every base: N = 15 is the smallest
> relevant odd semiprime, and all seven of its usable bases are free -- the
> canonical demonstration instance is degenerate because of the modulus, not an
> unlucky base. These are representation costs for a pre-QFT observable, not
> Shor output-sampling costs. Exact input-state contraction bounds reduced
> support independently of exponent width in both order branches.

---

## Abstract (draft)

Shor's algorithm is the canonical target for classical simulation studies of
quantum advantage, yet the difficulty of simulating it is usually discussed in
terms of qubit count and circuit depth. We show that for Pauli-path simulation
(also called sparse Pauli dynamics) the dominant cost parameter is instead an
arithmetic property of the instance: the 2-adic structure of the multiplicative
order r of the chosen base.

Our instrument is the exact correspondence, established separately, between the
Pauli support carried by a Heisenberg-picture propagation through a
basis-permutation circuit and the Walsh–Hadamard spectrum of the corresponding
output-bit Boolean function. Because a^e mod N depends only on e modulo r, and
because the Walsh basis is the character group of (Z/2)^t, the arithmetic of r
relative to powers of two controls the spectrum directly. Writing r = β·2^α with
β odd, the idealised bit function e ↦ bit_j(a^e mod N) is a function of the low α
bits of the exponent alone when β = 1, giving a spectrum of *constant* size
independent of exponent-register width. Any odd factor β > 1 does not by itself
force every selected bit to be dense: the relevant period is the selected
scalar function's minimal period (for N=13, a=4, r=6, the least-significant bit
is 1⊕e₀ and has sparsity 1). The β > 1 circuit branch is therefore an
empirical, observable-dependent result, with growth consistent with Θ(2^t) in
the measured family. At t = 24 the β=1 spectrum has 4 nonzero coefficients for
N = 15, a = 7 (r = 4). Density for β > 1 need not be
1 — it oscillates with period ord₂(β) and can sit at exactly one half: N = 7,
a = 3 (r = 6) gives exactly 0.500000 on even t, including t = 24, while N = 21,
a = 2 gives 1.000000 at t = 24.

We then verify the effect on complete, end-to-end verified reversible modular
exponentiation circuits in a controlled design that holds the modulus, ancilla
layout, gate structure and observable fixed and varies only the base, so that r
alone changes. For N = 7 with a = 6 (r = 2) the Pauli support is *exactly
invariant* at 15549 terms, and peak memory *exactly invariant* at 24369 terms,
across a 64- to 256-fold growth in Hilbert-space dimension; for a = 3 (r = 6) the
same quantities grow by a clean factor of 4.00 per two added qubits. The result
reproduces at N = 21 and N = 5. Consequently the exponent register — which sets
the precision of the continued-fractions post-processing, and which one would
expect to be the expensive resource — is free for Pauli-path simulation when r
is a power of two once the identity tail is present. In the reported β>1
controls, cost quadruples per two qubits; that branch is an empirical circuit
result rather than a universal odd-factor theorem.

The invariance has a proved sufficient threshold. Writing α = v₂(r), the
controlled-multiplier block attached to exponent bit i multiplies by a^(2^i)
mod N, which is the identity on the valid subspace exactly when i ≥ α; the
parity proof therefore guarantees locking for **n_exp ≥ α + 1**. Whether an
earlier width has already stabilised is observable-dependent. We
confirm this for α = 1 and α = 2 (where the support grows once, 15493 → 32143,
before freezing) and find α = 3 and α = 4 still growing at the largest width we
can reach, as predicted. The locked value is close to half the Hilbert-space
dimension at the lock point, so the support saturates to half density and then
freezes in absolute terms while density falls fourfold per two added qubits.

The invariance is not an artifact of exact arithmetic. Applying a coefficient
threshold to the finished spectrum preserves it at *every* threshold, since the
surviving magnitudes are identical across widths. Under the incremental
truncation that Pauli-path simulators actually perform, the peak term count —
the quantity that bounds memory — was likewise unchanged at every threshold
tested, an empirical result rather than a theorem for every threshold. Accuracy
is independent of width up to moderate thresholds but degrades
first for the wider circuit at aggressive ones, not because its exact spectrum
differs (it is identical) but because incremental truncation has more gates to
act upon. We therefore state the practical claim as: peak cost is free in the
exponent register, and accuracy is free up to moderate truncation.

For instances with an odd factor, which is the generic and cryptographically
relevant case, we find the Walsh density converges monotonically to one half
(0.473 → 0.498 over 15 to 24 qubits, growth 1.008 bits per qubit), so
Pauli-path cost is Θ(2^n): asymptotically no better than state-vector
simulation. The limiting constant of one half is not accidental and not
algorithmic: the support carries a linear structure w = b_msb ⊕ anc, pairing the
accumulator's sign bit with the modular-reduction comparison ancilla, which
confines it to a hyperplane. We verify this is observable-independent and shared
by both compilations, which use the same add/subtract/restore reduction, so the
Θ(2^n) conclusion is robust while the constant should be read as a property of
the reduction discipline rather than of modular exponentiation. The structure is
forced by three features of the standard construction acting together: flipping
the accumulator's most significant bit is the same operation as adding a power
of two and therefore commutes with the modular adder; the comparison ancilla is
coupled to that bit only through exclusive-or; and the swap network transferring
the accumulator into the multiplicand register excludes the most significant bit
entirely, so it never reaches the measured register. We confirm the diagnosis by
substituting reductions that couple the sign bit nonlinearly, which leaves the
computed function unchanged but destroys the linear structure and raises the
density from 0.473 to 0.716 — a cost increase of roughly one half arising purely
from a compilation choice, and the only such effect we observe. Reaching 24 qubits is possible only because the Walsh route computes
the exact cost without running the simulation, which stalls near 17.

Finally, we observe that the same decomposition r = β·2^α governs
matrix-product-state simulation of Shor's algorithm, where α is the number of
trailing zeros of r and the memory reduction is by a factor β² carried by the
odd part. Two structurally unrelated classical methods keying on the same
arithmetic invariant suggests a property of the algorithm rather than of either
simulator.

The practical corollary is stronger than a remark about one instance. Since every
order divides the Carmichael function λ(N), a modulus with λ(N) a power of two
lies in the free branch for *every* base — and λ(N) is a power of two exactly
when N is a power of two times a product of distinct Fermat primes. The odd
semiprimes with this property are precisely p·q with p and q both Fermat primes:
15 = 3×5, 51 = 3×17, 85 = 5×17, and so on. **N = 15 is the smallest relevant
odd semiprime, and all
seven of its usable bases lie in the free branch**, against 3 of 11 for N = 21.
The canonical demonstration instance is therefore degenerate for classical
simulation not because of an unlucky choice of base but because of the choice of
modulus, and the property that makes it the natural smallest example is the same
one that makes it uninformative as a benchmark.

---

## Honesty note that must survive into any submitted version

**Updated twice. The circuit-level invariance is now PROVED, not empirical.**

The original argument (r | 2^α ⟹ the bit function depends only on the low α bits
of e) is exact at *function* level but does **not** apply at circuit level: the
added exponent qubits remain live, because `u_a(ctrl, a^{2^i})` with
`a^{2^i} = 1` is the identity only on the valid subspace. That earlier claim is
withdrawn.

The correct mechanism is a **parity reduction**, and it is proved:

> `u_a(·,1) = A⁻¹SA` with S a product of disjoint transpositions, so S² = id and
> hence **V := u_a(·,1) is an involution**. All identity blocks apply this same V
> on distinct controls that V never modifies, so they commute and compose to
> `V^p` with `p = ⊕_{i≥α} e_i`. Averaging the Walsh character over that tail
> confines the support to `z_I ∈ {0, 1_I}` — two values regardless of tail
> length — so the support size carries no dependence on the number of exponent
> qubits.

This explains constancy and liveness together rather than in tension, predicts
the onset n_exp = α+1, predicts the branch relation h(y) = g(y ⊕ e_prev), and
fails for β>1 as required. Its sharpest consequence — support confined to
z_I ∈ {0, 1_I} — was derived before being tested and then confirmed with zero
violations across seven instances, with the two halves individually constant.

**Scope, now a criterion rather than a caveat.** The proof uses nothing about
the block beyond V² = id, the fact that V does not modify the controls, and —
this part matters — the fact that each block is controlled on **its own qubit**.
We confirm involutivity is operative by substituting synthetic blocks unrelated
to modular arithmetic: an involutive one leaves the support invariant across
four widths, while an order-3 one makes it grow by roughly a factor of two per
block, at two moduli, with a vacuity check confirming both genuinely act on the
observed bit. Multiply–swap–unmultiply qualifies because A⁻¹SA is a conjugate of
a product of disjoint transpositions.

> **Involutivity alone is not sufficient, and windowed arithmetic is where that
> shows.** Applying the criterion to windowed/table-lookup modular
> exponentiation (C36–C39) produced two constructions computing the same map
> whose identity-tail blocks are *both* involutions, exhaustively verified, and
> which sit on opposite sides of the result. Gidney-style **table lookup**
> always looks up a^(j·2^{kw}) and always multiplies; in the identity tail the
> table is all-ones, so the block never reads its window register and the tail
> exponent qubits go **dead** — a strictly sharper confinement than z_I ∈
> {0, 1_I}, and directly distinguishable from it, since the lookup support
> contains no z with any tail bit set. **Select-multiply**, which skips the
> j = 0 branch because multiplying by 1 does nothing, makes the block V gated on
> OR(window); support then grows by exactly 2^w per window and the β = 1
> advantage is lost entirely. The general condition is therefore that the
> identity-tail block's dependence on the exponent register be **affine** —
> the identity function (one fresh control) or a constant function (no control)
> both qualify, OR does not. Two further caveats belong in any statement of
> this: for the lookup construction the cost is **bounded and 2-periodic** in
> the number of tail windows rather than exactly constant, since the tail
> applies its involution unconditionally; and the analysis assumes a unitary
> unlookup, so Gidney's measurement-based uncomputation is outside it.

## Claims ledger (Paper B)

The claims ledger lives in **`CLAIMS.md`**, which is the single source of truth
for claim statuses. The claims this paper rests on are **C15, C21, C22, C18, C7,
C30, C31, C32, C27, C28, F12, C19, C20, C23, C24, C29 (narrowed), C36, C37,
C38, C39**, plus **C8** (the Walsh
identity, used as the measurement instrument) and **C17** (`perm_pps.py`, which
supplies the peak-memory numbers), whose rows sit in Paper A's section of that
file. The two earlier guesses disproved along the way — support confined to the
low exponent bits, and affineness of the identity block — are recorded in
`CLAIMS.md` under "Retracted / dead".

## Dependencies on Paper A

- The Walsh identity (Paper A, C8) is used as the measurement instrument.
- `perm_pps.py` (Paper A, C17) supplies the peak-memory numbers; it is what makes
  the n_exp = 8 peak measurement tractable (7.4 s at q = 21 versus rotation-level
  propagation failing to finish q = 17 in 20 minutes).

## Open problems

1. ~~Prove the circuit-level invariance.~~ **RESOLVED.** Proved: a parity
   reduction via an involution (C23/C24/C29). See the honesty note above.
2. **Intermediate 2-adic structure — RESOLVED, negative.** There is no
   intermediate law. The N = 323 reading (density 0.981 at t = 16 washing out
   to 1.000 by t = 24) was sampling aliasing: t = 16, 20, 24 hit residues
   4, 2, 0 mod ord₂(9) = 6, three points on an oscillation misread as a trend.
   The robust structure remains the binary β = 1 / β > 1 split. See
   `NOTES.md` §I.
3. **Does the criterion survive noise/truncation? — RESOLVED.** Peak memory
   (N_max) is unchanged at every δ tested; accuracy holds up to moderate δ but
   degrades first for the wider circuit at aggressive δ (C27/C28). See
   `NOTES.md` §T.
4. Whether the MPS correspondence extends to other methods (decision diagrams,
   tensor-network contraction orders) or is specific to these two.
5. **What is the residual non-linear structure?** Breaking the linear
   structure w = b_msb ⊕ anc (by coupling the msb nonlinearly) raises density
   to 0.716–0.721, not to 1.000 as a random function would give. Some
   non-linear structure survives the destruction of the linear one, and it is
   unidentified. See `NOTES.md` §L2, "Still open".

## Literature anchors

- Dang, Hill & Hollenberg, [arXiv:1712.07311](https://arxiv.org/abs/1712.07311) —
  MPS Shor, 60 qubits. §4 gives α = trailing zeros of r, β = r/2^α odd part;
  §5.2 gives the β² memory factor. **The 2-adic decomposition is explicit in the
  body even though the abstract says only "its factors".**
- Beauregard, [quant-ph/0205095](https://arxiv.org/abs/quant-ph/0205095) —
  Fourier modexp construction.
- Cuccaro et al., [quant-ph/0410184](https://arxiv.org/abs/quant-ph/0410184) —
  ripple-carry adder used for the Toffoli compilation.
- Gharibyan et al., [arXiv:2507.10771](https://arxiv.org/pdf/2507.10771) — the
  PPS resource framework whose extrapolation machinery this bypasses in this
  regime.

## Reproduction

```
uv run python -m experiments.experiment_c15    # controlled: fix N, vary a  (N=7 and N=21)
uv run python -m experiments.experiment_c15b   # exact constancy sweep, n_exp 2..10
uv run python -m experiments.experiment_c7     # scaling to 24 qubits, density -> 1/2
uv run python test_claims.py                   # headline rows, pinned, ~15 s
```
