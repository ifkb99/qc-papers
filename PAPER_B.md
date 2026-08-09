# The 2-Adic Structure of the Order Determines Pauli-Path Simulation Cost for Modular Exponentiation

**Ian Baker**

*Draft v1, 2026-08-08. Companion to `PAPER_A.md`, whose identity is used here as
an instrument (§3). Claim identifiers are kept out of the prose and collected in
Appendix A; the working ledger `CLAIMS.md` is supplementary material and is the
single source of truth for claim status. Numbers are post-bugfix; see §12.2.*

---

## Abstract

Shor's algorithm is the canonical target for classical simulation studies of
quantum advantage, yet its simulation difficulty is usually discussed in terms of
qubit count and circuit depth. We show that for Pauli-path simulation the
dominant cost parameter is instead an *arithmetic property of the instance*: the
2-adic structure of the multiplicative order r of the chosen base.

Writing r = β·2^α with β odd, we prove that the Pauli support carried by
propagation through a reversible modular-exponentiation circuit is **independent
of the exponent-register width when β = 1**, and measure it to be Θ(2ⁿ) with
density bounded away from zero otherwise. The invariance has an exact onset at
n_exp = α + 1. In a controlled design that fixes the modulus, ancilla layout,
gate structure and observable and varies only the base — so that r alone changes
— the support for N = 7, a = 6 (r = 2) is *exactly* 15549 terms across a 64- to
256-fold growth in Hilbert-space dimension, and peak memory *exactly* 24369
terms across a 64-fold growth, while a = 3 (r = 6) grows by a clean factor of
4.00 per two added qubits. It reproduces at N = 21 and N = 5.

The circuit-level mechanism is **proved**, not observed: every identity-tail
block applies the same permutation V = u_a(·,1), which is an involution because
it is conjugate to a product of disjoint transpositions, so the circuit depends
on the entire tail through a single parity bit. Averaging the Walsh character
over that tail confines the support to two values regardless of tail length.

The proof yields a checkable criterion rather than a family restriction — and we
report that our first statement of that criterion was **too strong**, corrected
here. Applying it to Gidney-style windowed arithmetic shows the invariance
survives, in a strictly sharper form, but that an obvious optimisation destroys
it entirely.

The practical corollary is uncomfortable for benchmarking: since every order
divides λ(N), a modulus with λ(N) a power of two lies in the free branch for
*every* base. **N = 15 is the smallest such modulus and all seven of its usable
bases are free**, against 3 of 11 for N = 21 — so the canonical demonstration
instance is degenerate for classical simulation because of the modulus, not an
unlucky base.

---

## 1. Introduction

The classical simulation cost of Shor's algorithm is usually parameterised by
size: how many qubits, how deep, how many T gates. For Pauli-path simulation
(PPS, also sparse Pauli dynamics) we find that size is the wrong parameter over a
large part of the instance space. Two instances with identical qubit counts,
identical gate structure, identical ancilla layout and identical observable can
differ by an unbounded factor in simulation cost, determined entirely by a
number-theoretic property of the base.

Specifically: let r be the multiplicative order of a mod N and write r = β·2^α
with β odd. Then

- **β = 1** (r a power of two): the exponent register is *free*. Adding exponent
  qubits — which is exactly what buys precision in the continued-fractions
  post-processing, and which one would expect to be the expensive resource —
  does not increase Pauli support at all, beyond a threshold we locate exactly.
- **β > 1**: cost is Θ(2ⁿ), quadrupling per two added exponent qubits, and
  asymptotically no better than state-vector simulation.

This is not a small effect at the margin; it is a dichotomy, and it is invisible
to a resource model that extrapolates from circuit size.

**What is proved and what is measured.** The dichotomy at the level of the
idealised bit function is elementary (§4). The circuit-level statement is the
substantive one, and its mechanism is proved (§7) rather than inferred from the
measurements (§5, §6). We are explicit throughout about which is which, and §12
records five claims withdrawn or narrowed along the way.

---

## 2. Setting

We use two independent, end-to-end verified reversible constructions of
|x⟩|1⟩ ↦ |x⟩|a^x mod N⟩:

- a **Toffoli compilation** — Cuccaro et al. ripple-carry addition taken mod 2^m,
  with the standard add / subtract-N / conditional-restore modular reduction;
- a **Fourier compilation** — Beauregard's 2n+3 qubit construction, with Draper
  phase-space addition.

Both are gated by layer-by-layer correctness suites and cross-checked against
each other on the valid subspace. The observable throughout is a
computational-basis Z on a single bit of the multiplicand register; §5 varies it
to confirm observable-independence.

The **controlled design** is the methodological core and is worth stating
explicitly, because an earlier version of this work was invalidated by failing
it: we fix N, the ancilla layout, the gate structure and the observable, and vary
**only the base a**. Then r is the single quantity that changes. Comparisons that
vary width and instance simultaneously are not admissible here.

---

## 3. The instrument (Paper A)

We use as a lemma the identity established in the companion paper.

**Lemma 1 (Paper A, Theorem 1).** *For a circuit implementing a basis
permutation π and a computational-basis observable Z_j, the pullback π†Z_jπ is
the diagonal
operator (−1)^{g(y)} with g(y) = bit j of π(y), its Pauli expansion is purely
Z-type, and its coefficient vector is exactly the Walsh–Hadamard transform of
(−1)^g. The number of Pauli terms PPS must carry is the Walsh sparsity of g.*

Two consequences make this paper possible.

1. **Cost is computable without simulating.** We reach 30 qubits by computing
   spectra directly, where rotation-level propagation stalls near 17.
2. **Peak memory is measurable at scale.** Propagating X, CNOT and Toffoli as
   atomic permutations (Paper A, §5) keeps the expansion Z-type throughout and
   makes the n_exp = 8 peak measurement tractable — 7.4 s at 21 qubits, against
   rotation-level propagation failing to finish 17 qubits in 20 minutes.

Because the Walsh basis is the character group of (Z/2)^t, and because
a^e mod N depends only on e mod r, the arithmetic of r *relative to powers of
two* controls the spectrum directly. That is the whole idea; the rest is
establishing it where it actually matters, at circuit level.

---

## 4. The dichotomy at function level

Consider the idealised bit function g(e) = bit_j(a^e mod N) on t exponent bits.

If β = 1 then r = 2^α divides 2^t for every t ≥ α, so g depends only on the low
α bits of e. Its Walsh spectrum is supported on those α coordinates and has
**constant size independent of t**. If β > 1 the period is incommensurate with
the Walsh basis and no such collapse occurs.

Measured: at t = 24, N = 15 with a = 7 (r = 4) has a spectrum of exactly **4
nonzero coefficients**; sparsity is constant at 4 across t = 11…22.

**Density for β > 1 need not be 1**, and this is worth stating because it is easy
to over-claim. It oscillates with period ord₂(β): N = 7, a = 3 (r = 6) sits at
exactly 0.500000 on even t including t = 24, while N = 21, a = 2 reaches
1.000000 at t = 24.

**There is no intermediate law — a negative result.** We looked for a
quantitative law in α interpolating the two branches and there is none. An
apparent one at N = 323 (density 0.981 at t = 16 washing to 1.000 by t = 24) was
**sampling aliasing**: t = 16, 20, 24 hit residues 4, 2, 0 mod ord₂(9) = 6 — three
points of an oscillation misread as a trend. The robust structure is the binary
β = 1 / β > 1 split, and the claim should be stated without hedging about
intermediate regimes.

**Caveat, checked.** The idealised function is not automatically a proxy for the
real one. Rechecked with the actual table h[c] = bit_j(a^c mod N): the constancy
holds, and α is confirmed *not* to organise the data (α = 0 at r = 3 matches
α = 1 at r = 6, while α = 1 spans r = 6, 10, 18, 22 with wildly different
behaviour). Separately, the real modexp table is measurably **non-generic** —
sparser than random at r = 6, denser at r = 10 — so random-table densities must
not be used as a proxy for real ones.

---

## 5. The dichotomy at circuit level

Function-level behaviour need not survive compilation into a circuit with
ancillas and scratch registers, and in this case the naive reason it might is
**false** (§12.2). We therefore measure it directly.

**Final support**, fixing N = 7 and varying only a:

| a | r | support vs n_exp |
|---|---|---|
| 6 | 2 (β=1) | **exactly 15549**, invariant over a 64–256× growth in dimension |
| 3 | 6 (β=3) | grows by a factor **4.00 per two added qubits** |

**Peak memory**, the quantity that actually bounds a run:

| a | r | peak vs n_exp = 2, 4, 6, 8 |
|---|---|---|
| 6 | 2 | **exactly 24369** at every width (64× dimension growth) |
| 3 | 6 | grows 4.02× per step |

Reproduced at **N = 21 and N = 5**, giving three moduli. The invariance is
observable-independent.

So the exponent register — the resource that buys post-processing precision — is
free for PPS when r is a power of two, and quadruples cost per two qubits
otherwise.

---

## 6. The onset is exact

The invariance does not hold from the smallest width; it switches on.

The controlled-multiplier block attached to exponent bit i multiplies by
a^(2^i) mod N, which is the identity on the valid subspace exactly when i ≥ α.
The first such block therefore appears at **n_exp = α + 1**, and the support
locks there and not before.

**Confirmed for α = 1, 2, 3 and 4.** The α = 2 case is the first discriminating
one: N = 5 with a = 2 (r = 4) **grows once**, 15493 → 32143, and then freezes —
exactly one step of growth before the lock. α = 3 and α = 4 were for a long
time only "still growing at the largest width we can reach, as predicted"; they
are now measured:

| N | a | r | α | \|support\| by n_exp = 1, 2, … | locks at |
|---|---|---|---|---|---|
| 17 | 2 | 8 | 3 | 255,104 · 1,037,405 · 2,093,137 · **4,188,525 · 4,188,525** | 4 = α+1 |
| 41 | 3 | 8 | 3 | 2,070,878 · 8,346,567 · 16,766,478 · **33,539,711 · 33,539,711** | 4 = α+1 |
| 17 | 3 | 16 | 4 | 255,356 · … · 4,188,537 · **8,379,626 · 8,379,626** | 5 = α+1 |
| 41 | 6 | 40 | 3 | 2,072,174 · 8,346,761 · 16,766,484 · 33,539,777 · 67,086,624 | never (β=5) |

Growth is strict at every one of the seven steps below the onset, so the lock
is not an artifact of a flat measurement. The last row is the matched control:
**same modulus, same α, same circuit width, only β differs.** N = 41 admits
both r = 8 and r = 40 = 5·2³, which makes the comparison exact rather than
merely careful.

> **A caution for anyone benchmarking.** At n_exp = 4 the β = 1 and β = 5 rows
> of that table differ by 66 terms in 3.4 × 10⁷ — 2 parts per million. **The
> invariant is invisible in a cost measurement at any single width.** Only the
> growth separates the two regimes, which is why every claim here is stated
> over a swept n_exp and not at a fixed size.

> **Why this was nearly missed, and what it cost.** Every earlier sweep used
> α = 1 and started at n_exp = 2 = α + 1 — exactly on the threshold. The
> invariance looked unconditional purely by luck of parameter choice. The rule
> was found only when a third modulus introduced an α = 2 instance. We record
> this because it is the strongest argument in this paper for stepping a sweep
> by one and starting below the expected threshold.

The locked value sits close to half the Hilbert-space dimension at the lock
point, so the support saturates to half density and then freezes in absolute
terms while density falls fourfold per two added qubits.

---

## 7. The mechanism, proved

### 7.1 The parity reduction

The natural explanation — that the added exponent qubits are inert because
a^(2^i) = 1 — is **false at circuit level**, and we withdrew it (§12.2). The
blocks are the identity only on the *valid subspace*; off it they act, and the
added qubits remain live. Constancy and liveness looked to be in tension.

They are not. The resolution:

**Theorem 2.** *Write the identity-tail block as V := u_a(·, 1). Then
V = A⁻¹SA, where A is the multiply-accumulate and S the controlled swap layer.
S is a product of **disjoint transpositions**, so S² = id, and therefore
V² = id: V is an involution.*

*Every identity block applies this same V, controlled on its own qubit, which V
does not modify. The blocks therefore commute and compose to V^p with*

> p = ⊕_{i ≥ α} e_i,

*a single parity bit. The circuit depends on the entire identity tail through one
effective variable, however many qubits that tail spans.*

**Corollary 3.** *Averaging the Walsh character over the tail confines the
support to z_I ∈ {0, 1_I} — two values regardless of tail length — so the support
size carries no dependence on the number of exponent qubits.*

This explains constancy and liveness *together*, predicts the onset n_exp = α+1
independently, predicts the branch relation h(y) = g(y ⊕ e_prev), and fails for
β > 1 as required.

**Evidence, each link checked.** V² = id verified exhaustively at 2¹⁵ and 2²¹.
Branch magnitudes |ĥ(z)| = |ĝ(z)| match bit-for-bit for β = 1 (error 0.00e+00)
against 14934/1036508 mismatches for β > 1. The sign pattern is a linear
character, exact on 100% of the support. Corollary 3's confinement was **derived
before being tested** and then confirmed with **zero violations across 7/7 β = 1
instances**, with the two halves individually constant (7770/7779), against 48189
and 1556046 violations for the β > 1 controls.

> **A caveat on how that confinement is tested, which cost us a vacuous
> measurement.** At |I| = 1 the statement "z_I ∈ {0, 1_I}" is a tautology: a
> single bit is either all-zeros or all-ones. Rows at the onset width therefore
> carry no evidence, and a β > 1 control evaluated there *passes*. Only widths
> with |I| ≥ 2 test anything. All counts quoted here are from such widths.

**Corollary 3 is in fact stronger than a statement about size.** Step
(iv) computes the surviving coefficient explicitly and the expression contains
no |I|, so the supports at consecutive widths must coincide as *sets*, not
merely in cardinality. Writing each element as (z_rest, tailflag) with the flag
recording whether z_I = 0 or 1_I, the sets are bit-for-bit identical at α = 3
(two moduli) and α = 4, over supports of 4.2, 8.4 and 33.5 million elements.
Predicted before measuring. The two halves are individually constant but
unequal to each other (2,094,285 against 2,094,240): the invariance is per
half, not a symmetry between halves.

**Two failed routes, recorded so they are not retried:** affineness of the
identity block (refuted, 14336 of 32768 violations) and the (z, z⊕e) pairing (an
n_exp = 1 artifact, not r-dependent).

### 7.2 The criterion, and a correction to it

Steps of the proof use nothing about V beyond V² = id and the fact that V does
not modify its controls. So the scope caveat becomes a *checkable criterion*
rather than a family membership. We verified involutivity is operative by
substituting synthetic blocks unrelated to modular arithmetic: an involutive one
(controlled swap, order 2) leaves the support constant across four widths at two
moduli, while an order-3 block grows by roughly 2× per block, with a vacuity
check confirming both genuinely act on the observed bit.

> **We first stated this criterion as "V² = id is the entire condition". That is
> too strong, and we correct it here.** Every synthetic block we
> tested was controlled on exactly one fresh qubit, so the experiment could not
> separate *"V is an involution"* from *"the control is a single qubit"*. The
> second hypothesis is load-bearing. The correct statement is that the
> identity-tail block's dependence on the exponent register must be **affine**:
> the identity function (one fresh control qubit) and the constant function (no
> control at all) both qualify; a nonlinear control such as OR does not.

### 7.3 Windowed arithmetic

The natural test of the criterion is windowed / table-lookup arithmetic, which is
what is actually proposed for hardware. It is also where the correction above
shows its teeth: two windowed constructions computing the same map have
identity-tail blocks that are **both** involutions, verified exhaustively, and
they sit on opposite sides of the result.

- **Table lookup (Gidney).** Always look up T[j] = a^(j·2^{kw}) mod N and always
  multiply. In the identity tail every table entry is 1, so the lookup
  permutation does not read the window register at all and the tail exponent
  qubits go **dead**. This is a *strictly sharper* confinement than Corollary 3
  and is directly distinguishable from it: the lookup support contains **no** z
  with any tail bit set, whereas Corollary 3's support contains the all-ones tail
  vector. Measured at the same modulus and base: standard construction contains
  1_I, lookup contains no tail bit at all.
- **Select-multiply.** Skip the j = 0 branch, because multiplying by 1 does
  nothing. The block becomes V gated on OR(window) — nonlinear — and support
  grows by exactly **2^w per window** (measured ×4.00 at w = 2), with density
  pinned at ½, i.e. the generic β > 1 behaviour. The β = 1 advantage is lost
  entirely: 2061438 terms against 32075 at four windows, a 64× gap growing as 4^K.

**This is validated against the source, not reconstructed.** Gidney states the
mechanism as his own design rationale: *"this also removes the need for the
multiplications to be controlled, because the table lookup can evaluate to the
factor 1 in cases where none of the exponent qubits are set"* (arXiv:1905.07682
§3.5). A choice made to reduce Toffoli count is exactly what preserves the PPS
advantage — designed by nobody for that purpose. His relabelling swap, moreover,
is what makes the block an involution at all, and his `if a is not target:
swap(a, b)` emits a physical swap exactly when the block count is odd.

**Two caveats belong in any statement of this.** For the lookup construction the
cost is **bounded and 2-periodic** in the number of tail windows rather than
exactly constant, because the tail applies its involution unconditionally rather
than under a live control. And the analysis assumes a **unitary** unlookup;
Gidney's measurement-based uncomputation is genuinely outside it, and that is the
one real remaining gap.

---

## 8. Robustness to truncation

An exact-arithmetic invariance is of limited interest if it evaporates under the
approximations a real run makes. It largely does not.

- **Terminal thresholding preserves it at every δ.** Counts are identical across
  n_exp = 3, 4, 5 at all seven δ values tested, at both β = 1 moduli, while the
  β = 3 control diverges. This follows from Corollary 3: identical magnitude
  multisets mean any magnitude threshold keeps identical counts.
- **Peak cost survives incremental truncation** — the truncation PPS actually
  performs. N_max is the *same* at every δ tested (24369 / 24369 / 1028 / 34).
  That is the quantity bounding memory.
- **Accuracy has a limit, and it must be stated.** The final term count
  drifts at aggressive δ (8 → 16 → 32) and ⟨O⟩ collapses to 0 at δ = 10⁻¹ for the
  wider circuits while n_exp = 3 stays exact. The wider circuit is *more fragile
  at the same δ* despite an identical exact spectrum, because it has more gates
  and therefore more incremental truncation events.

The honest practical claim: **peak cost is free in the exponent register at every
δ; accuracy is free up to moderate δ.** It must not be stated unqualified.

---

## 9. The generic case, and what sets its constant

For instances with an odd factor — the generic and cryptographically relevant
case — Walsh density converges monotonically to one half: 0.473 → 0.4994 over 15
to **30** qubits, growth **1.006 bits per qubit** (two series varying the
modulus at fixed exponent width). PPS cost is therefore Θ(2ⁿ), asymptotically no
better than state-vector simulation.

**The limiting constant of ½ is neither accidental nor algorithmic.**
The support carries a linear structure w = b_msb ⊕ anc, pairing the accumulator's
sign bit with the modular-reduction comparison ancilla, which confines it to a
hyperplane. That a linear structure confines the Walsh support to a coset is
classical — Carlet, *Boolean Functions for Cryptography and Coding Theory*,
Proposition 29 — and we **cite rather than claim** it; our contribution is
locating it here and reading it as a cost cap. It is observable-independent and
shared by both compilations, which use the same reduction discipline. So the
Θ(2ⁿ) conclusion is robust while **the constant should be read as a property of
the reduction, not of modular exponentiation**.

It is forced by three features acting together: flipping the accumulator's msb
*is* adding a power of two and so commutes with the modular adder (verified
exhaustively, 0/256 violations); the comparison ancilla is coupled to that bit
only by XOR; and the swap network excludes the msb entirely, so it never reaches
the measured register.

**Confirmed by substitution.** Replacing the reduction with one that
couples the sign bit *nonlinearly* leaves the computed function unchanged but
destroys the linear structure and raises density from 0.473 to 0.716 — a cost
increase of roughly one half from a pure compilation choice, and the only such
effect we observe. Linear couplings (CNOT, cswap) leave it intact. The framing is
the reverse of the usual one: the standard construction is already cheaper than
it needs to be, for a reason nobody designed.

---

## 10. Two corollaries about instances

### 10.1 N = 15 is degenerate by construction

Every order divides the Carmichael function λ(N). So if λ(N) is a power of two,
*every* base lies in the free branch — and λ(N) is a power of two exactly when N
is a power of two times a product of distinct Fermat primes. For odd semiprimes
this means p·q with both p and q Fermat: 15 = 3×5, 51 = 3×17, 85 = 5×17, …

**N = 15 is the smallest, and all seven of its usable bases are free**, against
3 of 11 for N = 21.

The canonical demonstration instance is therefore degenerate for classical
simulation not because of an unlucky base but because of the *modulus* — and the
property that makes 15 the natural smallest example is the same one that makes it
uninformative as a benchmark. Simulation results on N = 15 should not be
extrapolated.

### 10.2 Most of the support is inert for the real input state

Shor's actual initial state puts the exponent register in |+⟩. Since ⟨Z⟩ = 0
there, *every* Pauli term with support on the exponent register contributes
exactly zero to the expectation value. Combined with Corollary 3 this predicts a
useful fraction of exactly **2^−(α+1)**. Measured: **75.0%** of the support is
inert at α = 1 and **87.5%** at α = 2, across widths, with the jump occurring
exactly at n_exp = α + 1 — reproducing §6's onset from an independent direction.

**This is not a speedup and we do not present it as one.** PPS meets the input
state only at the end of propagation, and exponent support is not monotone under
back-propagation (a term carrying it can lose it and go on to contribute), so
these terms cannot simply be pruned early. What the observation does is separate
*cost* from *useful work*: at α = 2, seven of every eight Pauli terms carried are
irrelevant to the answer.

---

## 11. Cross-method corroboration

The same decomposition r = β·2^α governs matrix-product-state simulation of
Shor's algorithm. Dang, Hill & Hollenberg (arXiv:1712.07311) give α as the number
of trailing zeros of r and a memory reduction carried by β²; the 2-adic
decomposition is explicit in their §4 and §5.2 even though their abstract says
only "its factors".

Two structurally unrelated classical methods — Heisenberg-picture Pauli
propagation and Schrödinger-picture tensor networks — keying on the same
arithmetic invariant is suggestive of a property of the algorithm rather than of
either simulator. We state this as corroboration, not proof; a third method
would make it hard to argue with, and remains open (§13).

---

## 12. Limitations and the honesty record

### 12.1 Scope

- **Constructions.** The proof of Theorem 2 uses the multiply–swap–unmultiply
  form, so it covers both compilations here and Vedral/Beauregard-style
  constructions, plus anything satisfying the corrected criterion of §7.2. It
  does not transfer automatically to a modexp built otherwise.
- **Observables.** Computational-basis only. Real Shor measures after an inverse
  QFT, which leaves the diagonal and is outside both this paper and Paper A.
- **Sizes.** 30 qubits via the Walsh route, n_exp = 8 for peak memory. The α = 3
  onset is confirmed at two moduli and the α = 4 onset at one (§6); α ≥ 5 needs
  r = 32, whose smallest instance is out of reach.
- **No implication for factoring.** Efficient classical simulation of Shor's
  algorithm on general inputs would be a classical factoring algorithm. Nothing
  here bears on that; these are diagnostic results about where PPS breaks, and
  the free branch is exactly the branch where the order is already easy to find.

### 12.2 Retractions and corrections

1. **"The added exponent qubits are inert because a^(2^i) = 1."** Withdrawn.
   True at function level, false at circuit level — the blocks are the identity
   only on the valid subspace. Replaced by the parity reduction (§7.1).
2. **"Support is confined to the low exponent bits."** Disproved: the added
   qubits are live in about half the support terms. Superseded by Corollary 3.
3. **"The (z, z⊕e) pairing is the mechanism."** Disproved: an n_exp = 1 artifact.
4. **"V² = id is the entire condition".** Narrowed, §7.2. Not wrong about
   anything it tested; wrong about what it had tested — every block tested had
   exactly one control, so the experiment could not separate the two hypotheses.
5. **An intermediate 2-adic law in α.** There is none; the apparent case was
   sampling aliasing (§4).

Two methodology errors are recorded because they produced confident-looking
numbers first: degenerate random tables at small r, and drawing a fresh table per
width, which confounded width-dependence with table variance — a direct violation
of "vary exactly one parameter".

### 12.3 Reproducibility

```
uv run python -m experiments.experiment_c15    # controlled: fix N, vary a (N=7, N=21)
uv run python -m experiments.experiment_c15b   # exact constancy sweep, n_exp 2..10
uv run python -m experiments.experiment_c7     # scaling to 24 qubits, density -> 1/2
LAB_GPU=1 uv run python -m experiments.experiment_c21_onset  # the onset at alpha = 3, 4
LAB_GPU=1 uv run python -m experiments.experiment_c7_scale   # the same series to 30 qubits
uv run python -m experiments.experiment_windowed   # the windowed criterion of §7.3
uv run python test_claims.py                   # headline rows, pinned, ~15 s
```

Every claim resolves, via the claim map of Appendix A, to a row in `CLAIMS.md`.
Headline rows are re-verified by an executable regression suite
(`test_claims.py`) inside a nine-suite correctness gate. Experiments declare
predictions *before* measurement and carry must-fail controls; the harness
reports when a control fails to fail, which has caught two vacuous results that
had already produced confident-looking numbers — one of them a confinement test
in this paper, discussed in §7.1.

---

## 13. Open problems

1. **A third simulation method on the same invariant.** Decision diagrams (MQT
   DDSIM) are the natural candidate: rerun the controlled design and see whether
   memory keys on β. Mostly integration work, but it would move §11 from
   suggestive to hard to argue with.
2. **Through the inverse QFT.** Everything here stops where the identity stops.
   Real Shor measures the exponent register after an inverse QFT; the Walsh basis
   is the character group of (Z/2)^t, the QFT diagonalises translation on Z/2^t,
   and the cost story is already 2-adic. Whether there is an exact
   characterisation in a mixed character basis is the natural next question, and
   it may simply be dense and structureless.
3. **Measurement-based uncomputation.** Gidney's unlookup is not unitary, so Lemma 1
   does not apply to it as written (§7.3). This is the one real gap in the
   windowed analysis.
4. **Whether the MPS correspondence extends** to other contraction orders and
   tensor-network methods, or is specific to these two.

---

## 14. Conclusion

For Pauli-path simulation of modular exponentiation, the dominant cost parameter
is not size but arithmetic. Writing r = β·2^α, the exponent register is free when
β = 1 — provably, with an exact onset at n_exp = α+1 and a proved mechanism — and
Θ(2ⁿ) otherwise. The free branch is reachable for *every* base of certain
moduli, which makes the field's canonical benchmark instance degenerate for
reasons that have nothing to do with the base chosen.

The criterion we extract is checkable rather than a family restriction, and
applying it to the arithmetic people actually propose to run reveals both that
the invariance survives — in a sharper form than we proved — and that a single
obvious optimisation would destroy it.

---

## Appendix A — claim map

The working ledger `CLAIMS.md` records every claim's status, evidence and
location, and is supplied as supplementary material. Inline identifiers are
kept out of the prose; this table is the mapping. Identifiers prefixed
*Paper A* have their rows in the companion paper's section of the ledger.

| section | claims |
|---|---|
| 3. The instrument | Paper A: C8, C17 |
| 4. The dichotomy at function level | F12, C33 |
| 5. The dichotomy at circuit level | C15, C18 |
| 6. The onset is exact | C21 |
| 7.1 The parity reduction | C23, C24, C43 |
| 7.2 The criterion, and a correction to it | C29 (narrowed) → C36 |
| 7.3 Windowed arithmetic | C36, C37, C38, C39 |
| 8. Robustness to truncation | C27, C28 |
| 9. The generic case, and what sets its constant | C7, C30, C31, C32 |
| 10.1 N = 15 is degenerate by construction | C20, C22 |
| 10.2 Most of the support is inert for the real input state | C42 |
| 11. Cross-method corroboration | C19 |
| 12.2 Retractions and corrections | C29, and the retracted rows of the ledger |
