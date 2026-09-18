# The 2-Adic Structure of the Order Determines Pauli-Path Simulation Cost for Modular Exponentiation

**Ian Baker**

ifkb99@gmail.com

*Working consolidation, 2026-09-10. This file is the Paper B manuscript and
tracking document: the arithmetic mechanism in the full representation, how
state-aware contraction changes its consequences, and the conditional-sampling
frontier. New synthesis is in §10.2–10.3 and §11.1–11.4. Detailed evidence and
open tasks remain in the research ledger. Companion to Paper A (*Walsh–Hadamard Sparsity Exactly
Determines Pauli-Path Simulation Cost for Reversible Quantum Arithmetic*), whose
identity is used here as an instrument (§3). Claim identifiers are kept out of
the prose and collected in Appendix A; the working ledger `CLAIMS.md` is
supplementary material and is the single source of truth for claim status.
Numbers are post-bugfix; see §12.2.*

---

## Abstract

Shor's algorithm is the canonical target for classical simulation studies of
quantum advantage, yet its simulation difficulty is usually discussed in terms of
qubit count and circuit depth. For the full, uncontracted diagonal Pauli
representation of the studied modular-exponentiation construction, an
*arithmetic property of the instance* controls a striking width-scaling split:
the 2-adic structure of the multiplicative order r of the chosen base.

Writing r = β·2^α with β odd, we prove that the Pauli support carried by
propagation through the studied reversible modular-exponentiation construction is
**independent of the exponent-register width when β = 1**, and measure its
density to be bounded away from zero in the sampled β > 1 circuit families —
0.473 → 0.4994, monotone, up to 30 qubits — with growth consistent with Θ(2ⁿ).
The parity proof guarantees the invariant once the identity tail is present;
the first-lock location is measured for the observables and instances reported.
In a controlled design that fixes the modulus, ancilla layout, arithmetic
template and observable and varies the base — while constants and emitted gates
necessarily vary with the base, so r is the intended arithmetic parameter rather
than literally the only circuit change — the support for N = 7, a = 6 (r = 2) is *exactly* 15549 terms across a 64- to
256-fold growth in Hilbert-space dimension, and peak memory *exactly* 24369
terms across a 64-fold growth, while a = 3 (r = 6) grows by a clean factor of
4.00 per two added qubits. It reproduces at N = 21 and N = 5.

These costs describe the full, uncontracted Heisenberg representation, not a
lower bound for state-aware simulation. For the pre-QFT work-register observable,
exact contraction of finished exponent controls bounds retained support
independently of exponent width in both branches (§10.2).
For independent control windows acting through one involution, the reduced
tail depends only on activation bias, with an exact residual formula (§10.3).
Separately, validated semiclassical instruments sample actual post-QFT outputs.
On-demand replay removes full-scratch tables but remains orbit-sized; charged
order-informed classical baselines are stronger on the tested instances
(§11.2–11.4). These are not new general simulation methods.

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
*every* base. **N = 15 is the smallest relevant odd semiprime with this property
and all seven of its usable
bases are in this width-invariant full-operator branch**, against 3 of 11 for
N = 21. The canonical instance therefore fails to probe the odd-order-factor
behavior; this is a limitation of that benchmark, not a lower bound on other
instances or other simulators.

---

## 1. Introduction

The classical simulation cost of Shor's algorithm is usually parameterised by
size: how many qubits, how deep, how many T gates. For Pauli-path simulation
(PPS, also sparse Pauli dynamics) we find that size is the wrong parameter over a
large part of the instance space. Two instances with identical qubit counts,
the same arithmetic template, identical ancilla layout and identical observable can
differ in simulation cost by a factor that grows without bound in every sweep
we can run, tracking a number-theoretic property of the base.

Specifically: let r be the multiplicative order of a mod N and write r = β·2^α
with β odd. Then

- **β = 1** (r a power of two): the exponent register is *free*. Adding exponent
  qubits — which is exactly what buys precision in the continued-fractions
  post-processing, and which one would expect to be the expensive resource —
  does not increase full-operator Pauli support beyond the proved sufficient
  threshold; the first-lock widths are measured in §6.
- **β > 1**: in the circuit families measured here, cost quadruples per two
  added exponent qubits, with density bounded away from zero at every size
  reached — growth consistent with Θ(2ⁿ). This is an empirical branch result,
  not a theorem for every selected output bit or compilation.

This is not a small effect at the margin; it is a dichotomy, and it is invisible
to a resource model that extrapolates from circuit size.

**What is proved and what is measured.** The β=1 statement at the level of the
idealised bit function is elementary; scalar-period counterexamples rule out
a universal β>1 density theorem (§4). The circuit-level β=1 mechanism is proved
(§7), while β>1 growth is measured (§5, §6). We are explicit about which is which, and §12
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
computational-basis Z on a single bit of the multiplicand register (the register
that ends holding a^x mod N); §5 varies it to confirm observable-independence.

The **controlled design** is the methodological core and is worth stating
explicitly, because an earlier version of this work was invalidated by failing
it: we fix N, the ancilla layout, the arithmetic template and the observable,
and vary **only the base a**. The loaded constants and emitted gate sequence can
also change with a; r is the intended arithmetic control, not literally the
only circuit-level change. Comparisons that vary width and instance
simultaneously are not admissible here.

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
**constant size independent of t**. If β > 1 there is no corresponding
universal conclusion for every selected bit: the relevant period is the minimal
period of that scalar output function, which can be a proper divisor of r. For
example, N=13, a=4 has r=6 but its least-significant-bit sequence on the orbit
is 1,0,1,0,1,0, so g(e)=1⊕e₀ and has Walsh sparsity 1. The β>1 circuit
branch below is therefore explicitly empirical and observable-dependent.

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

| n_exp | qubits | dimension | a = 6 (r = 2, β = 1) | a = 3 (r = 6, β = 3) | growth (a = 3) |
|---|---|---|---|---|---|
| 2 | 15 | 32,768 | **15,549** | 15,539 | — |
| 4 | 17 | 131,072 | **15,549** | 64,353 | 4.14× |
| 6 | 19 | 524,288 | **15,549** | 258,691 | 4.02× |
| 8 | 21 | 2,097,152 | **15,549** | 1,037,728 | 4.01× |
| 10 | 23 | 8,388,608 | **15,549** | 4,155,634 | 4.00× |

**Peak memory**, the quantity that actually bounds a run
(permutation-native propagation, Paper A §5):

| n_exp | qubits | a = 6 (r = 2) | a = 3 (r = 6) |
|---|---|---|---|
| 2 | 15 | **24,369** | 24,412 |
| 4 | 17 | **24,369** | 98,018 (4.02×) |
| 6 | 19 | **24,369** | — |
| 8 | 21 | **24,369** | — |

(The peak control was run at two widths; the final-support control at five.)

Reproduced at **N = 21 and N = 5**, giving three moduli — at N = 21, a = 8
(r = 2) gives 1,037,174 at both n_exp = 2 and 4 while the dimension quadruples,
against a = 2 (r = 6) growing 1,037,322 → 4,186,980 and a = 4 (r = 3) growing
512,784 → 4,152,181. The invariance is observable-independent.

So the exponent register — the resource that buys post-processing precision — is
free for PPS when r is a power of two once the identity tail is present. In the
reported β>1 circuit controls, cost quadruples per two added qubits; that growth
is empirical rather than a universal odd-factor theorem.

---

## 6. The proved threshold and measured first-lock rows

In the reported rows the invariance does not hold from the smallest width; it
switches on at the measured first-lock width.

For β=1, the controlled-multiplier block attached to exponent bit i multiplies by
a^(2^i) mod N, which is the identity on the valid subspace exactly when i ≥ α.
The parity proof therefore guarantees width-invariance once the first such block
is present, at **n_exp ≥ α + 1**. Whether an earlier width happens to have
already stabilised is an additional observable- and construction-dependent
question; the reported rows measure first locking at α + 1, but do not by
themselves prove a universal lower bound.

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

Growth is strict at every one of the seven measured steps below the reported
onset, so those locks are not artifacts of a flat measurement. The last row is
the matched control:
**same modulus, same α and same circuit width, with the base varied.** N = 41
admits both r = 8 and r = 40 = 5·2³. The loaded constants and emitted gates
also vary with the base, as noted in §2.

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

This explains constancy and liveness *together*, guarantees the invariant for
n_exp ≥ α+1, predicts the branch relation h(y) = g(y ⊕ e_prev), and fails for
β > 1 in the tested controls as required. The measured first-lock rows are
reported separately in §6.

Full proofs of Theorem 2 and Corollary 3 are given in Appendix B.

**Evidence, each link checked.** V² = id verified exhaustively at 2¹⁵ and 2²¹.
Branch magnitudes |ĥ(z)| = |ĝ(z)| match bit-for-bit for β = 1 (error 0.00e+00)
against 14934/1036508 mismatches for β > 1. The sign pattern is a linear
character, exact on 100% of the support. Corollary 3's confinement was **derived
before being tested** and then confirmed with **zero violations across 7/7 β = 1
instances**, with the two halves individually constant (7770 terms with z_I = 0
and 7779 with z_I = 1_I, summing to the invariant 15549), against 48189 and
1556046 violations for the β > 1 controls.

> **A caveat on how that confinement is tested, which cost us a vacuous
> measurement.** At |I| = 1 the statement "z_I ∈ {0, 1_I}" is a tautology: a
> single bit is either all-zeros or all-ones. Rows at the onset width therefore
> carry no evidence, and a β > 1 control evaluated there *passes*. Only widths
> with |I| ≥ 2 test anything. All counts quoted here are from such widths.

**Corollary 3 is in fact stronger than a statement about size.** Its proof
computes the surviving coefficient explicitly (Appendix B), and the expression
contains no |I|, so the supports at consecutive widths must coincide as *sets*, not
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
does real work: two windowed constructions computing the same logical modular
map on the valid subspace have
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
  nothing. The block becomes V gated on OR(window) — nonlinear — and the
  reported experiment grows by **2^w per window** (×4.00 at w = 2), with density
  pinned at ½ in that construction. This is an empirical windowed result, not a
  general consequence of nonlinearity. The β = 1 advantage is lost entirely in
  the measured case: 2061438 terms against 32075 at four windows, a gap that
  grows fourfold per added window and stands at 64× there.

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
than under a live control. And the unlookup may be unitary or measurement-based.
For a computational-basis observable, measuring the lookup register in the X basis,
applying the diagonal phase fixup and returning the register to |0⟩ act in the
Heisenberg picture as the pullback by the map that clears that register. Propagation
stays diagonal and branch-free, and Lemma 1 holds with the basis permutation replaced
by this non-injective population map. The two circuits compute the same full-space
function on every input whose lookup register starts clean, a set that contains the
valid subspace, and the dead tail holds for both. Off that set the two functions
need not agree, and do not at an instance traced by hand; the sizes reported in this section are for the unitary unlookup, and the
measurement-based sizes are not reported. The activation
ancillae are still uncomputed unitarily here.

---

## 8. Robustness to truncation

An exact-arithmetic invariance is of limited interest if it evaporates under the
approximations a real run makes. It largely does not.

- **Terminal thresholding preserves it at every δ.** Counts are identical across
  n_exp = 3, 4, 5 at all seven δ values tested, at both β = 1 moduli, while the
  β = 3 control diverges. This follows from Corollary 3: identical magnitude
  multisets mean any magnitude threshold keeps identical counts.
- **Peak cost was unchanged in the tested incremental truncation runs** — the
  truncation PPS actually performs. At each tested δ, N_max was identical across
  n_exp = 3, 4 and 5; the four δ levels gave N_max = 24369 (δ = 0), 24369
  (10⁻⁴), 1028 (10⁻²) and 34 (10⁻¹). This is the quantity bounding memory;
  equality at all thresholds is an empirical result, not a theorem at every δ.
- **Accuracy has a limit, and it must be stated.** The final term count
  drifts at aggressive δ (8 → 16 → 32) and ⟨O⟩ collapses to 0 at δ = 10⁻¹ for the
  wider circuits while n_exp = 3 stays exact. The wider circuit is *more fragile
  at the same δ* despite an identical exact spectrum, because it has more gates
  and therefore more incremental truncation events.

The honest practical claim: **exact peak cost is free in the exponent register;
the tested truncation runs preserve that equality, and accuracy is free up to
moderate δ in those runs.** The incremental peak and accuracy observations
should not be generalized to every threshold or circuit. Terminal thresholded
count invariance, in contrast, follows from the proved coefficient invariant.

---

## 9. The generic case, and what sets its constant

For the measured circuit instances with an odd factor — the generic and
cryptographically relevant cases in this construction — Walsh density rises
monotonically toward one half: 0.473 → 0.4994 over 15 to **30** qubits, with
growth **1.006 bits per qubit** (two series varying the modulus at fixed
exponent width). The measured growth is consistent with Θ(2ⁿ) — asymptotically
no better than state-vector simulation — with one caveat stated plainly: the
½ ceiling is proved (below), while "bounded away from zero" and the growth law
are measured at every size we can reach, not proved for every selected bit or
compilation. We state them as measurements, not extrapolation.

**The limiting constant of ½ is neither accidental nor algorithmic.**
The support carries a linear structure w = e_msb ⊕ e_anc, the 0/1 mask pairing
the accumulator's sign bit with the modular-reduction comparison ancilla, which
confines it to a hyperplane. That a linear structure confines the Walsh support
to a coset is classical — Carlet, *Boolean Functions for Cryptography and Coding
Theory*, Proposition 29 — and we **cite rather than claim** it; our contribution
is locating it here and reading it as a cost cap. It is observable-independent and
shared by both compilations, which use the same reduction discipline. So the
growth conclusion is robust while **the constant should be read as a property of
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

**N = 15 is the smallest relevant odd semiprime with this property, and all
seven of its usable bases are free**, against 3 of 11 for N = 21.

The canonical demonstration instance is therefore degenerate for classical
simulation not because of an unlucky base but because of the *modulus* — and the
property that makes 15 the natural smallest example is the same one that makes it
uninformative as a benchmark. Simulation results on N = 15 should not be
extrapolated.

### 10.2 Most of the support is inert for the real input state

Shor's actual initial state puts the exponent register in |+⟩. Since ⟨Z⟩ = 0
there, *every* Pauli term with support on the exponent register contributes
exactly zero to the expectation value. Combined with Corollary 3 this gives an
exact characterization of which terms contribute, but not an exact cardinality
fraction: the surviving sectors need not have equal support sizes. Measured
counts are **3883/15549 ≈ 0.24973** useful at N=7, a=6, α=1 and
**4014/32143 ≈ 0.12488** at N=5, a=2, α=2 (with analogous width-invariant
rows), close to but not exactly 1/4 and 1/8. The reported α=2 change occurs at
n_exp = α + 1 in these instances; this is a measured onset, not an independent
proof of a universal fraction or lower bound.

**This is not a speedup and we do not present it as one.** In the uncontracted
PPS implementation, the input state is applied only at the end of propagation,
and exponent support is not monotone under back-propagation (a term carrying it
can lose it and go on to contribute), so naive within-block term pruning is
unsafe. That does not exclude a sound contraction after a completed block. What
the observation does is separate *cost* from *useful work*: at α = 2, most of the
carried terms are irrelevant to the answer.

An exact state-aware extension now implements that completed-control
contraction. After the first forward gate touching a control has been passed
in reverse, all remaining gates avoid it, so contracting its independent |+⟩
input commutes with the remaining propagation. In the diagonal representation
this removes terms carrying Z on that control. For m work qubits and one fresh
control per contiguous block, at most m+1 qubits are active, giving a retained
support bound of 2^(m+1) in either order branch, including each block's temporary
pre-contraction peak. This is a different, reduced-operator task; it leaves the
full-operator spectra above unchanged. Repeated identity-tail involutions also
give an idempotent reduced map, allowing all but one such block to be omitted.
The implementation, coefficient-level checks and measurements are supplied in
`experiments/experiment_control_trace.py`.

The operation is ordinary last-use variable elimination, related to the
contraction-order viewpoint of Markov and Shi. Its significance here is the
explicit safe schedule and its effect on the interpretation of the arithmetic
split, not novelty of contraction itself. The support bound does not bound
gate count, coefficient precision, dictionary overhead or all allocated buffers.

### 10.3 Balanced controls and an exact tail residual

Let a unitary involution V act on the work register, controlled by fresh,
independent, identically distributed windows. If the Boolean activation
function is true with probability p, contraction gives
E=(1−p)Id+p Ad_V. Split the observable into O±=(O±V†OV)/2. These are the
two eigenspaces of conjugation by V, so direct iteration gives

> E^K(O) = O+ + (1−2p)^K O−.

The error on replacing the tail by its invariant part is exactly
|1−2p|^K ||O−|| in any homogeneous norm. At p=1/2 the map is an idempotent
projection after one window. Balanced nonlinear controls and balanced parity
controls therefore give the same reduced operator, despite potentially very
different full Walsh support. Constant controls need not converge.

The balanced-control experiment checks this formula on toy work maps and the
actual arithmetic identity-tail involution, with biased and constant controls.
It is a block/function-level comparison; it does not compile every Boolean
window into the arithmetic gate set. The result also does not cover correlated
or reused windows. This elementary two-eigenspace calculation complements the
full-operator criterion of §7.2: Boolean complexity can matter for the full
representation even when only its bias matters for the contracted task.

---

## 11. Cross-method corroboration

### 11.1 Shared arithmetic is not a shared lower bound

The same decomposition r = β·2^α appears in matrix-product-state simulation of
Shor's algorithm. Dang, Hill & Hollenberg (arXiv:1712.07311) give α as the
number of trailing zeros of r and discuss a β² memory factor; the decomposition
is explicit in their §4 and §5.2. Their tensor-partition resource statement is
not the same asymptotic assertion as PPS growth with exponent width, so this is
corroboration of a shared arithmetic parameter, not a proof of the PPS claim.

Two structurally unrelated classical methods — Heisenberg-picture Pauli
propagation and Schrödinger-picture tensor networks — keying on the same
arithmetic invariant is suggestive of a property of the algorithm rather than of
either simulator. This does not establish method-independent hardness: at fixed
N and r the pre-QFT state has exponent/work Schmidt rank at most r, and the
contraction in §10.2 also changes PPS width scaling. Further methods can test
the role of the arithmetic parameter, but agreement alone would not prove a
universal simulation lower bound (§13).

Paper A's tensor diagnostic makes the same distinction in another basis:
normalized Walsh transformation preserves cut singular values, while a scalar
periodic output can have much smaller rank than the full scratch-space
function. Scalar output rank, joint-state rank, and the memory needed to obtain
conditional probabilities must not be substituted for one another.

### 11.2 Conditional sampling through the inverse QFT

A measured |+⟩ control, with branch work permutations P0 and P1 and phase
feedback θ, induces the unnormalized work-state update

> M_b ψ = (P0ψ + (−1)^b exp(iθ)P1ψ)/2.

The norm squared is the outcome probability. Keeping the measured outcome
preserves the interference that the averaged channel in §10 loses. Processing
the largest exponent power first and applying inverse-QFT feedback gives a
one-history sampler without an exponent-register state vector. This is the
known semiclassical/iterative route of Griffiths and Niu, not a new Fourier
transform algorithm.

Reordering the compiled multiplications is valid on the reachable clean
arithmetic subspace because modular multiplications commute there and preserve
that subspace. The arbitrary full-scratch maps need not commute. An explicit
noncommuting-gate control detects the error made by ignoring this qualification.

The conditional-order-finding experiment validates complete small output
distributions against an independent FFT and Fourier-compiled state-vector
circuits, and individual wider sampled probabilities against a geometric-series
reference. The sampler is not supplied the order; only the independent
reference is. Its first implementation extracts full-work transition tables
and stores a dense work vector, so both setup and work memory are exponential
in work width. Output-history enumeration is used only for small validation
cases. This crosses the inverse QFT as a correctness baseline, not as a
scalability result. The open research task is reducing total setup and
conditional-state cost while comparing against existing orbit and tensor
methods on the same sampling task.

### 11.3 On-demand compiled arithmetic and its orbit-sized limit

The reachable-state extension replays the existing classical gate kernel only
on queried basis labels. It constructs neither a full work-space transition
table nor a dense scratch-state vector, and receives neither the order nor an
enumerated orbit. Both control branches are replayed, preserving the distinction
between clean arithmetic and arbitrary scratch behavior.

After s descending-power controls of a t-bit experiment, all reached work
labels have the form a^(2^(t−s)j) mod N for 0≤j<2^s. The union before outcome
selection therefore has size at most

> min(2^s, r/gcd(r,2^(t−s))).

This follows directly from the processed binary subset sums and the order of
the generated subgroup. Intermediate classical gates preserve each branch's
number of labels, even while scratch is dirty; completed blocks return to the
clean subspace. The bound is orbit-sized, not polynomial in modulus bit length.
Memoizing d distinct block maps can additionally store O(dr) transitions.

The reachable-order-finding experiment checks complete small distributions,
wide sampled paths, forced rare outputs, block correctness on all clean inputs,
and the per-step bound. It separately records compilation, cold replay, warm
sampling, and transition-cache costs. Ordinary high-level modular arithmetic
and orbit-enumeration baselines are faster in the reported pilot. Removing a
full scratch allocation is a useful implementation improvement, but not an
advantage over those classical methods. Numerical amplitude counts are also
not measurements of total process memory.

### 11.4 A latent eigenphase is enough once the order is known

A stronger classical baseline follows from the familiar spectral account of
order finding (Cleve, Ekert, Macchiavello and Mosca, §6). On the length-r orbit,
modular multiplication has eigenvalues exp(2πik/r), and |1⟩ has equal squared
overlap with their eigenvectors. Every measured work instrument is a polynomial
in this same unitary, so all output-history probabilities equal a uniform
mixture over k. Sample k once; conditional on it, each bit is just a scalar
Bernoulli trial with probability

> |1 + (−1)^b exp(i(θ + 2πk2^i/r))|² / 4.

No work-state vector is needed. The benchmark explicitly pays for discovering
r by repeated modular multiplication until returning to 1, using O(r)
arithmetic steps but no orbit table, then performs O(t) scalar updates per
sample. This is not efficient order discovery in general, nor a constant
bit-complexity or precision assertion. The conditional probability for a fixed
k must not be confused with the marginal output probability; the experiment
checks the entire mixture on small cases against independent references.

The coherent pre-QFT state's exponent/work Schmidt rank is nevertheless
min(r,2^t): each distinct work label is paired with a disjoint, nonempty
exponent class. Large state rank and cheap order-informed sampling coexist.
This identifies the distinction more sharply than a sparsity comparison:
compact prediction can depend on spectral information that is costly to obtain.
The argument is standard phase estimation, not a new classical factoring
algorithm. It also relies on the commuting ideal arithmetic; it is not a
replacement work state for arbitrary later operations.

---

## 12. Limitations and the honesty record

### 12.1 Scope

- **Constructions.** The proof of Theorem 2 uses the multiply–swap–unmultiply
  form, so it covers both compilations here and Vedral/Beauregard-style
  constructions, plus anything satisfying the corrected criterion of §7.2. It
  does not transfer automatically to a modexp built otherwise. The dead-tail
  statement of §7.3, and its 2-periodicity under the same involution hypothesis,
  hold with either a unitary or a measurement-based unlookup of the lookup register;
  the sizes in §7.3 are for the unitary unlookup.
- **Observables.** The full-support theorems concern computational-basis work
  observables before the inverse QFT. The conditional baseline of §11.2 is a
  different, state-aware sampling calculation, not a non-diagonal PPS theorem.
- **Sizes.** We distinguish modulus bit length, exponent-register width, and
  total physical qubits: the Walsh series reaches 30 total qubits, while peak
  memory is measured at n_exp = 8. The α = 3 onset is confirmed at two moduli
  and the α = 4 onset at one (§6); α ≥ 5 needs r = 32, whose smallest instance
  is out of reach.
- **No implication for factoring.** Efficient classical simulation of Shor's
  algorithm on general inputs would be a classical factoring algorithm. Nothing
  here bears on that; these are diagnostic results about where PPS breaks, and
  the free branch is exactly the branch where the order is already easy to find.

### 12.2 Claims withdrawn, and corrections

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
uv run python -m experiments.experiment_control_trace
uv run python -m experiments.experiment_balanced_controls
OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_conditional_order_finding
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_reachable_order_finding
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
   suggestive to settled.
2. **Where the cheap spectral description stops working.** The on-demand
   baseline in §11.3 removes full-scratch tables but does not outperform the
   charged classical baselines, including the latent-eigenphase route of §11.4.
   A sharper next experiment inserts one precisely specified noncommuting
   operation between arithmetic blocks, keeps its location fixed, and varies
   only its strength. Determine which spectral coherences the chosen output
   task now needs before designing a larger tensor framework. This is a
   perturbed-circuit question, not ideal Shor order finding. Larger ideal-orbit
   runs alone do not resolve it; the broader compressed-prediction problem
   remains open.
3. **Measurement-based uncomputation beyond the lookup register.** The
   measurement-based unlookup of the lookup register is covered (§7.3). The same
   Heisenberg rule applies to a measurement-based uncomputation of the activation
   ancillae, which would change the full-space function again and is not modelled or
   measured here, and the support sizes under the measurement-based unlookup are
   not reported.
4. **Whether the MPS correspondence extends** to other contraction orders and
   tensor-network methods, or is specific to these two.

---

## 14. Conclusion

For the full, uncontracted diagonal Pauli representation of the studied modular
exponentiation, arithmetic governs width scaling. Writing r = β·2^α, the exponent register is free when
β = 1 — provably once n_exp ≥ α+1, with first locking at that width measured in
the reported rows — and shows measured Θ(2ⁿ)-type growth otherwise. The free branch is reachable for *every* base of certain
moduli, which makes the field's canonical benchmark instance degenerate for
reasons that have nothing to do with the base chosen.

The criterion we extract is checkable rather than a family restriction, and
applying it to the arithmetic people actually propose to run reveals both that
the invariance survives — in a sharper form than we proved — and that a single
obvious optimisation would destroy it.

State-aware contraction changes the conclusion for the reduced pre-QFT task:
retained support is bounded independently of exponent width in both branches,
and balanced independent controls project repeated involutions immediately.
The conditional sampler then provides a separate, validated post-QFT task.
Together these findings distinguish a representation-cost mechanism from
unavoidable classical prediction difficulty. They strengthen the diagnostic
claim; they do not establish a simulation breakthrough.

---

## References

- Beauregard, *Circuit for Shor's algorithm using 2n+3 qubits*,
  quant-ph/0205095.
- Carlet, *Boolean Functions for Cryptography and Coding Theory*, Cambridge
  University Press, 2021.
- Cuccaro, Draper, Kutin & Moulton, *A new quantum ripple-carry addition
  circuit*, quant-ph/0410184.
- Cleve, Ekert, Macchiavello & Mosca, *Quantum Algorithms Revisited*,
  [quant-ph/9708016](https://arxiv.org/abs/quant-ph/9708016), §6.
- Dang, Hill & Hollenberg, *Optimising Matrix Product State Simulations of
  Shor's Algorithm*, arXiv:1712.07311.
- Draper, *Addition on a Quantum Computer*, quant-ph/0008033.
- Gidney, *Windowed quantum arithmetic*, arXiv:1905.07682.
- Griffiths & Niu, *Semiclassical Fourier Transform for Quantum Computation*,
  [quant-ph/9511007](https://arxiv.org/abs/quant-ph/9511007).
- Markov & Shi, *Simulating quantum computation by contracting tensor networks*,
  [quant-ph/0511069](https://arxiv.org/abs/quant-ph/0511069).
- Vedral, Barenco & Ekert, *Quantum Networks for Elementary Arithmetic
  Operations*, quant-ph/9511018.

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
| 6. The proved threshold and measured first-lock rows | C21 |
| 7.1 The parity reduction | C23, C24, C43 |
| 7.2 The criterion, and a correction to it | C29 (narrowed) → C36 |
| 7.3 Windowed arithmetic | C36, C37, C38, C39, C105 |
| 8. Robustness to truncation | C27, C28 |
| 9. The generic case, and what sets its constant | C7, C30, C31, C32 |
| 10.1 N = 15 is degenerate by construction | C20, C22 |
| 10.2 Most of the support is inert for the real input state | C42, C45, C46 |
| 10.3 Balanced controls and an exact tail residual | C47 |
| 11.1 Shared arithmetic is not a shared lower bound | C19; Paper A: C48 |
| 11.2 Conditional sampling through the inverse QFT | C49 |
| 11.3 On-demand compiled arithmetic and its orbit-sized limit | C51 |
| 11.4 A latent eigenphase is enough once the order is known | C52 |
| 12.2 Claims withdrawn, and corrections | C29, and the retracted rows of the ledger |

---

## Appendix B — proofs of Theorem 2 and Corollary 3

**Setting.** The controlled-multiplier block attached to exponent bit i is
built in the multiply–swap–unmultiply form

> u(ctrl, c) = cmult(ctrl, c) ; cswap layer ; cmult(ctrl, c⁻¹)⁻¹,

with c = a^{2^i} mod N. For i ≥ α and β = 1 we have c = 1, so the two multiply
stages are mutual inverses and every such block applies the same operator

> V := u(·, 1) = A⁻¹ S A,

where A is the controlled multiply-accumulate and S the controlled-swap layer,
each controlled on the block's own exponent qubit. Write I = {i : i ≥ α} for
the identity tail; the tail blocks are the final |I| blocks of the circuit, so
their product is a suffix of it.

*Proof of Theorem 2.* **(i) V is an involution.** S is the product over
register position k of cswap(ctrl, x_k, b_k). The target pairs (x_k, b_k) are
pairwise disjoint, so the factors commute and each squares to the identity;
hence S² = id. Conjugation preserves order:
V² = (A⁻¹SA)(A⁻¹SA) = A⁻¹S²A = id.

**(ii) The tail composes to a parity.** V acts as the identity on every
exponent qubit — exponent qubits enter the block only as controls — so a
controlled application C_{e_i}(V) neither reads nor writes any e_j with j ≠ i.
On a computational-basis state with tail bits (e_i)_{i∈I}, the product
Π_{i∈I} C_{e_i}(V), in any order, therefore acts as V^{Σ_{i∈I} e_i}, and since
V² = id this is V^p with p = ⊕_{i∈I} e_i. The circuit depends on the entire
tail through that single bit. ∎

*Proof of Corollary 3.* By Theorem 2 the pulled-back bit function factors
through the parity: writing the full input as (y, e) with e = (e_i)_{i∈I} the
tail bits and y all remaining qubits (n_y of them, exponent bits below α
included), g(y, e) = F(y, p(e)) for some Boolean function F. The Walsh
coefficient at (z_y, z_I) is

> ĉ = 2^{−n} Σ_y (−1)^{y·z_y} Σ_{e ∈ F₂^{|I|}} (−1)^{F(y, p(e))} (−1)^{e·z_I}.

Group the inner sum by parity. Let K = {e : p(e) = 0}, the kernel of the parity
form — a subgroup of index 2 — and K¹ its coset. The character
χ(e) = (−1)^{e·z_I} sums to zero over K, and over K¹, unless χ is trivial on K,
i.e. unless z_I lies in the annihilator K^⊥; since K is the kernel of the
all-ones form, K^⊥ = {0, 1_I}. Hence ĉ = 0 whenever z_I ∉ {0, 1_I}, and for
the two surviving values the inner sum is

> z_I = 0:  2^{|I|−1} [ (−1)^{F(y,0)} + (−1)^{F(y,1)} ]
> z_I = 1_I: 2^{|I|−1} [ (−1)^{F(y,0)} − (−1)^{F(y,1)} ]

(for z_I = 1_I, use (−1)^{e·1_I} = (−1)^{p(e)}). With n = n_y + |I|,

> ĉ_{(z_y, 0)}  = 2^{−n_y−1} Σ_y (−1)^{y·z_y} [ (−1)^{F(y,0)} + (−1)^{F(y,1)} ]
> ĉ_{(z_y, 1_I)} = 2^{−n_y−1} Σ_y (−1)^{y·z_y} [ (−1)^{F(y,0)} − (−1)^{F(y,1)} ]

— expressions in which |I| does not appear. Writing each support element as
(z_y, flag), with the flag recording whether z_I = 0 or 1_I, the support is
therefore the same *set* at every tail length, each of its two halves is
individually constant, and its size carries no dependence on the number of
tail qubits.

Where the hypotheses fail, so does the confinement: for i < α the block
multiplies by a^{2^i} ≠ 1 and is not V, and for β > 1 no block is the identity
on the valid subspace at all, so distinct blocks apply distinct permutations
and no parity reduction is available. This is the growth below the onset and
in the generic branch. ∎
