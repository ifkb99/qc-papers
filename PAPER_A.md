# Walsh–Hadamard Sparsity Exactly Determines Pauli-Path Simulation Cost for Reversible Quantum Arithmetic

**Ian Baker**

*Draft v1, 2026-08-08. Status markers are load-bearing and must survive to
submission; see §11. Every quantitative claim carries its ledger ID from
`CLAIMS.md`, which is the single source of truth for claim status. Numbers
here are post-bugfix (see §11.2).*

---

## Abstract

Pauli Path Simulation (PPS) — also called sparse Pauli dynamics or Pauli
propagation — has become a leading method for classically simulating
utility-scale quantum circuits. Its practical deployment depends on predicting,
before committing to an expensive run, how many Pauli terms a circuit will
require: a question currently answered by empirical power-law extrapolation
calibrated on brickwork and Trotterised circuits with generic rotation angles.

We show that for circuits implementing a permutation of the computational basis
— all reversible arithmetic, and hence the bulk of Shor's algorithm and of
quantum algorithms generally — this quantity is not merely predictable but
*exactly computable in closed form*, with no extrapolation and no fitting. For a
circuit implementing basis permutation π, the Heisenberg pullback π†Z_jπ is the
diagonal operator (−1)^{g(y)} with g(y) = bit j of π(y); expanding a diagonal
operator in the Pauli basis is precisely the Walsh–Hadamard transform of
(−1)^g. The Pauli support PPS must carry is therefore *exactly* the Walsh
spectrum of g, and its size is the Walsh sparsity of that Boolean function.

Four consequences follow. **(i)** PPS cost for reversible arithmetic is a
property of the **full-space basis permutation implemented** — including its
action on ancillas — and not of the gate set implementing it: we verify exact
Z-closure for both Toffoli-compiled and Fourier-compiled (Beauregard) modular
exponentiation, which implement the same permutation by very different means.
**(ii)** The tractability of linear arithmetic is explained rather than observed
— a ripple-carry adder's low output bit is XOR-affine, Walsh sparsity 1
characterises affineness, and PPS collapses to a single term. **(iii)** Peak
memory is a distinct and larger quantity than final support, and its excess is an
artifact of Clifford+T decomposition; propagating X, CNOT and Toffoli as atomic
permutations keeps the expansion Z-type at every step, reduces peak memory by a
measured factor of 2.000 (adders) to 1.9997 (modular exponentiation), and makes
the peak itself a Walsh quantity. **(iv)** Walsh sparsity and nonlinearity are
the same object linear cryptanalysis studies, giving a transfer: any published
nonlinearity lower-bounds PPS cost for *every* circuit computing that function,
with no simulation.

The quantity the model computes is the size of the **Heisenberg representation
PPS maintains** — its memory footprint — not the difficulty of the expectation
value it is used to estimate; §1.1 states the target task and input-state regime
explicitly, since for some inputs the expectation is obtainable by other means
entirely.

These results are diagnostic, not a simulation speedup. The Walsh transform is
itself exponential, and nothing here bears on the classical hardness of
factoring.

---

## 1. Introduction

Pauli propagation simulates a quantum circuit in the Heisenberg picture by
expanding an observable in the Pauli basis and pushing it backwards through the
circuit gate by gate. Each non-Clifford gate can split a Pauli string into
several, so the represented operator grows, and practical runs control this with
a truncation threshold δ that discards small coefficients. The method's success
on utility-scale experiments has made a single practical question urgent: *given
a circuit, how many Pauli terms will it need?*

The existing answer is empirical. Resource models are calibrated on brickwork
and Trotterised Hamiltonian dynamics with generic rotation angles, where
coefficient magnitudes follow a power law and the cost/accuracy trade-off is
predictable by extrapolation. Reversible arithmetic circuits differ from that
family on every axis that the calibration assumes: their rotation angles are
highly correlated (exactly ±π/4 and ±π/2 in a Clifford+T compilation), their
connectivity is long-range, and their structure is deterministic rather than
random.

We show that for this family the question has an exact answer rather than a
better extrapolation.

### 1.1 What is being predicted, and for which task

A cost model must say what it costs. Ours predicts the **size of the Heisenberg
representation PPS maintains** — the number of Pauli terms held in memory — at
δ = 0. That is a statement about a data structure, and it is deliberately
separate from two things it could be confused with.

**It is not a claim that the expectation value is hard to obtain.** PPS estimates
⟨ψ|π†Z_jπ|ψ⟩ for a specific input state, and for some input states that number is
available by other means. For a computational-basis input |y⟩ the operator is
diagonal with ±1 entries, so ⟨O⟩ = (−1)^{g(y)} *exactly and trivially* — one
evaluation of π settles it. At the opposite extreme, for |+⟩^n every Z^z with
z ≠ 0 has zero expectation, so only ĉ_0 survives and the sum collapses to a
single coefficient. In both regimes the *answer* is cheap while the
*representation* is not, and it is the representation that determines whether a
run fits in memory.

**Why that is the useful quantity.** Permutation circuits occur as subcircuits of
larger circuits that are not permutations — Shor's algorithm surrounds its
arithmetic with Hadamards and an inverse QFT — and a PPS run through the whole
thing must carry the arithmetic block's representation whether or not the
surrounding structure eventually collapses it. Memory is the binding constraint
in practice, and it is set by the peak representation size, not by the difficulty
of the final number. A companion result quantifies exactly this gap: for the
genuine Shor initial state, a determinate fraction of the carried terms
contribute nothing to the expectation value, yet cannot be discarded early
because the property is not monotone under back-propagation.

Readers who want a sharper "so what" should read §5 (peak memory, where the model
changes what one should actually do), §7 (truncation, where it exposes a
pathology), and §8 (a bound requiring no simulation at all).

**Contribution.** The technical core is a single identity (§3) whose ingredients
are individually standard and whose composition, as far as we can find, is not
stated: for permutation circuits with computational-basis observables, the
Pauli-spectrum/Boolean-Fourier analogy becomes an *exact identity*, and that
identity is a cost model. What we claim is that composition and what follows
from it — compilation invariance, the affine explanation of adder collapse,
permutation-native propagation and its near-exact factor of two, the cryptanalytic
transfer, and the structural caps of §6. We are explicit in §10 about which
ingredients are prior art.

---

## 2. Preliminaries

### 2.1 Pauli propagation

We write a Pauli string as a pair of bitmasks (x, z) with

> P(x,z) = i^{|x ∧ z|} X^x Z^z,

the phase chosen so that P is Hermitian; all coefficients in a Heisenberg
expansion are then real. Every gate is a Pauli rotation U = exp(−iθσ/2), with a
single conjugation rule:

> U† P U = P                              if [P, σ] = 0
> U† P U = cos θ · P + sin θ · (iσP)      if {P, σ} = 0.

Clifford gates have θ = ±π/2, so cos θ = 0 and a Pauli maps to a *single* Pauli
— no branching. T gates have θ = ±π/4, giving cos = sin = 1/√2, the maximum
branching weight. A circuit's PPS cost is driven by how much of that branching
survives cancellation.

### 2.2 Walsh–Hadamard transform

For a Boolean function g on n bits, the Walsh–Hadamard transform of its
character form (−1)^g is

> ĉ_z = 2^{−n} Σ_y (−1)^{g(y)} (−1)^{y·z}.

The **Walsh support** is {z : ĉ_z ≠ 0} and the **Walsh sparsity** S is its
cardinality. Parseval gives Σ_z ĉ_z² = 1. The **nonlinearity** of g is
NL = 2^{n−1}(1 − max_z |ĉ_z|). Sparsity 1 characterises affine functions; a flat
spectrum (all |ĉ_z| = 2^{−n/2}) characterises bent functions, which have maximal
nonlinearity and full support.

---

## 3. The identity

### 3.1 Statement and proof

**Theorem 1 (C8).** *Let a circuit implement a permutation π of the
computational basis, so that the unitary acts as |y⟩ ↦ |π(y)⟩. Let Z_j be a
computational-basis observable. Then*

> π† Z_j π = Σ_z ĉ_z Z^z,   with   ĉ_z = 2^{−n} Σ_y (−1)^{g(y)} (−1)^{y·z},
> *where* g(y) = bit j of π(y).

*The pullback is diagonal, expands purely in Z-type Pauli strings, and its
coefficient vector is exactly the Walsh–Hadamard transform of (−1)^g. In
particular the number of Pauli terms PPS must carry is the Walsh sparsity of g.*

*Proof.* π†Z_jπ|y⟩ = π†Z_j|π(y)⟩ = (−1)^{bit_j(π(y))} π†|π(y)⟩ = (−1)^{g(y)}|y⟩.
The operator is therefore diagonal with entries (−1)^{g(y)}. Any diagonal
operator D expands as D = Σ_z d̂_z Z^z with d̂_z = 2^{−n} Σ_y D_{yy} (−1)^{y·z},
since ⟨y|Z^z|y⟩ = (−1)^{y·z} and the Z^z are orthogonal under the normalised
trace inner product. Substituting D_{yy} = (−1)^{g(y)} gives the claim. ∎

The proof is elementary. Its content is not the derivation but the
identification: the quantity PPS practitioners estimate by extrapolation is a
named, exactly computable invariant of a Boolean function.

### 3.2 Which Boolean function — a definition that must be stated precisely

Throughout, g is bit j of π(y) where **y ranges over the entire register,
ancillas and scratch included**, and π is the **full-space** permutation. This is
not a technicality; getting it wrong makes two of our results look contradictory.

Two circuits can compute the same arithmetic result on the subspace one cares
about — the "valid subspace" where scratch registers start and end at |0⟩ — while
implementing *different* full-space permutations, because they leave different
intermediate garbage off that subspace. Their g's differ, so their Walsh spectra
differ, and Theorem 1 assigns them different costs. That is a consequence of the
theorem, not an exception to it.

We therefore state compilation invariance in the only form that is true:

> **PPS cost is invariant among circuits implementing the same full-space
> permutation π, and is not invariant among circuits that merely agree on the
> valid subspace.**

§4.1 exhibits the first case: two compilations of modular exponentiation, built
from entirely different gate sets, implementing the same π and having the same
cost. §6.3 exhibits the second: two circuits computing the same a^e mod N whose
π's diverge off the valid subspace, with a ~51% cost difference. Both are
predictions of Theorem 1.

**Verification (C8).** Confirmed to machine precision on 6/6 instances spanning
modular exponentiation and ripple-carry addition, across both compilations, with
**maximum error ≤ 6.7 × 10⁻¹⁶**. We check not merely that the counts agree but
that the *support sets are identical* — the stronger statement, and the one that
would fail first under a coincidence.

### 3.3 Scope: exactly where this holds and where it stops

The theorem needs two things: the unitary is a basis permutation, and the
observable is diagonal.

- **Diagonal observables of any weight are covered (C13).** Multi-qubit Z-type
  observables work identically; verified exact on 4/4 Z-type observables tested.
- **X and Y observables are not covered (C13, F10).** Their pullbacks leave the
  diagonal and the expansion is fully non-diagonal. This is a hard boundary, not
  a gap in the analysis, and it is where the present paper stops. Real Shor
  measures the exponent register after an inverse QFT, which is outside this
  scope; we do not address it.
- **Permutation-ness is required of the *unitary*, not of the gate set (C6).**
  A circuit whose individual gates are not permutations may still implement one,
  and the identity applies. This is the point of §4.1.

---

## 4. Consequences

### 4.1 Compilation invariance (C1, C2, C6)

Because the identity depends only on π and j, PPS cost is a property of the
full-space permutation implemented, not of the gates used to implement it (§3.2).
This section is the invariance half; §6.3 is the non-invariance half, and they are
the same statement applied to circuits that do and do not share a π.

This contradicts a natural intuition — that Clifford+T compilation of reversible
arithmetic is distinguished, because such a circuit is a permutation matrix
*gate by gate* and so trivially preserves Z-type strings, whereas Fourier-space
arithmetic (Draper addition, as in Beauregard's 2n+3 qubit construction) is
built from controlled-phase and Z-rotation gates that are individually not
permutations. We verify that **both** compilations exhibit exact closure of the
Z-type Pauli subalgebra (0 non-Z terms in each), because both implement the same
basis permutation.

> **Retraction note.** An earlier version of this work claimed the opposite —
> that compilation determines simulability, with a four-order-of-magnitude gap
> between the two compilations. That claim was false, killed by a propagator bug
> (θ = π gates treated as no-ops) compounded by an ancilla-count confound in the
> matched comparison. It is retracted in full. See §11.2.

### 4.2 Affine collapse is exactly sparsity one (C10)

It is folklore that adders are easy for Heisenberg-picture methods. The identity
makes this exact: sparsity 1 characterises affine functions, so a PPS collapse
to a single term is *equivalent* to the output bit being XOR-affine.

For a 4-bit Cuccaro ripple-carry adder, the low output bit b₀ = a₀ ⊕ b₀ is
affine and its Walsh sparsity is **exactly 1**. The control is the discriminating
half: bit b₂ carries carries, is not affine, and has sparsity **10**. The
collapse is therefore not about permutation-ness — an earlier explanation we
retract — but about algebraic degree.

### 4.3 A cost model at O(2ⁿ n)

Computing the Walsh spectrum requires one pass to extract π and one fast
Walsh–Hadamard transform, i.e. **O(2ⁿ n)** time. This is exponential and no
speedup; the point is that it is *cheaper than the run it predicts*, by roughly
three orders of magnitude on our instances (C11), because it does no branching
bookkeeping and no coefficient arithmetic. It answers "will this run fit in
memory" before the run.

Two things must be stated with it. First, the model gives the **final** support
exactly; peak memory is a different and larger quantity (§5). Second, it is
exact only at δ = 0; under truncation the relationship is more subtle (§7).

---

## 5. Peak versus final cost, and permutation-native propagation (C17, C18)

The final Pauli support is what Theorem 1 computes. What bounds memory in
practice is the *peak* over the propagation.

The two differ for a reason that is entirely an artifact of compilation. The
standard Clifford+T decomposition of a Toffoli passes through Hadamard gates, so
a Z-type string leaves the diagonal mid-circuit and returns to it only at the
end of the gadget. The excursion is real, costs memory, and cancels exactly.

Propagating X, CNOT and Toffoli as **atomic permutation primitives** removes the
excursion. The expansion is then Z-type at *every* intermediate step, and:

- peak memory drops by a measured factor of **2.000000** for ripple-carry adders
  (128→64, 512→256, 2048→1024) and **1.9997** for modular exponentiation
  (13666→6834, 16386→8194, 128138→64070). The modexp instances satisfy
  `rot = 2·perm − 2` exactly, in 3 of 3 cases;
- the peak becomes a Walsh quantity in its own right — the maximum over circuit
  suffixes of the corresponding sparsity — so it is predictable by the same
  model;
- propagation runs about **10× faster** than rotation-level propagation on the
  same circuits, since no branch is created only to be cancelled.

**Is the factor of two a theorem or an observation? Currently an observation,
and we state it as one.** The mechanism is clear enough to suggest a proof: the
gadget's Hadamard on the target qubit c exchanges the Z-sector with an X-sector,
and a Z-type string containing Z_c is carried into a mirrored partner while a
string without Z_c is untouched — so the doubling applies to the Z_c-containing
subset and bounds the ratio above by 2. That accounts for the adders hitting
2.000000 exactly and for modexp falling just short. It does **not** yet account
for the deficit being exactly 2 in all three modexp instances, which is a clean
regularity we have not explained. Until it is proved we claim only:
*the ratio is at most 2, is attained for the adders tested, and is
2 − O(1/peak) for the modular exponentiation instances tested.*

This is a genuine practical recommendation regardless: for permutation circuits,
do not decompose to Clifford+T before propagating.

---

## 6. Structure in the spectrum: why arithmetic pullbacks are not generic

The identity converts PPS cost questions into Boolean-function questions, which
lets known structure do work. Reversible-arithmetic pullbacks turn out to be
markedly non-generic, and the deviations are exactly characterisable.

### 6.1 Linear structures cap the density (C30, C31)

A **linear structure** of g is a w with g(y ⊕ w) = g(y) for all y. The relevant
classical fact is due to Carlet (*Boolean Functions for Cryptography and Coding
Theory*, **Proposition 29**): D_e f is null (resp. constant 1) iff supp(W_f) is
contained in {0,e}^⊥ (resp. in its complement). A linear-structure space of
dimension k therefore confines the Walsh support to a coset of codimension k,
capping density at 2^{−k}. Functions whose Walsh support is an affine subspace
are Carlet's **partially bent** functions.

**We cite this rather than claim it.** Our contribution is locating it in
modexp pullbacks and reading it as a cost cap. Modular exponentiation carries
the linear structure w = b_msb ⊕ anc: GF(2) rank n−1 in every instance,
g(y ⊕ w) = g(y) verified pointwise, present in *both* compilations, and absent
from random f. Consequently density is capped at exactly ½ — which is why
measured density converges to 0.498 *from below* and never crosses (C7). A
ripple-carry adder has a kernel of dimension 5, capping its density at 2⁻⁵.

**Mechanism (C31).** The structure is forced by three ingredients acting
together: flipping the accumulator's msb *is* adding 2^{m−1}, which commutes
with mod-2^m addition (verified exhaustively, 0/256 violations); the ancilla is
coupled to the msb only by XOR, so flipping both restores it; and the msb is
excluded from the swap network, so it never reaches the observed register.

### 6.2 A conditional generalisation (C40, C41)

The above has a natural extension we did not find stated, and which we present
at corollary altitude: it follows from Proposition 29 plus the standard
decomposition of the Walsh transform over a coset partition.

**Proposition 2 (C40).** *Let a coset partition split F₂ⁿ into cells H_u, and
suppose g restricted to H_u has linear structure w_u. Then the Walsh support
avoids E = {z : w_u·z = 1 for every u}, and hence*

> density ≤ 1 − 2^{−d},   d = dim span{w_u},

*provided the system {w_u·z = 1} is consistent.*

*Proof sketch.* Decomposing over the partition gives
c_z = 2^{−n} Σ_u (−1)^{u·z_C} A_u(z′) with A_u the transform of the restriction;
Proposition 29 kills A_u(z′) whenever w_u·z′ = 1. If every A_u vanishes, so does
c_z. ∎

**Proposition 3 (C41).** *That system is consistent iff every linear dependency
among the w_u has even support. An odd dependency makes E empty and removes the
cap entirely.*

Verified on planted instances at d = 1, 1, 2, 3 with caps 0.5, 0.5, 0.75, 0.875
and **zero violations**; shown invariant under a random GF(2) change of basis
(tested, not asserted), so it covers arbitrary subspace partitions. The parity
condition is sharp: at the same cell count, an odd dependency gives density
0.9845 against 0.8671 for an even one.

The cap is governed by the **span dimension** of the per-cell structures, not by
the depth of conditioning. This matters practically: it says which modifications
to a construction can possibly change its cost — only those adding an
independent structure vector, or introducing an odd dependency.

### 6.3 A compilation choice with a real cost (C32)

§4.1 showed cost is invariant among circuits sharing a full-space permutation.
This section is the other half of §3.2: two circuits that agree on the *valid
subspace* but not off it, and therefore do not share a π.

The standard Beauregard/Vedral modular reduction is *accidentally* friendly to
PPS: its ½ cap is a byproduct of the msb being excluded from the swap network
while commuting with the adder. Replace the reduction with one that touches the
msb *nonlinearly* — a single Toffoli conjugating the reduction, controlled on
scratch qubits that are |0⟩ on the valid subspace — and the linear structure is
destroyed, raising density from 0.473 to 0.716: a **~51% cost increase**. Linear
couplings (CNOT, cswap) leave it intact; only nonlinearity breaks it.

**Where the two permutations diverge, precisely.** The modified circuit computes
the same a^e mod N — verified for every exponent, with all scratch returned to
|0⟩ — and the two circuits agree on every valid input. But they are **different
full-space permutations**: at N = 5, a = 2, n_exp = 2 they differ on **8554 of
32768 basis states (26.1%)**, all of them off the valid subspace, where the
modified reduction leaves different garbage in scratch. Their g's therefore
differ, and Theorem 1 correctly assigns them different Walsh spectra — 15493
against 23464 terms.

This is worth stating explicitly because the loose version of the sentence —
"a modification that changes nothing about the computed function" — appears to
contradict §4.1 and does not survive contact with the definition in §3.2. It
changes nothing about the *arithmetic result on the valid subspace*; it changes
π, and hence g, and hence the cost. Note also that the quantity differing here is
**final support**, exactly what Theorem 1 computes, not a peak or intermediate
quantity — so this is a prediction of the theorem rather than a measurement of
something else.

The framing worth keeping is the reverse of the usual one: not "compile to X to
make PPS cheaper", but "the standard construction is already cheaper than it
needs to be, for a reason nobody designed". This is a constant-factor effect;
the asymptotics are unaffected.

The practical reading is sharper than it first appears: since scratch behaviour
off the valid subspace is invisible to correctness testing, two implementations
that pass identical verification suites can differ by ~51% in PPS cost. Ancilla
discipline is a cost parameter, not merely a hygiene concern.

---

## 7. Truncation: what the model does and does not say

The identity is exact at δ = 0. Practical PPS truncates, and the model's
relationship to truncated runs is not the naive one.

- **Coefficient (δ) truncation is the right knob (C16).** The spectrum is
  heavy-tailed in a specific way: on a representative instance the 4 terms with
  |c| > 0.1 (of 3086) already reproduce ⟨O⟩ exactly, and the remaining 3082 sum
  to **exactly** 0.
  *Read this carefully — the total is not the surprising part.* For a
  computational-basis input the operator is diagonal with ±1 entries, so
  ⟨O⟩ = ±1 is guaranteed a priori (§1.1) and "−1.00000000" is not evidence of
  anything. The content is the **concentration**: that 4 coefficients of 3086
  suffice, and that the other 3082 cancel to zero rather than merely to
  something small. That is a structural statement about the spectrum, and it is
  what makes δ-truncation viable here.
- **Weight truncation is the wrong knob**, and structurally so. Signed sums per
  weight level are large and alternating, cancelling only over all levels, so any
  cutoff slices the cancellation; the dominant coefficients sit at weights 1, 2,
  8 and 9, not all low-degree. The deeper reason is that Hamming weight is *not
  an affine invariant* while the Walsh support and |c_z| are — so a weight-graded
  truncation is basis-dependent by construction. (Independently, Gangopadhyay et
  al., *J. Appl. Math. Comput.* 69:3337–3357, 2023, show the analogous
  weight-graded spectrum is not invariant under extended affine equivalence.)
- **Truncation error is non-monotonic in δ (C14).** On one instance δ = 10⁻¹ is
  *exact* with 34 terms while δ = 10⁻³ is off by 0.285 with 23482 terms. Small
  coefficients cancel as a set; removing some of them is worse than removing all.
- **Truncated estimates can be inadmissible (C5).** We observed |⟨O⟩| > 1 in
  4 of 18 runs, all at *mild* δ. We propose the operator-norm bound as a
  zero-cost admissibility check on PPS output.

The honest summary: the exact model describes the untruncated support, and
truncation on this family behaves badly enough that it should be monitored
rather than trusted.

---

## 8. The cryptanalytic bridge (C25, C26, C12)

Walsh sparsity and nonlinearity are the central quantities of linear
cryptanalysis. The identity turns that coincidence into a transfer.

**Proposition 4 (C25).** *For any Boolean function g with nonlinearity NL,*

> S ≥ (1 − NL/2^{n−1})^{−2}.

*Proof.* NL = 2^{n−1}(1 − max_z|ĉ_z|) gives max_z|ĉ_z| = 1 − NL/2^{n−1}.
Parseval gives 1 = Σ_z ĉ_z² ≤ S · max_z ĉ_z². Rearranging gives the bound. ∎

The consequence is that **any published nonlinearity is a lower bound on PPS
cost for every circuit computing that function**, obtained with no simulation
and independent of compilation.

- **Tight at both extremes.** Affine functions give S ≥ 1 (attained); bent
  functions give S ≥ 2ⁿ (attained). We verify the bent end **end to end**,
  building reversible circuits computing an inner product and propagating them:
  support is exactly 2^{2m} — 64, 256, 1024 — matching the prediction. When every
  coefficient has the same magnitude, no truncation is available at all; this is
  the clean worst-case statement.
- **External validation.** The AES S-box nonlinearity comes out at exactly
  **112** over all 255 nonzero linear combinations, reproducing the published
  constant. This is an independent check on the whole Walsh pipeline, not just
  on the bound.
- **Honest limitation (C26).** The bound is loose away from the extremes,
  because it uses only max|ĉ_z| and discards the rest of the spectrum: for AES it
  certifies 64 against an actual 239; for modexp, 4 against 3086. For the
  cryptographic families whose *complete* Walsh value/multiplicity distribution
  is published, S is determined **exactly** rather than bounded — a strictly
  stronger import that we flag as available but do not carry out here.

---

## 9. Demonstration of reach

Applied to modular exponentiation, the exact model yields a sharp arithmetic
criterion: PPS cost is controlled by the 2-adic structure of the multiplicative
order r. Writing r = β·2^α with β odd, cost is independent of exponent-register
width when β = 1 and Θ(2ⁿ) otherwise, with the invariance switching on exactly
at n_exp = α + 1.

That analysis is developed separately and **is not claimed here**. We note it
only as evidence that an exact cost model buys structural results that
extrapolation cannot: the criterion is a statement about arithmetic, and it is
not visible to a fitted power law.

---

## 10. Related work, and what is prior art

**The ingredient is standard and we do not claim it.** That a diagonal
operator's Pauli-Z expansion is the Walsh–Hadamard transform of its diagonal is
stated outright by Welch et al., *Efficient Quantum Circuits for Diagonal
Unitaries Without Ancillas* (arXiv:1306.3991): "the diagonals of Pauli basis
operators correspond to Walsh functions". The Pauli-spectrum ↔
Boolean-Fourier-spectrum *analogy* is likewise established (*On the Pauli
Spectrum of QAC0*, arXiv:2311.09631). Our claim is the specialisation in which
the analogy becomes an identity, and its use as a cost model.

**Adjacent, must be distinguished:**

- **Quipu / stabilizer frames** (García & Markov; arXiv:1712.03554) simulates
  reversible ripple-carry adders and QFT circuits efficiently *for specific input
  states*, using superpositions of stabilizer states. Different mechanism —
  Schrödinger-picture stabilizer frames, not Heisenberg propagation of a diagonal
  observable, and no Fourier characterisation. It establishes "quantum arithmetic
  is efficiently simulable" by another route, so our novelty must be phrased
  against it rather than around it.
- **Cîrstoiu**, *A Fourier analysis framework for approximate classical
  simulations of quantum circuits* (arXiv:2410.13856), is closest on the name and
  genuinely different — but **not** for the reason one might guess: it explicitly
  covers finite as well as compact groups. The distinction is the *domain*. Its
  circuits are C(g) = U₁(g₁)W₁ … with the U_i forming a representation of a group
  G, and the transformed function is g ↦ ⟨O⟩_{C(g)}: the expectation as a
  function of **circuit parameters**, over an **ensemble**, approximately and on
  average. Ours is y ↦ bit_j(π(y)): a function of the **input state**, for a
  **single fixed** circuit, exactly.
- **Vidal & Ballarin et al.** (NJP 2025) compute general Pauli decompositions via
  the FWHT in O(N² log N). Adjacent tooling; no diagonal special case, no
  propagation, no Boolean sparsity.
- **"Characterizing Pauli Propagation via Operator Complexity"**
  (arXiv:2510.22311) is the closest competitor in intent — it is about PPS cost —
  but uses operator stabilizer Rényi entropy and gives approximate asymptotic
  bounds, with no exact formula for any circuit class.
- **Implementations.** No available implementation propagates permutation gates
  natively. Qiskit `pauli-prop` accepts only Pauli rotation gates and rejects
  Toffoli; PauliPropagation.jl's Clifford map contains no Toffoli; stim is
  Clifford-only and cannot express one. (Checked at source; documentation was
  misleading in two cases.) See also the framework paper arXiv:2505.21606.

**On C1/C6.** That permutation circuits preserve the diagonal Pauli subalgebra
is elementary, and its ingredients are folklore in the stabilizer literature. A
targeted search found no prior statement of it *as a simulation cost mechanism*.
We present it as elementary-but-unstated rather than as a discovery.

**Context:** Gharibyan et al. (arXiv:2507.10771) is the PPS resource framework
whose boundary this work probes; their Appendix C already rules out "correlated
angles alone" as an explanation, and our claim is structural rather than angular.
Begušić & Chan (arXiv:2306.16372) for sparse Pauli dynamics at scale;
Beauregard (quant-ph/0205095) and Cuccaro et al. (quant-ph/0410184) for the two
arithmetic constructions used.

---

## 11. Limitations, and the honesty record

### 11.1 What this does not do

- **No simulation speedup.** The Walsh transform is O(2ⁿ n). This is a
  diagnostic, not an algorithm.
- **No implication for factoring.** Efficient classical simulation of Shor's
  algorithm on general inputs would be a classical factoring algorithm. Nothing
  here bears on that, and the standing constraint bounds the whole programme.
- **Scope stops at the diagonal.** X/Y observables are not covered (§3.2). Real
  Shor's measurement follows an inverse QFT and is outside this analysis.
- **Truncation is not covered by the exact statement** (§7), and behaves badly
  enough on this family to warrant the admissibility check we propose.
- **Instance sizes are modest.** Results reach 24 qubits via the Walsh route.

### 11.2 Retractions

Two substantial claims were made and withdrawn during this work; both are
recorded because the reasons are instructive.

1. **"Compilation, not algorithm, determines Pauli-path simulability."** The
   former title. Killed twice over: a propagator bug (θ = π gates treated as
   no-ops) and, independently, an ancilla-count confound in the matched
   comparison. Both compilations have exact Z-closure and the identity is
   compilation-invariant.
2. **"The QFT is the PPS bottleneck in Shor."** Killed by a bit-reversal bug
   plus a toy circuit that was not a proxy for Shor — the observable never met
   the arithmetic.

3. **"Permutation-native propagation halves peak memory exactly."** The factor
   is 2.000000 for adders but 1.9997 for modular exponentiation; the phrase was
   written from a table rounded to one decimal place. Corrected in §5, where the
   claim is now an upper bound of 2 with the measured values given. The
   `rot = 2·perm − 2` regularity in 3/3 modexp instances remains unexplained.

A further correction is methodological and worth stating: the collapse of adders
was originally explained by permutation-ness. That conclusion was right and the
mechanism wrong; the real reason is affineness (§4.2). A right conclusion via a
wrong mechanism is a failure, and is logged as one.

### 11.3 Reproducibility

Every quantitative claim in this paper carries an ID resolving to a row in
`CLAIMS.md`, which records status, evidence and location. The headline rows are
re-verified by an executable regression suite (`test_claims.py`) that runs as
part of an eight-suite correctness gate. Experiments declare their predictions
*before* measurement and carry must-fail controls; the harness reports a test in
which the control failed to fail, which has caught at least one vacuous result.

---

## 12. Conclusion

For circuits that implement a permutation of the computational basis, the number
of Pauli terms Pauli-path simulation must carry is not an empirical quantity to
be extrapolated. It is the Walsh sparsity of the Boolean function the circuit
computes — exactly, compilation-independently, and computable before the run.

The identity is elementary. What it buys is not: an exact account of why adders
collapse and general arithmetic does not, an exact factor of two from
propagating permutations natively, structural caps on achievable density with a
parity condition governing when they apply, a transfer from published
cryptanalytic constants to simulation cost with no simulation, and — developed
separately — an arithmetic criterion for modular exponentiation that a fitted
power law cannot see.

The natural next question is the one the method itself asks and we do not
answer: what happens through the inverse QFT, where the observable leaves the
diagonal and the exact model stops.
