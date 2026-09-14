# Walsh–Hadamard Sparsity Exactly Determines Pauli-Path Simulation Cost for Reversible Quantum Arithmetic

**Ian Baker**

ifkb99@gmail.com

*Working consolidation, 2026-09-10. This file is the Paper A manuscript and
tracking document: representation cost, its arithmetic applications, and its
operational limits. The new scratch-equivalence and tensor-rank results are
integrated in §6.4–6.5; state-aware contraction and conditional sampling belong
in Paper B. Detailed evidence and open tasks remain in the linked research
ledger, not in a second manuscript copy. The limitations and withdrawn-claims material in §10 is
load-bearing and should survive to submission. Claim identifiers have been moved
out of the prose into Appendix A; the working ledger `CLAIMS.md` is supplementary
material. Numbers here are post-bugfix (see §10.2).*

---

## Abstract

Pauli Path Simulation (PPS) — also called sparse Pauli dynamics or Pauli
propagation — has become a leading method for classically simulating
utility-scale quantum circuits. Its practical deployment depends on predicting,
before committing to an expensive run, how many Pauli terms a circuit will
require: a question currently answered by empirical power-law extrapolation
calibrated on brickwork and Trotterised circuits with generic rotation angles.

We show that for circuits implementing a permutation of the computational basis
— all reversible arithmetic, and hence the arithmetic core of Shor's algorithm
and similar algorithms — this quantity is not merely predictable but
*exactly characterized by a Walsh transform*, with no extrapolation and no fitting. For a
circuit implementing basis permutation π, the Heisenberg pullback π†Z_jπ is the
diagonal operator (−1)^{g(y)} with g(y) = bit j of π(y); expanding a diagonal
operator in the Pauli basis is precisely the Walsh–Hadamard transform of
(−1)^g. The Pauli support of the full pullback is therefore *exactly* the Walsh
spectrum of g, and its size is the Walsh sparsity of that Boolean function.

Four consequences follow. **(i)** PPS cost for reversible arithmetic is a
property of the **full-space basis permutation implemented** — including its
action on ancillas — and not of the gate set implementing it. We verify exact
Z-closure for both Toffoli-compiled and Fourier-compiled (Beauregard) modular
exponentiation; the supplied implementations agree on the valid arithmetic
subspace but have different scratch layouts and therefore are not a matched
same-permutation comparison.
**(ii)** The tractability of linear arithmetic is explained rather than observed
— a ripple-carry adder's low output bit is XOR-affine, Walsh sparsity 1
characterises affineness, and PPS collapses to a single term. **(iii)** Peak
memory is a distinct and larger quantity than final support; Clifford+T
decomposition can enlarge the peak, but the gap is not entirely a compilation
artifact. Propagating X, CNOT and Toffoli as atomic
permutations keeps the expansion Z-type at every step and makes the peak itself a
Walsh quantity. For the tested arithmetic rows, the measured relation is:
N_max^rot = 2·N_max^perm − |B|, where B is the set of peak-time Pauli strings
missing the target qubit of the Toffoli gadget in which the peak falls. B is
empty for ripple-carry adders, giving exactly 2; for modular exponentiation with
the standard observable it has two elements, and those two are precisely the
dominant Walsh coefficients of the pulled-back function. **(iv)** Walsh sparsity
and nonlinearity are the same object linear cryptanalysis studies, giving a
transfer: any published nonlinearity lower-bounds PPS cost for *every* circuit
computing that function, with no simulation.

The quantity the model computes is the size of the **Heisenberg representation
PPS maintains** — its memory footprint — not the difficulty of the expectation
value it is used to estimate; §1.1 states the target task and input-state regime
explicitly, since for some inputs the expectation is obtainable by other means
entirely.

These results are diagnostic, not a simulation speedup. The Walsh transform is
itself exponential, and nothing here bears on the classical hardness of
factoring. Same-layout circuits that agree on the entire clean logical code
can have different full Walsh support, while equal support counts can accompany
different measured output distributions. Moreover, Walsh support is not tensor
rank: local Walsh transforms preserve every fixed-cut singular spectrum. These
distinctions delimit the representation cost being characterized.

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
arithmetic with Hadamards and an inverse QFT. A full-operator PPS implementation
carries intermediate expansions until they cancel or are contracted with the
input. Their peak size is a memory cost of that implementation, not a lower
bound on all contraction schedules. Paper B, §10.2, distinguishes terms that
contribute to the pre-QFT work observable from terms eliminated by the actual
exponent input. Naive within-block pruning is unsafe because exponent support
is not monotone, but exact completed-control contraction is possible and changes
the width scaling of that reduced task in both order branches. Non-diagonal
observables introduced by an inverse QFT require a separate analysis.

Readers who want a sharper "so what" should read §5 (peak memory, where the model
changes what one should actually do), §7 (truncation, where it exposes a
pathology), and §8 (a bound requiring no simulation at all).

**Contribution.** The technical core (§3) applies the established identity
between diagonal Pauli expansions and classical Boolean Fourier expansions
(Montanaro and Osborne, Proposition 9) to circuit output functions. We do not
claim that identity itself as new. Our contribution is its use as an arithmetic
PPS representation-cost model and the resulting analysis: full-permutation
compilation invariance, the affine explanation of adder collapse,
permutation-native propagation with its exact suffix-Walsh characterization,
the measured arithmetic peak relations, the cryptanalytic transfer, and the
structural caps of §6. Section 9 separates these applications from prior art.

---

## 2. Preliminaries

### 2.1 Pauli propagation

We write a Pauli string as a pair of bitmasks (x, z) with

> P(x,z) = i^{|x ∧ z|} X^x Z^z,

the phase chosen so that P is Hermitian; all coefficients in a Heisenberg
expansion are then real. Every gate is compiled to a product of Pauli rotations
U = exp(−iθσ/2), σ a Pauli string, governed by a single conjugation rule:

> U† P U = P                              if [P, σ] = 0
> U† P U = cos θ · P + sin θ · (iσP)      if {P, σ} = 0.

Clifford rotations have θ = ±π/2, so cos θ = 0 and a Pauli maps to a *single*
Pauli — no branching. T gates have θ = ±π/4, giving cos = sin = 1/√2, the maximum
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

**Theorem 1.** *Let a circuit implement a permutation π of the
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

§4.1 states the first case abstractly: any two verified implementations of the
same full-space π have the same cost. The supplied Fourier and Toffoli
modular-exponentiation implementations instead agree only on the valid
subspace and have different full-space actions (for N=5, a=2, n_exp=2 and the
same measured bit, their final supports are 451 on 10 qubits and 15493 on 15
qubits). §6.3 exhibits the same distinction deliberately: two circuits
computing the same a^e mod N whose π's diverge off the valid subspace, with a
~51% cost difference. All are predictions of Theorem 1 once the full-space π
is specified.

**Verification.** Confirmed to machine precision on 6/6 instances spanning
modular exponentiation and ripple-carry addition, across both compilations, with
**maximum error ≤ 6.7 × 10⁻¹⁶**. We check not merely that the counts agree but
that the *support sets are identical* — the stronger statement, and the one that
would fail first under a coincidence.

### 3.3 Scope: exactly where this holds and where it stops

The theorem needs two things: the unitary is a basis permutation, and the
observable is diagonal.

- **Diagonal observables of any weight are covered.** Multi-qubit Z-type
  observables work identically; verified exact on 4/4 Z-type observables tested.
- **X and Y observables are not covered.** Their pullbacks leave the
  diagonal and the expansion is fully non-diagonal. This is a hard boundary, not
  a gap in the analysis, and it is where the present paper stops. Real Shor
  measures the exponent register after an inverse QFT, which is outside this
  scope; we do not address it.
- **Permutation-ness is required of the *unitary*, not of the gate set.**
  A circuit whose individual gates are not permutations may still implement one,
  and the identity applies. This is the point of §4.1.

---

## 4. Consequences

### 4.1 Compilation invariance

Because the identity depends only on π and j, PPS cost is a property of the
full-space permutation implemented, not of the gates used to implement it — the
invariance half of §3.2's statement; §6.3 is the other half. The Fourier and
Toffoli arithmetic examples below are independent applications of the identity,
not two decompositions of one common full-space π.

This contradicts a natural intuition — that Clifford+T compilation of reversible
arithmetic is distinguished, because such a circuit is a permutation matrix
*gate by gate* and so trivially preserves Z-type strings, whereas Fourier-space
arithmetic (Draper addition, as in Beauregard's 2n+3 qubit construction) is
built from controlled-phase and Z-rotation gates that are individually not
permutations. We verify that **both** supplied implementations exhibit exact
closure of the Z-type Pauli subalgebra (0 non-Z terms in each), as follows from
each complete unitary being a basis permutation; their different scratch
layouts mean they need not implement the same full-space permutation.

> An earlier version of this work claimed the opposite — that compilation
> determines simulability. It is withdrawn in full; see §10.2 for what killed
> it.

### 4.2 Affine collapse is exactly sparsity one

It is folklore that adders are easy for Heisenberg-picture methods. The identity
makes this exact: sparsity 1 characterises affine functions, so a PPS collapse
to a single term is *equivalent* to the output bit being XOR-affine.

For a 4-bit Cuccaro ripple-carry adder, the low output bit — the in-place update
b₀ ↦ a₀ ⊕ b₀ — is affine and its Walsh sparsity is **exactly 1**. The control
is the discriminating half: bit b₂ carries carries, is not affine, and has
sparsity **10**. The collapse is therefore not about permutation-ness — an
earlier explanation we retract — but about algebraic degree.

### 4.3 A cost model at O(2ⁿ n)

Computing the Walsh spectrum requires one pass to extract π and one fast
Walsh–Hadamard transform. Given an explicitly available permutation table, the
FWHT itself costs **O(2ⁿ n)**; replaying a circuit with L constant-size logical
gates to obtain that table adds **O(L·2ⁿ)**, so the end-to-end cost is
**O((L+n)2ⁿ)**. Extracting a permutation from an arbitrary compiled unitary
requires a separate representation or cost assumption. This is exponential and
no speedup; the point is that it is *cheaper than the run it predicts* — measured
**143× to 219×** on the modular-exponentiation instances of Table 1 — because it
does no branching bookkeeping and no coefficient arithmetic. It answers "will
this run fit in memory" before the run.

Two things must be stated with it. First, the model gives the **final** support
exactly; peak memory is a different and larger quantity (§5). Second, it is
exact only at δ = 0; under truncation the relationship is more subtle (§7).

### 4.4 The instances

Every claim in this paper is measured on the following set. Walsh sparsity is
computed by the identity; the peak columns come from propagation; the three
timing columns are wall-clock on one core.

**Table 1 — instances, exact costs, and wall-clock.**

| instance | qubits | final *S* | density | peak (perm) | peak (rot) | ratio | Walsh (s) | perm-PPS (s) | rot-PPS (s) |
|---|---|---|---|---|---|---|---|---|---|
| 3-bit adder, Z(b₀) | 8 | 1 | 0.0039 | 64 | 128 | 2.0000 | <0.001 | <0.01 | <0.01 |
| 4-bit adder, Z(b₀) | 10 | 1 | 0.0010 | 256 | 512 | 2.0000 | <0.001 | <0.01 | 0.01 |
| 5-bit adder, Z(b₀) | 12 | 1 | 0.0002 | 1,024 | 2,048 | 2.0000 | <0.001 | <0.01 | 0.04 |
| 4-bit adder, Z(b₂) | 10 | 10 | 0.0098 | 16 | 32 | 2.0000 | <0.001 | <0.01 | <0.01 |
| modexp N=5, a=2 | 14 | 3,086 | 0.1884 | 6,834 | 13,666 | 1.9997 | 0.011 | 0.05 | 1.57 |
| modexp N=7, a=3 | 14 | 3,206 | 0.1957 | 8,194 | 16,386 | 1.9998 | 0.011 | 0.07 | 2.24 |
| modexp N=15, a=7 | 17 | 31,176 | 0.2379 | 64,070 | 128,138 | 2.0000* | 0.215 | 1.27 | 46.98 |

\* 1.999969 — displayed to four places. **All three modular-exponentiation rows
satisfy `rot = 2·perm − 2` exactly**, while all four adder rows satisfy
`rot = 2·perm`. See §5.

Four things to read off it. **(a)** Final *S* matches the propagated term count
on every row — this is Theorem 1, and it is the support *set* that matches, not
merely its size. **(b)** The three adder rows with *S* = 1 are the affine collapse
of §4.2, with Z(b₂) = 10 as the control that must not collapse. **(c)** The
Walsh route is 143–219× faster than rotation-level propagation on the modexp
rows, and the gap widens with size. **(d)** Density stays well below ½ — the
structural cap of §6, which is why these numbers are not those of a random
function.

The largest instance reached by the Walsh route elsewhere in this work is **30
qubits** (|S| = 536,271,623 for N = 143, exactly counted); propagation stalls
near 17, which is why the table's rot-PPS column stops there. The asymmetry is
noted in §10.1 rather than dressed up: the structural results are proved and do
not depend on instance size, and the circuit series now spans a 32768-fold
range of Hilbert-space dimension.

---

## 5. Peak versus final cost, and permutation-native propagation

The final Pauli support is what Theorem 1 computes. What bounds memory in
practice is the *peak* over the propagation.

The two differ partly because of compilation. The
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

For the tested arithmetic rows, the peak relation is an exact measured
regularity. Let *c* be the target qubit of the Toffoli gadget in which the
rotation-level peak falls, and let *S* be the permutation-native peak set. The
local gadget analysis explains why a doubling is possible, but it does not by
itself determine the global peak after intervening CNOTs and interference.

First, a Z-type string commutes with every Z-rotation, so it cannot branch until
an H turns a Z into an X; the gadget's only H acts on *c*, so **only strings
carrying Z_c ever leave the diagonal.** Second, once such a string is X_c-type,
the gadget's four T gates on *c* rotate it *within* the two-dimensional space
spanned by {X_c, Y_c},

  X_c ↦ cos θ · X_c − sin θ · Y_c,   Y_c ↦ cos θ · Y_c + sin θ · X_c,

which is closed for a fixed rest-label component. But this target-local fact
does not imply that the full strings branch only once: intervening CNOTs can
change the other Pauli labels and subsequent terms can interfere. Indeed, the
four-qubit circuit
`CCX(3,1,0); CCX(1,0,2)` with observable `Z_2` has atomic peak 4 and
rotation-level peak 10, whereas the former proposed formula predicts 6 (and its
implied 2× upper bound predicts at most 8). A dense-matrix check confirms the
rotation peak, and the final Walsh and atomic coefficients agree exactly.

The arithmetic rows retain a useful empirical statement: in 9/9 tested cases
(six modular exponentiations, three ripple-carry adders), folding the
X_c/Y_c partners at the observed peak recovers the atomic peak set with the
measured multiplicities, and `N_max^rot = 2·N_max^perm − |B|`. This is a
family-specific observation, not a universal proposition. The exact general
statement that survives is that the atomic peak equals the maximum Walsh
sparsity over the corresponding circuit suffixes.

So the ratio is 2 for the adders because B is empty there, and the unexplained
"deficit of exactly 2" in modular exponentiation is the statement |B| = 2. **The
two strings are the same in every instance: Z on the measured
multiplicand-register bit alone, and that bit together with one exponent qubit —
which are precisely the two Walsh coefficients of magnitude ½** identified
independently in §6 as the dominant Fourier modes of the pulled-back function
(6/6 against the full spectrum). They are supported off the scratch register,
never acquire Z_c, and so never double.

**The constant belongs to the observable, not to the circuit family.** Holding
the circuit fixed and moving the observable changes it: measuring the low
accumulator bit gives a deficit of 4004, the reduction ancilla 4014, and a
scratch-register bit 0 — the last recovering a ratio of exactly 2. Any statement
of the form "modular exponentiation has deficit 2" is therefore a statement
about the standard computational-basis observable, and we state it that way.

This is a genuine practical recommendation regardless: for permutation circuits,
do not decompose to Clifford+T before propagating.

---

## 6. Structure in the spectrum: why arithmetic pullbacks are not generic

The identity converts PPS cost questions into Boolean-function questions, which
lets known structure do work. Reversible-arithmetic pullbacks turn out to be
markedly non-generic, and the deviations are exactly characterisable.

### 6.1 Linear structures cap the density

A **linear structure** of g is a w with g(y ⊕ w) = g(y) for all y. The relevant
classical fact is due to Carlet (*Boolean Functions for Cryptography and Coding
Theory*, **Proposition 29**): D_e f is null (resp. constant 1) iff supp(W_f) is
contained in {0,e}^⊥ (resp. in its complement). A linear-structure space of
dimension k therefore confines the Walsh support to a coset of codimension k,
capping density at 2^{−k}. Functions whose Walsh support is an affine subspace
are not, by themselves, a definition of Carlet's **partially bent** functions:
partial bentness also imposes the required flatness condition on the nonzero
Walsh spectrum. For example, OR on three variables has full (hence affine)
support but unequal coefficient magnitudes, so support geometry alone is
insufficient.

**We cite this rather than claim it.** Our contribution is locating it in
modexp pullbacks and reading it as a cost cap. Modular exponentiation carries
the linear structure w = e_msb ⊕ e_anc — the 0/1 mask supported on the
accumulator's most significant bit and the comparison ancilla: GF(2) rank n−1
in every instance, g(y ⊕ w) = g(y) verified pointwise, present in *both*
compilations, and absent from random f. Consequently density is capped at
exactly ½ — which is why measured density converges to 0.498 *from below* and
never crosses. The ripple-carry adder's linear-structure space has dimension 5,
capping its density at 2⁻⁵ (measured: 0.0098).

**Mechanism.** The structure is forced by three ingredients acting
together: flipping the accumulator's msb *is* adding 2^{m−1}, which commutes
with mod-2^m addition (verified exhaustively, 0/256 violations); the ancilla is
coupled to the msb only by XOR, so flipping both restores it; and the msb is
excluded from the swap network, so it never reaches the observed register.

### 6.2 A conditional generalisation

The above has a natural extension we did not find stated, and which we present
at corollary altitude: it follows from Proposition 29 plus the standard
decomposition of the Walsh transform over a coset partition.

**Proposition 2.** *Let a coset partition split F₂ⁿ into cells H_u, and
suppose g restricted to H_u has linear structure w_u. Then the Walsh support
avoids E = {z : w_u·z = 1 for every u}, and hence*

> density ≤ 1 − 2^{−d},   d = dim span{w_u},

*provided the system {w_u·z = 1} is consistent.*

*Proof sketch.* Decomposing over the partition gives
c_z = 2^{−n} Σ_u (−1)^{u·z_C} A_u(z′) with A_u the transform of the restriction;
Proposition 29 kills A_u(z′) whenever w_u·z′ = 1. If every A_u vanishes, so does
c_z. ∎

**Proposition 3.** *That system is consistent iff every linear dependency
among the w_u has even support. An odd dependency makes E empty and removes the
cap entirely.*

Verified on planted instances at d = 1, 1, 2, 3 with caps 0.5, 0.5, 0.75, 0.875
and **zero violations**. The statement is covariant under affine changes of
variables (Lemma 5), so it covers arbitrary subspace partitions, not only
coordinate-aligned ones — additionally confirmed by testing under a random
GF(2) change of basis. Full proofs of both propositions are in Appendix B. The
parity condition is sharp: at the same cell count, an odd dependency gives
density 0.9845 against 0.8671 for an even one.

**At scale.** The propositions are proved, so size demonstrates rather than
establishes — but the check is cheap and forecloses the obvious objection.
Repeating the decisive cases up to **n = 30 (2³⁰ = 1.07 × 10⁹ points)**:

| case | d | cap | density at n=26 | n=28 | n=30 | violations |
|---|---|---|---|---|---|---|
| d=1 (linear structure) | 1 | 0.5000 | 0.4978 | 0.4953 | 0.4908 | **0** |
| d=2, two cells | 2 | 0.7500 | 0.7455 | 0.7411 | 0.7323 | **0** |
| d=3, four cells | 3 | 0.8750 | 0.8696 | 0.8640 | 0.8531 | **0** |
| odd dependency | 2 | *none* | 0.9927 | 0.9851 | 0.9703 | **0** |
| even dependency | 3 | 0.8750 | 0.8696 | 0.8640 | 0.8531 | **0** |

**The densities drift below their caps with n, and that drift is a measurement
artifact, not structure.** With an absolute coefficient cut of 10⁻⁶ and
random-function coefficients distributed as N(0, 2^(−n/2)), the expected erased
fraction is erf(tol/(σ√2)) = 0.7%, 1.3%, 2.6% at n = 26, 28, 30. The
odd-dependency row is the clean test, since its true density is exactly 1: it
should read 0.9935 / 0.9869 / 0.9739 and reads 0.9927 / 0.9851 / 0.9703. The
decisive quantity is unaffected — **zero support elements inside E at every
size**, because E-membership is combinatorial rather than a magnitude test.

The cap is governed by the **span dimension** of the per-cell structures, not by
the depth of conditioning. This matters practically: it says which modifications
to a construction can possibly change its cost — only those adding an
independent structure vector, or introducing an odd dependency.

### 6.3 A compilation choice with a real cost

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

### 6.4 Same physical computation, different full-space support

A same-layout construction removes the qubit-count confound entirely. Let J
embed all legitimate logical inputs with clean scratch, and suppose U returns
them to a clean-output code. Any extension W that is identity on that code
satisfies WUJ=UJ, so the restricted observable and all subsequent measurement
statistics are unchanged. This is an identity of logical isometries, not just
agreement on one classical input. It need not imply U†W†OWU=U†OU off the code.

For the fixed modular-exponentiation layout, append a variable number of
Toffoli gates with both controls on clean scratch bits and target on the
observed work bit. Each is identity on the entire legitimate output code.
The experiment verifies every legitimate exponent/work basis input, coherent
code states, and full order-finding output distributions. The full Walsh
support nevertheless changes; the clean-restricted function does not. Exact
rows and validation scope are recorded in the scratch-equivalence experiment
and its ledger entry (Appendix A).

The converse is equally important. An appended CNOT from an unchanged exponent
bit to the observed work bit multiplies the sign function by that exponent's
Walsh character. It therefore translates every Walsh mask and preserves both
the support count and coefficient magnitudes. Yet it changes legitimate
arithmetic, and in the tested odd-order example also changes order-finding
output statistics. Equal counts do not certify physical equivalence; unequal
counts do not certify different physical computations. Full-operator cost is
exactly what Theorem 1 describes, and neither inference follows from it.

### 6.5 Walsh sparsity is not tensor memory

For any partition L|R of input bits, write the normalized truth tensor as a
matrix F and the coefficient tensor as C. The Walsh transform factors locally:

> C = H_L F H_Rᵀ,

with orthogonal normalized Walsh matrices. Thus every cut rank and singular
value is unchanged by this basis change. This is a standard local-basis fact,
not a new tensor theorem. A coefficient vector can be dense yet have small
tensor ranks.

For an ideal scalar function of exponent e with period r, the tensor factors
through the running residue modulo r, giving cut rank at most r. For the
natural low-k/high-(t−k) cut, the high bits access at most r/gcd(r,2^k)
residues, sharpening the bound to min(2^k, 2^(t−k), r/gcd(r,2^k)). The scalar
period can be smaller than the multiplicative order. These are upper bounds,
not claims of minimal realizations or efficient period discovery.

The controlled measurements exhibit growing ideal-function Walsh support at
constant small cut rank, but substantially larger ranks for the full
scratch-space pullback. Selected intermediate suffixes also have larger ranks
than the ideal scalar function. Detailed ranks, tolerances and controls live
in the tensor-memory experiment and its ledger entry (Appendix A). Neither
the ideal rank nor a small final rank bounds a full simulation's setup,
intermediate memory, or conditional output-sampling cost. The relevant
finite-state/tensor connection is already studied by Li, Precup and Rabusseau.

Paper B makes the distinction concrete for actual output sampling. The coherent
exponent/work Schmidt rank is min(r,2^t), even when a chosen scalar bit has rank
one. Yet an order-informed latent-eigenphase sampler needs no work vector at
all. Thus even joint-state rank is not by itself a memory lower bound for this
restricted measurement task; obtaining the spectral information has a separate
computational cost.

---

## 7. Truncation: what the model does and does not say

The identity is exact at δ = 0. Practical PPS truncates, and the model's
relationship to truncated runs is not the naive one.

- **Coefficient (δ) truncation is the right knob.** The spectrum is
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
  truncation is basis-dependent by construction.
- **Truncation error is non-monotonic in δ.** Under incremental (per-gate)
  truncation — the truncation PPS actually performs — on a Toffoli-compiled
  N = 5 modular exponentiation, δ = 10⁻¹ is *exact* while carrying a peak of 34
  terms, whereas δ = 10⁻³ is off by 0.285 while carrying 23482. Small
  coefficients cancel as a set; removing some of them is worse than removing all.
- **Truncated estimates can be inadmissible.** We observed |⟨O⟩| > 1 in
  4 of 18 runs, all at *mild* δ. We propose the operator-norm bound as a
  zero-cost admissibility check on PPS output.

The honest summary: the exact model describes the untruncated support, and
truncation on this family behaves badly enough that it should be monitored
rather than trusted.

---

## 8. The cryptanalytic bridge

Walsh sparsity and nonlinearity are the central quantities of linear
cryptanalysis. The identity turns that coincidence into a transfer.

**Proposition 4.** *For any Boolean function g with nonlinearity NL,*

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
- **Honest limitation.** The bound is loose away from the extremes,
  because it uses only max|ĉ_z| and discards the rest of the spectrum: for AES it
  certifies 64 against an actual 239; for modexp, 4 against 3086. For the
  cryptographic families whose *complete* Walsh value/multiplicity distribution
  is published, S is determined **exactly** rather than bounded — a strictly
  stronger import that we flag as available but do not carry out here.

---

## 9. Related work, and what is prior art

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
- **Georges, Berntson, Sünderhauf & Ivanov**, *Pauli decomposition via the fast
  Walsh–Hadamard transform*, New Journal of Physics **27**, 033004 (2025),
  DOI [10.1088/1367-2630/adb44d](https://doi.org/10.1088/1367-2630/adb44d),
  compute general Pauli decompositions via the FWHT in O(N² log N). Adjacent
  tooling; no diagonal special case, propagation, or Boolean-sparsity result.
- **Montanaro & Osborne**, *Quantum Boolean Functions*, Proposition 9
  ([arXiv:0810.2435](https://arxiv.org/abs/0810.2435)), explicitly identifies
  the Pauli expansion of a diagonal Boolean operator with its classical Fourier
  expansion. Applying it to a permutation pullback gives the ingredient of
  Theorem 1; our contribution is the circuit-function specialisation as a PPS
  cost account and the arithmetic applications.
- **"Characterizing Pauli Propagation via Operator Complexity"**
  (arXiv:2510.22311) is the closest competitor in intent — it is about PPS cost —
  but uses operator stabilizer Rényi entropy and gives approximate asymptotic
  bounds, with no exact formula for any circuit class.
- **Stabilizer rank, stabilizer extent, and magic monotones.** The
  Clifford+Toffoli / CNOT-dihedral simulation literature (Bravyi–Gosset and
  successors) measures a *different* cost in a *different* picture:
  Schrödinger-picture decomposition of the state into stabilizer terms, with
  cost exponential in non-Clifford count. The contrast is instructive rather
  than competitive, and sharpens our claim: **the circuits that are expensive
  for stabilizer-rank methods — Toffoli-heavy reversible arithmetic — are
  precisely the ones for which Heisenberg propagation of a diagonal observable
  is exactly characterised.** Neither cost measure bounds the other, and nothing
  in that literature computes a per-function Pauli support. Magic-monotone work
  on permutation-plus-diagonal classes is adjacent for the same reason and
  equally does not state the identity.
- **Boolean-function side.** Work computing complete Walsh spectra for
  structured families (e.g. permutation-inverse families) is the right
  neighbourhood for §6. We state Propositions 2 and 3 at corollary altitude for
  that reason: they follow from Carlet's Proposition 29 together with the
  standard coset decomposition of the Walsh transform, and we would not be
  surprised to find them folklore in that literature. We claim their application
  as a PPS cost cap, not the underlying combinatorics.
- **Implementations.** In the versions inspected for this draft, Qiskit
  `pauli-prop` accepted only Pauli rotation gates and rejected Toffoli;
  PauliPropagation.jl's Clifford map contained no Toffoli; and stim was
  Clifford-only. This is a dated account of the inspected implementations, not
  a universal absence claim. See also the framework paper arXiv:2505.21606.

**On Z-closure.** That permutation circuits preserve the diagonal Pauli subalgebra
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

## 10. Limitations, and the honesty record

### 10.1 What this does not do

- **No simulation speedup.** The Walsh transform is O(2ⁿ n). This is a
  diagnostic, not an algorithm.
- **No implication for factoring.** Efficient classical simulation of Shor's
  algorithm on general inputs would be a classical factoring algorithm. Nothing
  here bears on that, and the standing constraint bounds the whole programme.
- **Scope stops at the diagonal.** X/Y observables are not covered (§3.3). Real
  Shor's measurement follows an inverse QFT and is outside this exact diagonal
  model. Paper B includes a separately validated conditional-sampling baseline;
  it is not an extension of Theorem 1 to non-diagonal observables.
- **Truncation is not covered by the exact statement** (§7), and behaves badly
  enough on this family to warrant the admissibility check we propose.
- **Instance sizes are asymmetric, and we say so.** Circuit-level results reach
  **30 qubits** via the Walsh route (q = 15…30 across two series varying the
  modulus at fixed exponent width; density 0.4728 → 0.4994, growth exponent
  1.006 bits/qubit) and stall near 17 for propagation — bounded by gate count,
  not memory, so pushing that further is a matter of engineering rather than
  insight. The *structural* results of §6.2, being function-level, also reach
  **n = 30 (10⁹ points)**. The remaining asymmetry is between the Walsh route
  and rotation-level propagation, not between circuit and function level.

### 10.2 Claims withdrawn during this work

Three substantial claims were made and withdrawn during this work; all are
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
   written from a table rounded to one decimal place. Corrected in §5. The
   `rot = 2·perm − 2` regularity in 3/3 modexp instances was recorded here as
   unexplained through several revisions of this paper remains a measured,
   family-specific regularity. The local Toffoli argument does not prove the
   global peak relation; the four-qubit counterexample in §5 disproves that
   universal derivation. The exact atomic-peak/Walsh-suffix characterization
   remains valid.

A further correction is methodological and worth stating: the collapse of adders
was originally explained by permutation-ness. That conclusion was right and the
mechanism wrong; the real reason is affineness (§4.2). A right conclusion via a
wrong mechanism is a failure, and is logged as one.

### 10.3 Reproducibility

```
uv run python -m experiments.experiment_c8            # the identity, 6/6 instances
uv run python -m experiments.experiment_scope         # scope: Z-type yes, X/Y no; truncation
uv run python -m experiments.experiment_c17_deficit   # the peak relation of §5, 9/9
uv run python -m experiments.experiment_crypto        # the nonlinearity bound of §8
uv run python -m experiments.experiment_gf2law        # the conditional-structure law of §6.2
uv run python -m experiments.experiment_reduction2    # the compilation cost of §6.3
LAB_GPU=1 uv run python -m experiments.experiment_scratch_equivalence
OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_tensor_memory
uv run python test_claims.py                          # headline rows, pinned, ~15 s
```

Every quantitative claim in this paper resolves, via the claim map of Appendix
A, to a row in `CLAIMS.md`, which records status, evidence and location. The
headline rows are re-verified by an executable regression suite
(`test_claims.py`) that runs as part of a nine-suite correctness gate.
Experiments declare their predictions *before* measurement and carry must-fail
controls; the harness reports a test in which the control failed to fail, which
has caught two vacuous results that had already produced confident-looking
numbers.

---

## 11. Conclusion

For circuits that implement a permutation of the computational basis, the number
of Pauli terms Pauli-path simulation must carry is not an empirical quantity to
be extrapolated. It is the Walsh sparsity of the Boolean function the circuit
computes — exactly, compilation-independently, and computable before the run.

The identity is elementary. What it buys is not: an exact account of why adders
collapse and general arithmetic does not, an exact suffix-Walsh characterization
of permutation-native peaks (with decomposition gaps measured, not universally
fixed), structural caps on achievable density with a
parity condition governing when they apply, a transfer from published
cryptanalytic constants to simulation cost with no simulation, and — developed
separately — an arithmetic criterion for modular exponentiation that a fitted
power law cannot see.

The clean-code and tensor comparisons sharpen the central limitation: full
Walsh support is neither a physical-equivalence invariant nor a bound on every
compressed representation. Paper B shows how input-aware contraction changes
the consequences and validates a separate route through the inverse QFT.
Characterizing useful compressed conditional simulation remains open.

---

## References

- Beauregard, *Circuit for Shor's algorithm using 2n+3 qubits*,
  quant-ph/0205095.
- Begušić & Chan, *Fast classical simulation of evidence for the utility of
  quantum computing before fault tolerance*, arXiv:2306.16372.
- Bravyi & Gosset, *Improved classical simulation of quantum circuits
  dominated by Clifford gates*, arXiv:1601.07601.
- Carlet, *Boolean Functions for Cryptography and Coding Theory*, Cambridge
  University Press, 2021.
- Cîrstoiu, *A Fourier analysis framework for approximate classical simulations
  of quantum circuits*, arXiv:2410.13856.
- Cuccaro, Draper, Kutin & Moulton, *A new quantum ripple-carry addition
  circuit*, quant-ph/0410184.
- Draper, *Addition on a Quantum Computer*, quant-ph/0008033.
- García & Markov, *Simulation of Quantum Circuits via Stabilizer Frames*,
  arXiv:1712.03554.
- Gharibyan, Hariprakash, Mullath & Su, *A Practical Guide to using Pauli Path
  Simulators for Utility-Scale Quantum Experiments*, arXiv:2507.10771.
- Nadimpalli, Parham, Vasconcelos & Yuen, *On the Pauli Spectrum of QAC0*,
  arXiv:2311.09631.
- Rudolph, Jones, Teng, Angrisani & Holmes, *Pauli Propagation: A Computational
  Framework for Simulating Quantum Systems*, arXiv:2505.21606.
- Shao, Cheng & Liu, *Characterizing Pauli Propagation via Operator
  Complexity*, arXiv:2510.22311.
- Georges, Berntson, Sünderhauf & Ivanov, *Pauli decomposition via the fast
  Walsh–Hadamard transform*, New Journal of Physics 27, 033004 (2025),
  DOI 10.1088/1367-2630/adb44d.
- Montanaro & Osborne, *Quantum Boolean Functions*, arXiv:0810.2435.
- Li, Precup & Rabusseau, *Connecting Weighted Automata, Tensor Networks and
  Recurrent Neural Networks through Spectral Learning*,
  [arXiv:2010.10029](https://arxiv.org/abs/2010.10029).
- Welch et al., *Efficient Quantum Circuits for Diagonal Unitaries Without
  Ancillas*, arXiv:1306.3991.

---

## Appendix A — claim map

The working ledger `CLAIMS.md` records every claim's status, evidence and
location, and is supplied as supplementary material. Inline identifiers are
kept out of the prose; this table is the mapping. Identifiers prefixed
*Paper B* have their rows in the companion paper's section of the ledger; those
listed against §10.2 are **retracted** rows, cited there as the record of what
was withdrawn and not as support for anything.

| section | claims |
|---|---|
| 1.1 What is being predicted, and for which task | Paper B: C42, C45 |
| 3.1 Statement and proof | C8 |
| 3.2 Which Boolean function — a definition that must be stated precisely | C8 |
| 3.3 Scope: exactly where this holds and where it stops | C6, C13 |
| 4.1 Compilation invariance | C1, C2, C6 |
| 4.2 Affine collapse is exactly sparsity one | C10 |
| 4.3 A cost model at O(2ⁿ n) | C11 |
| 4.4 The instances (Table 1) | C8, C10, C17, C44 |
| 5. Peak versus final cost, and permutation-native propagation | C17, C44; Paper B: C18 |
| 6.1 Linear structures cap the density | Paper B: C7, C30, C31 |
| 6.2 A conditional generalisation | Paper B: C40, C41 |
| 6.3 A compilation choice with a real cost | Paper B: C32 |
| 6.4 Same physical computation, different full-space support | C50 |
| 6.5 Walsh sparsity is not tensor memory | C48; Paper B: C52 |
| 7. Truncation: what the model does and does not say | C5, C14, C16 |
| 8. The cryptanalytic bridge | C12, C25, C26 |
| 9. Related work, and what is prior art | C1, C6 |
| 10.1 What this does not do | Paper B: C7, C40 |
| 10.2 Claims withdrawn during this work | C3, C4, C9, and the F- and H-rows of the ledger's retracted section |

---

## Appendix B — proofs for §6.2

Notation as in §2.2 and §6.2. Throughout, · is the GF(2) inner product, and a
*linear structure* of a Boolean function f is a vector w with f(y ⊕ w) = f(y)
for all y (the null-derivative case of Carlet's Proposition 29).

**Lemma 5 (affine covariance of the Walsh spectrum).** *Let L be an invertible
linear map on F₂ⁿ, t ∈ F₂ⁿ, and g′(y) = g(Ly ⊕ t). Then*

> ĉ′_z = (−1)^{(L⁻¹t)·z} · ĉ_{(Lᵀ)⁻¹z}.

*In particular the Walsh support of g′ is Lᵀ applied to the support of g, and
sparsity, density, and the coefficient magnitude multiset are invariant under
affine changes of variables.*

*Proof.* Substitute x = Ly ⊕ t, so y = L⁻¹(x ⊕ t), in the transform:
ĉ′_z = 2^{−n} Σ_x (−1)^{g(x)} (−1)^{(L⁻¹(x⊕t))·z}. The exponent splits as
(L⁻¹x)·z ⊕ (L⁻¹t)·z, and (L⁻¹x)·z = x·(Lᵀ)⁻¹z. The constant factor
(−1)^{(L⁻¹t)·z} comes out of the sum and what remains is ĉ evaluated at
(Lᵀ)⁻¹z. ∎

*Proof of Proposition 2.* By Lemma 5 we may assume the coset partition is
coordinate-aligned: there is a set C of m coordinates such that the cells are
H_u = {y : y_C = u}, indexed by u ∈ F₂^m, with the remaining n − m coordinates
written y′; each w_u is supported off C, so w_u·z = w_u·z′. Write z = (z_C, z′)
and let g_u(y′) be the restriction of g to H_u. Splitting the transform over
the cells,

> ĉ_{(z_C, z′)} = 2^{−n} Σ_u (−1)^{u·z_C} Σ_{y′} (−1)^{g_u(y′)} (−1)^{y′·z′}
>              = 2^{−m} Σ_u (−1)^{u·z_C} Â_u(z′),

where Â_u is the normalised Walsh transform of g_u on F₂^{n−m}. Each g_u has
linear structure w_u, so by Carlet's Proposition 29 the support of Â_u is
contained in {z′ : w_u·z′ = 0}. If z ∈ E — that is, w_u·z′ = 1 for *every*
u — then every Â_u(z′) vanishes, hence ĉ_{(z_C, z′)} = 0 for every z_C. So the
Walsh support of g avoids E.

For the count: when the affine system {w_u·z′ = 1 for all u} is consistent, its
solution set is a coset of the subspace {z′ : w_u·z′ = 0 for all u}, which has
codimension d = dim span{w_u}. Hence |E| = 2^m · 2^{(n−m)−d} = 2^{n−d}, and

> |supp(ĝ)| ≤ 2ⁿ − 2^{n−d},   i.e.   density ≤ 1 − 2^{−d}. ∎

*Proof of Proposition 3.* Arrange the w_u as the rows of a matrix W over GF(2);
the system in question is Wz′ = 1, with 1 the all-ones vector. A linear system
over a field is solvable iff the right-hand side is orthogonal to every vector
of the left kernel: Wz′ = 1 has a solution iff v·1 = 0 for every v with
vᵀW = 0. A left-kernel vector v is exactly a linear dependency
Σ_{u ∈ T} w_u = 0 with T = supp(v), and v·1 = |T| mod 2. So the system is
consistent iff every dependency among the w_u has even support. If some
dependency has odd support, the system is inconsistent, E is empty, and the
argument above constrains nothing — the cap disappears. ∎
