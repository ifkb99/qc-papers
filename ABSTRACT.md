# Draft abstract — Paper A (cost model)

**Working title:** *Walsh–Hadamard Sparsity Exactly Determines Pauli-Path
Simulation Cost for Reversible Quantum Arithmetic*

Baker, I.

> **Split note.** The modular-exponentiation / 2-adic results (C15, C18, C7,
> F12, C19, C20) now live in **`ABSTRACT_SHOR_2ADIC.md` ("Paper B")**. They use
> this paper's identity as an instrument but make a separate claim to a separate
> audience. Their ledger rows sit in Paper B's section of `CLAIMS.md`.

---

> **STATUS — read before circulating.** Revised after a bug in the propagator
> (θ=π gates treated as no-ops) invalidated the previous thesis. The former
> title, *"Compilation, Not Algorithm, Determines Pauli-Path Simulability"*, is
> **retracted**: the apparent compilation effect was an ancilla-count confound,
> and both compilations tested have exact Z-closure. The current claim (Walsh
> sparsity = PPS term count) is verified to machine precision on 6/6 instances
> across both compilations, with supports identical rather than merely counts
> matching. C7 (pre-asymptotic) is since **refuted** — the Walsh route reached
> 24 qubits and density converges to ½ with slope 1.008 bits/qubit. See
> `NOTES.md` STATUS 2 and the honesty log.

---

## Abstract (draft, v2)

Pauli Path Simulation (PPS), also called sparse Pauli dynamics, has become a
leading method for classically simulating utility-scale quantum circuits. Its
practical deployment depends on predicting, before committing to an expensive
run, how many Pauli terms a given circuit will require — a question currently
answered by empirical power-law extrapolation calibrated on brickwork and
Trotterised circuits with generic rotation angles.

We show that for the class of circuits implementing a permutation of the
computational basis — that is, all reversible arithmetic, and hence the bulk of
Shor's algorithm and of quantum algorithms generally — this quantity is not
merely predictable but *exactly computable in closed form*, with no
extrapolation and no fitting.

For a circuit implementing basis permutation π, the Heisenberg-picture pullback
of a computational-basis observable, π†Z_jπ, is the diagonal operator
(−1)^{g(y)} where g(y) is bit j of π(y). Expanding a diagonal operator in the
Pauli basis is precisely the Walsh–Hadamard transform of (−1)^g. It follows that
the Pauli support PPS must carry is exactly the **Walsh spectrum** of g, and its
cardinality is the **Walsh sparsity** of that Boolean function. We confirm this
numerically to machine precision across modular exponentiation and ripple-carry
addition, recovering not only the term count but the identical support set.

Three consequences follow. First, PPS cost for reversible arithmetic is a
property of the Boolean function being computed, not of the gate set used to
compute it: we verify that Toffoli-compiled and Fourier-compiled (Beauregard)
modular exponentiation both exhibit exact closure of the Z-type Pauli
subalgebra, contradicting the intuition that Clifford+T compilation is
distinguished. Second, the special tractability of linear arithmetic is
explained exactly rather than empirically — a ripple-carry adder's low output
bit is XOR-affine, its Walsh spectrum is a single point, and PPS collapses to
one term, since sparsity one characterises affineness. Third, the cost model is
computable in O(2ⁿn) time, three orders of magnitude faster than the PPS run it
predicts in our instances.

We distinguish two cost measures that the identity treats differently. The Walsh
sparsity gives the *final* Pauli support exactly; peak memory during propagation
is a separate, larger quantity, because the standard Clifford+T decomposition of
a Toffoli passes through Hadamard gates and so leaves the diagonal mid-circuit
even though the gate is a permutation. We show this is an artifact of
compilation rather than of the method: propagating with X, CNOT and Toffoli as
atomic permutation primitives keeps the expansion Z-type at every step, halves
peak memory exactly, and makes the peak itself a Walsh quantity — the maximum
over circuit suffixes of the corresponding sparsity — while running an order of
magnitude faster than rotation-level propagation on the same circuits.

As a demonstration of the model's reach we apply it to modular exponentiation,
where it yields a sharp arithmetic criterion: Pauli-path cost is controlled by
the 2-adic structure of the multiplicative order, being independent of
exponent-register width when the order is a power of two and Θ(2^n) otherwise.
That analysis is developed separately (Paper B, `ABSTRACT_SHOR_2ADIC.md`) and is
not claimed here; we note it only as evidence that an exact cost model buys
structural results that extrapolation-based estimates cannot.

Walsh sparsity and nonlinearity are the central quantities of linear
cryptanalysis, and the identity turns that coincidence into a transfer. Since
nonlinearity satisfies NL = 2^(n-1)(1 - max_z |c_z|) and Parseval forces
sum_z c_z^2 = 1, the support obeys S >= (1 - NL/2^(n-1))^(-2): any published
nonlinearity is a lower bound on Pauli-path cost for *every* circuit computing
that function, obtained without simulation. The bound is attained exactly at
both extremes — affine functions (support one) and bent functions (support 2^n,
flat spectrum) — and we verify the latter end to end, propagating through
reversible circuits computing an inner product and recovering full support with
no truncation available. Applied to the AES S-box, whose nonlinearity of 112 at
n = 8 we reproduce independently, the bound certifies at least 64 Pauli terms
for any circuit computing an output bit. We note the bound is loose away from
the extremes, since it uses only the largest coefficient; for the several
cryptographic families whose complete Walsh spectrum is published, the support
is determined exactly rather than bounded.

These results are diagnostic, not a simulation speedup: the Walsh transform is
itself exponential, and nothing here bears on the classical hardness of
factoring, since efficient simulation of Shor's algorithm on general inputs
would constitute a classical factoring algorithm.

---

## Superseded abstract (v1 — retracted, kept for the record)

Pauli Path Simulation (PPS), also known as sparse Pauli dynamics, has become a
leading method for classically simulating quantum circuits at utility scale,
notably in reproducing IBM's 127-qubit kicked-Ising experiments. Its cost model
is well understood for the circuit family on which it was developed: brickwork
and Trotterised Hamiltonian dynamics with generic rotation angles, where the
distribution of Pauli coefficients follows a power law and a truncation
threshold δ trades accuracy against the number of retained Pauli terms in a
predictable way.

We ask whether that cost model transfers to *structured arithmetic* circuits —
modular addition, multiplication and exponentiation — which dominate the gate
count of Shor's algorithm and of quantum algorithms more broadly, and which
differ from the brickwork family on every relevant axis: highly correlated
rotation angles, long-range connectivity, and deterministic rather than random
structure.

Our central observation is that PPS cost for arithmetic is not a property of the
arithmetic being computed, but of how it is compiled. A classical reversible
circuit built from Toffoli, CNOT and X gates is a permutation matrix *gate by
gate*; conjugating a diagonal operator by a permutation yields a diagonal
operator, so Z-type Pauli strings are closed under back-propagation and the
exponentially many intermediate branches induced by non-Clifford T gates cancel
exactly. The same logical map compiled to Fourier-space arithmetic (Draper
addition, as used in Beauregard's 2n+3 qubit construction) is built from
controlled-phase and Z-rotation gates which are *not* individually permutations;
the closure is destroyed and the Pauli support grows by four orders of magnitude
on identical instances. We quantify this gap across matched circuit pairs and
identify the structural invariant — closure of the Z-type subalgebra under
gate-level conjugation — that predicts which side of it a given compilation
falls on.

We further report that coefficient truncation on these circuits produces
expectation estimates that violate the operator norm bound |⟨O⟩| ≤ 1, and are
therefore not merely inaccurate but inadmissible, indicating that truncation
error is uncontrolled in this regime rather than small. We propose the norm
bound as a zero-cost admissibility check for PPS output.

These results are negative for simulation capability and carry no implication
for the classical hardness of factoring: efficient simulation of Shor's
algorithm on general inputs would constitute a classical factoring algorithm,
and nothing here bears on that. Their content is instead diagnostic. They
suggest that the applicability of PPS should be assessed against the compiled
circuit rather than the algorithm, that arithmetic-heavy circuits should be
compiled to permutation form before PPS is attempted, and that resource
extrapolation procedures validated on brickwork circuits should not be applied
to structured arithmetic without re-validation.

---

## Claims ledger

The claims ledger lives in **`CLAIMS.md`**, which is the single source of truth
for claim statuses; each row there must be independently reproducible before
this abstract goes anywhere, and all figures are post-bugfix (anything citing
pre-fix numbers is void). The claims this paper rests on are **C8, C1, C2, C6,
C10, C11, C17, C12, C25, C26, C16, C13, C14, C5**, plus **C7, C15, C18, C30,
C31, C32** whose rows sit in Paper B's section of that file. The dead ones
retired here — **C3, C4, C9**, and the F/H findings that preceded them — are in
`CLAIMS.md` under "Retracted / dead".

### Prior-art check on C8 — DONE, result is favourable but qualified

Searched and read. Nothing found that states the claim. What exists nearby:

- **Vidal & Ballarin et al., "Pauli decomposition via the fast Walsh–Hadamard
  transform"** ([NJP 2025](https://iopscience.iop.org/article/10.1088/1367-2630/adb44d)).
  General 2ⁿ×2ⁿ matrix → Pauli coefficients in O(N² log N) via FWHT. Explicitly
  claims novelty for its formulas. But: **no diagonal special case, no Pauli
  propagation, no reversible circuits, no Boolean sparsity.** Adjacent tooling,
  different question.
- **"Characterizing Pauli Propagation via Operator Complexity"**
  ([arXiv:2510.22311](https://arxiv.org/html/2510.22311)). The closest
  competitor — it *is* about PPS cost. But it uses **operator stabilizer Rényi
  entropy**, gives **approximate bounds** (one O(s²) asymptotic for the 1D XY
  model), and touches none of Walsh/permutation/diagonal-pullback. No exact
  formula for any circuit class.
- **PauliPropagation.jl v0.7.3** (read the source): Clifford map is
  `:H :X :Y :Z :SX :SY :S :CNOT :CZ :ZZpihalf :SWAP` — **no Toffoli**, no
  permutation or Walsh handling. Not implemented there.
- **stim**: Clifford-only API (`Tableau`, `PauliString`, `TableauSimulator`).
  Cannot propagate Toffoli at all, so it cannot express this. Ruled out.

**Honest caveat that must go in the paper.** The *ingredient* — a diagonal
operator's Pauli-Z expansion is the WHT of its diagonal — is standard, and now
has a specific citation: Welch et al., *Efficient Quantum Circuits for Diagonal
Unitaries Without Ancillas* ([arXiv:1306.3991](https://arxiv.org/pdf/1306.3991)),
states that "the diagonals of Pauli basis operators correspond to Walsh
functions". The Pauli-spectrum ↔ Boolean-Fourier-spectrum *analogy* is also
established (*On the Pauli Spectrum of QAC0*,
[arXiv:2311.09631](https://arxiv.org/pdf/2311.09631)).

The claimed contribution is therefore narrower and must be stated as: for
permutation circuits with computational-basis observables the analogy becomes an
**exact identity**, which turns it into a **cost model for Pauli propagation** —
plus what follows (compilation invariance, the affine/sparsity-1 explanation of
adder collapse, permutation-native propagation, the cryptanalysis bridge).
Frame it that way or a referee will, less kindly.

**Second prior-art round (source-level, not docs).** No implementation found
that propagates permutation gates natively: Qiskit `pauli-prop` accepts *only*
Pauli rotation gates and rejects Toffoli; PauliPropagation.jl's `clifford_map`
has no Toffoli; stim is Clifford-only. Yao.jl's docs list Toffoli but under
"two-qubit gates", which is wrong on its face — **source-level check still
owed**. Must also cite and distinguish **Quipu / stabilizer frames**
([arXiv:1712.03554](https://arxiv.org/pdf/1712.03554)), which simulates
reversible ripple-carry adders efficiently by a different mechanism
(Schrödinger-picture stabilizer superpositions, specific input states, no
Fourier characterisation), and **Cîrstoiu**
([arXiv:2410.13856](https://arxiv.org/pdf/2410.13856)), whose harmonic analysis
is over U(1)/U(4) with approximate truncated guarantees rather than over GF(2)ⁿ
with exact support counts.

### What to do next, in order

This list is stale: the C12 crypto import was completed (TODO step 5, see the
X section of `NOTES.md` and C25/C26 in `CLAIMS.md`); the C3 redo was mooted by C8 and
never done; and C7's bit-packed core became unnecessary once the Walsh route
reached 24 qubits. `TODO.md` is the live ranking of remaining work.

## Related work to position against

- Gharibyan et al., [arXiv:2507.10771](https://arxiv.org/pdf/2507.10771) — the
  PPS resource framework this work probes the boundary of. Their App. C already
  rules out "correlated angles alone" as the explanation; our claim is
  structural, not angular, and must be distinguished from theirs explicitly.
- Begušić & Chan, [arXiv:2306.16372](https://arxiv.org/pdf/2306.16372) — sparse
  Pauli dynamics vs IBM 127-qubit.
- Rudolph et al. — Pauli propagation foundations and the Julia implementation;
  need to check whether the Z-closure observation is already folded in there.
- Beauregard, [quant-ph/0205095](https://arxiv.org/abs/quant-ph/0205095) — the
  Fourier modexp construction used here.
- Cuccaro et al., [quant-ph/0410184](https://arxiv.org/abs/quant-ph/0410184) —
  the ripple-carry adder used for the Toffoli side.
- Dang, Hill & Hollenberg, [arXiv:1712.07311](https://arxiv.org/abs/1712.07311) —
  MPS simulation of Shor; the natural comparison point for a different method
  meeting the same circuits.

**Prior-art check still owed:** C1/C6 are elementary enough that they are
plausibly folklore in the stabilizer-simulation literature. Before any
submission, search specifically for prior statements that permutation circuits
preserve the diagonal Pauli subalgebra in a Heisenberg-picture simulator.
