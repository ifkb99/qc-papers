# TODO — ranked

Ranked by expected value, not by order of discovery. Rationale included so a
future session can re-rank rather than follow blindly.

Status key: `[ ]` not started · `[~]` in progress · `[x]` done · `[!]` blocked

---

## 1. `[x]` Prior-art sweep on C17 and C1/C6  — DONE, result favourable

**Why first:** cheap (hours), decisive, and could have invalidated a claim before
any writing effort went into it. Produces no new science; pure de-risking.

**Result: no prior art found for the claims. Two adjacent works to cite.**

Checked at source level (not just docs — docs were misleading in both cases):

- **Qiskit `pauli-prop` 0.2.0** — supports *only* Pauli rotation gates
  (`rx/rxx`, `ry/ryy`, `rz/rzz`, `PauliEvolutionGate`) plus Pauli-Lindblad
  noise. Feeding it a Toffoli raises. **No permutation-gate support at all.**
  (A web summary claimed Toffoli support; reading `circuit_to_rotation_gates`
  shows otherwise. Do not trust the summary.)
- **PauliPropagation.jl v0.7.3** — `clifford_map` is
  `:H :X :Y :Z :SX :SY :S :CNOT :CZ :ZZpihalf :SWAP`. **No Toffoli.**
- **Yao.jl Pauli propagation backend** — docs list "Toffoli" but under
  *two-qubit* gates, which is wrong on its face; likely a doc error and it
  probably wraps PauliPropagation.jl. **Worth a source-level check before
  submission** — this is the one loose end.
- **stim** — Clifford-only (`Tableau`, `PauliString`, `TableauSimulator`).
  Cannot express Toffoli. Ruled out.

**Adjacent work that must be cited and distinguished:**

- **Quipu / stabilizer frames** (García & Markov, Michigan). Simulates reversible
  ripple-carry adders and QFT circuits in polynomial time/space *for specific
  input states*, using superpositions of stabilizer states.
  **Different mechanism** — Schrödinger-picture stabilizer frames, not
  Heisenberg propagation of a diagonal observable, and no Walsh/Fourier
  characterisation. But it establishes "quantum arithmetic circuits are
  efficiently simulable" by another route, so the novelty claim must be phrased
  against it, not around it.
  - [Hybrid Techniques for Simulating Quantum Circuits using the Heisenberg
    Representation](https://deepblue.lib.umich.edu/handle/2027.42/107198)
  - [Simulation of Quantum Circuits via Stabilizer
    Frames](https://arxiv.org/pdf/1712.03554)
- Earlier sweep (already logged in `ABSTRACT.md`): the NJP FWHT
  Pauli-decomposition paper and "Characterizing Pauli Propagation via Operator
  Complexity" are both adjacent, neither states the claim.

**Concept-level hits — the ingredient is more thoroughly standard than assumed.
Cite all three; do not let a referee find them first.**

- **"Efficient Quantum Circuits for Diagonal Unitaries Without Ancillas"**
  ([arXiv:1306.3991](https://arxiv.org/pdf/1306.3991)) — states outright that
  *"the diagonals of Pauli basis operators correspond to Walsh functions"* and
  builds diagonal-unitary synthesis on it. **This is the citation for the
  standing caveat.** The diagonal ↔ Walsh correspondence is textbook; our
  contribution cannot include it.
- **"On the Pauli Spectrum of QAC0"**
  ([arXiv:2311.09631](https://arxiv.org/pdf/2311.09631)) — the Pauli spectrum as
  the natural quantum generalisation of the Fourier spectrum of a Boolean
  function. The *analogy* is established; ours is the specialisation to
  permutation circuits where the analogy becomes an identity.
- **Cîrstoiu, "A Fourier analysis framework for approximate classical
  simulations of quantum circuits"**
  ([arXiv:2410.13856](https://arxiv.org/pdf/2410.13856)) — closest on the name,
  but checked and genuinely different: harmonic analysis over **U(1)/U(4)**
  (continuous rotation parameters), **not GF(2)ⁿ**; gives **approximate,
  truncated** estimates of expectation values with asymptotic guarantees, not
  exact Pauli-support counts. Must be distinguished explicitly — the titles are
  close enough that a referee will ask.

**Net: the claims survive, with narrowed framing.**

- Not found anywhere: the identity used as an **exact PPS cost model** for
  permutation circuits (C8); **permutation-native propagation** and its
  diagonal-closure consequence (C17); the **2-adic order** result (Paper B).
- Already standard and must be cited, not claimed: diagonal ↔ Walsh
  (1306.3991); Pauli-spectrum ↔ Fourier-spectrum analogy (2311.09631);
  efficient simulation of reversible arithmetic by other means (Quipu).
- The correct one-line positioning: *the ingredients are known; what is new is
  recognising that for permutation circuits they compose into an exact cost
  model for Pauli propagation, and what follows from that.*

---

## 2. `[ ]` Weight-truncation as Fourier tail mass

**The biggest unexplored idea.** The Pauli weight of `Z^z` is `popcount(z)`,
which *is* the Fourier degree of that coefficient. So for permutation circuits,
**weight-truncation error is literally the Fourier tail mass** of the pulled-back
bit function.

Current work covers only coefficient-truncation (δ). Weight-truncation is the
other standard PPS knob and is completely untouched — this roughly doubles the
model's reach.

Concrete: measure mass-by-weight profiles for adder vs modexp vs the binomial
profile of a random function; check whether the exact identity extends to
weight-truncated PPS. Effort: days. Logged as R3 in `experiment_review.py`.

## 3. `[ ]` Prove the circuit-level C15 invariance

Paper B's soft centre. The valid-subspace argument **demonstrably does not
cover it** — the added exponent qubits are live (each set in ≈7779/15549 support
terms), so the support is relabelled, not confined.

Most promising route: is the support confined to an *affine subspace*? That
would also explain density → ½ exactly. Until this lands, Paper B reports an
empirical regularity over five widths and two moduli — publishable but weaker
than it reads.

## 4. `[ ]` Third modulus for C15

N=5, a=2 (r=4) was never swept across n_exp. Cheap insurance moving the
controlled design from two moduli to three. Hours, not days.

## 5. `[ ]` Import actual cryptanalytic results (C12 follow-through)

Feed functions with *published* nonlinearity bounds — multiplicative inverse /
AES S-box core, bent functions — through as reversible circuits. Turns "modexp
sits at 0.74 of the bent bound" from an observation into a genuine transfer of a
large existing literature. Highest upside, most speculative.

Also worth building: a **degree ladder** (affine → quadratic → dense) of
functions with known algebraic degree, as a controlled knob validating the
sparsity ↔ nonlinearity relation. The adder bits already hint at it
(b0:1, b1:4, b2:10, b4:46).

## 6. `[ ]` The intermediate 2-adic law

N=323 (r=144=16·9) shows partial sparsity at small t (density 0.981 at t=16)
washing out to 1.000 by t=24. Is there a quantitative law in α versus t, rather
than the current binary power-of-two / not split?

## 7. `[ ]` Does C15 survive truncation?

Every C15 figure is δ=0. Given C14's non-monotonicity, the practical claim needs
checking at realistic δ before it can be offered as guidance.

## 8. `[ ]` Decode the dominant coefficients

Which qubits/registers carry the 4 of 15493 Walsh coefficients that reproduce
⟨O⟩ exactly (F13)? Would turn the aggressive-truncation rule from empirical to
structural. Logged as R4 in `experiment_review.py`.

---

## Housekeeping

- `[ ]` Source-level check of Yao.jl's Pauli propagation backend (loose end from
  step 1).
- `[ ]` Decide whether `.env` should be tracked. It holds only an `LD_PRELOAD`
  path (no secret) but hardcodes an absolute path specific to this machine;
  currently committed.
- `[ ]` `experiment.py` §3, `experiment2.py` §§2–3 and all of `experiment3.py`
  are retracted. They carry warning headers; consider deleting once the papers
  are drafted and nothing references them.
