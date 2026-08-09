---
id: 1
state: done
title: Prior-art sweep on C17 and C1/C6
outcome: DONE, result favourable
claims: [C8, C17]
---
# Prior-art sweep on C17 and C1/C6

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
  ([arXiv:2410.13856](https://arxiv.org/pdf/2410.13856)) — closest on the name.
  Still genuinely different, but **the reason first given here was wrong and is
  corrected** (2026-08-08, body read): it does *not* require continuous
  parameters — §II explicitly covers "compact **or finite** groups (or
  homogeneous spaces)". The real distinction is **what is being transformed**:
  their circuits are `C(g) = U₁(g₁)W₁ … U_D(g_D)W_D` with the `U_i` forming a
  representation of a group G, and the function Fourier-analysed is
  `g ↦ ⟨O⟩_{C(g)}` — the expectation value as a function of **circuit
  parameters**, over an **ensemble**. Ours is `y ↦ bit_j(perm(y))`, a function
  of the **input basis state**, for a **single fixed** circuit. Their results
  are approximate, average-case, mean-square-convergent, and often lean on noise
  to supply a spectral gap; ours are exact and worst-case. Different domain,
  different object, different regime. Must still be distinguished explicitly —
  the titles are close enough that a referee will ask — but distinguish it on
  *domain*, not on discreteness.

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
