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

## 2. `[x]` Weight-truncation as Fourier tail mass — DONE, negative result

**Outcome: the exact model does NOT extend to weight truncation, and weight
truncation should not be used on reversible arithmetic.** See `NOTES.md` §W.

- **W1** The naive identity is false. PPS truncates *incrementally*, so a
  discarded term never branches and the result is not the truncated final
  operator. The Fourier tail describes **terminal** truncation only. Both are
  non-monotonic in k; incremental can be 400× better (k=3) or far worse (k=8).
  This also retroactively sharpens the δ story (C14): non-monotonicity is not
  only cancellation, it is truncation changing what subsequently branches.
- **W2** Signed sums per weight level are large and alternating
  (−0.500, −0.516, +0.117, …, −0.516, +0.297), cancelling only over all levels.
  Any cutoff slices the cancellation. No usable k exists.
- **W3** By contrast δ is exactly right: the `|c|>0.1` set (4 of 3086 terms)
  sums to **exactly** ⟨O⟩ = −1.00000000 and the remaining 3082 sum to
  **exactly** 0. Clean split.
- **W4** The four dominant coefficients sit at weights **1, 2, 8, 9** — not all
  low-degree, which is precisely why weight truncation cannot substitute for δ.
  (Also corrects a misreading: the "292× low-weight enrichment" is just those two
  0.5-magnitude coefficients; `k(99% mass)=11`, so the bulk is high-weight.)

Files: `experiment_weight.py`, `experiment_weight2.py`; `perm_pps.py` gained a
`max_weight` argument.

**Subsumes step 8** (decode the dominant coefficients) — done as part of W4.

## 3. `[x]` Prove the circuit-level C15 invariance — MECHANISM SOLVED

Paper B's soft centre, now explained. Full write-up in `NOTES.md`
"C15 MECHANISM SOLVED"; files `experiment_c15_proof{,2,3,4}.py`.

**The mechanism (C23).** For β=1 and i ≥ α, every block applies the *same*
permutation V = u_a(·,1), controlled on its own qubit, and **V is an
involution**. So the circuit depends on the entire identity tail through a
single parity bit ⊕_{i≥α} e_i — one effective variable however many qubits it
spans. That explains constancy and liveness *together* (they had looked in
tension), predicts the onset n_exp = α+1, predicts the exact branch relation
h(y) = g(y ⊕ e_prev), and fails for β>1 as required.

Evidence chain, each link verified:
- branch symmetry |ĥ(z)|=|ĝ(z)|: BOTH=0 for β=1, magnitudes match bit-for-bit
  (err 0.00e+00); BOTH=14934/1036508 for β>1;
- sign pattern is a linear character, exact on 100% of support, with v a single
  bit = the previous exponent qubit;
- V² = id exhaustively (2^15 and 2^21);
- flipping any two identity-block controls leaves the function pointwise
  unchanged.

**Two failed routes, do not retry:** affineness of the identity block (14336 of
32768 violations) and the (z, z⊕e) pairing (an n_exp=1 artifact, not
r-dependent).

**FORMALISED — C15 is now a theorem for this circuit family.** Both owed links
closed once the block structure was read properly: `u_a(·,1) = A⁻¹SA` with S a
product of *disjoint* transpositions, so S²=id and V is an involution by
conjugation; the blocks then commute and compose to `V^p`; and averaging the
Walsh character over the tail confines the support to `z_I ∈ {0, 1_I}`, two
values regardless of tail length, so the size carries no dependence on it.

Its sharpest consequence was **derived before being tested** and then confirmed:
0 violations in 7/7 β=1 instances with the two halves individually constant
(7770/7779 at |I| = 2, 3, 4), against 48189 and 1556046 for the β>1 controls.

Scope limit to state in the paper: the proof uses only the
multiply–swap–unmultiply form, so it covers both compilations here and any
Vedral/Beauregard-style construction, but does not automatically transfer to a
modexp built otherwise.

## 4. `[x]` Third modulus for C15 — DONE, and it caught an over-claim

Cheap, and it paid for itself twice. See `NOTES.md` "C15 REFINED" and
"C20 STRENGTHENED"; files `experiment_c15c.py`, `experiment_c15d.py`.

- **C21 (new).** C15 as previously stated was too strong. N=5 with a=2 has r=4
  (α=2), where every earlier sweep used r=2 (α=1) — and it is *not* constant
  from n_exp=2: it grows 15493 → 32143, then freezes. **The support locks at
  n_exp = v₂(r) + 1.** Clean mechanism: the block for exponent bit i multiplies
  by a^(2^i), which is the identity on the valid subspace iff i ≥ α, so the
  first identity block appears at n_exp = α+1. Verified exactly for α=1,2;
  α=3,4 still growing at the largest reachable width, as predicted.
  **Why it was missed: every earlier sweep used α=1 and started at n_exp=2 =
  α+1, exactly on the threshold. Luck.**
- **C22 (new).** Every order divides λ(N), so λ(N) a power of two ⟹ *every* base
  is free. That happens exactly when N = 2^a × (product of distinct Fermat
  primes). For odd semiprimes: p·q with both Fermat, i.e. 15=3×5, 51=3×17,
  85=5×17, … **15 is the smallest, and all 7 of its bases are free** (vs 3/11
  for N=21). So C20 strengthens from "N=15, a=7 is degenerate" to "N=15 is
  degenerate for every base, forced by the modulus" — and the property making it
  the natural smallest demo is the same one making it uninformative.

Also raised the controlled design from two moduli to three (N = 5, 7, 21).

## 5. `[x]` Import actual cryptanalytic results — DONE, bound is real but weak

See `NOTES.md` §X; file `experiment_crypto.py`.

- **C25.** Derived and verified **S ≥ (1 − NL/2ⁿ⁻¹)⁻²** from Parseval. Converts
  any published nonlinearity into a PPS cost lower bound for *every* circuit
  computing the function, with no simulation. Tight at both extremes.
- **Validation.** AES S-box nonlinearity comes out at exactly **112** over all
  255 nonzero linear combinations — the published constant. An external check on
  the whole Walsh pipeline.
- **End-to-end.** Built reversible circuits computing the inner product (bent)
  and propagated them: support exactly 2^(2m) (64/256/1024), matching the bent
  prediction. No truncation is available when every coefficient has the same
  magnitude — the clean worst-case statement.
- **C26, the honest limit.** Loose away from the extremes: AES bound 64 vs 239
  actual, modexp 4 vs 3086. NL uses only `max|c|` and throws away the rest of the
  spectrum. Follow-up worth doing: for crypto families whose **full** Walsh
  value/multiplicity distribution is published (the AES inverse among them), S is
  determined *exactly* rather than bounded — a much stronger import.
- **C12 correction.** The "0.74 of the bent bound" figure is stable across N at
  fixed n_exp, **not** across widths (0.50 at n_exp=1). Fixed in the abstract.

## 6. `[x]` The intermediate 2-adic law — DONE, negative result

**There is no intermediate law.** See `NOTES.md` §I; files
`experiment_c15_intermediate{,2,3}.py`.

- β=1: sparsity exactly constant in t (4, 4, 16 across t=11..22). Matches the
  proof.
- β>1: sparsity is Θ(2^t), density bounded away from 0 — but it does **not**
  organise by α. With a fixed table and consecutive t it shows a strong
  period-`ord₂(β)` oscillation plus a slow upward drift, and neither component
  is a function of α.
- The periodicity hypothesis was refuted too: residue classes drift rather than
  staying constant (r=7, t≡2 mod 3: 0.834, 0.831, 0.850, 0.875).
- **The N=323 "intermediate" observation was sampling aliasing** — t=16,20,24
  hit residues 4,2,0 mod ord₂(9)=6, so three points of an oscillation were read
  as a trend.

Net: the binary β=1/β>1 dichotomy is the robust structure and C15 should be
stated without hedging about intermediate regimes. A negative result that
protects the paper rather than complicating it.

**Two methodology bugs of my own, caught by the protocol:** degenerate random
tables at small r (rejected constants), and drawing a new table per t, which
confounded t-dependence with table variance — a direct violation of "vary
exactly one parameter". Both produced confident-looking numbers first.

**Scope caveat:** used `g(e)=h[e mod r]` with random h to sweep r freely, not
the real modexp bit function. The β=1 conclusion is unaffected (proved
independently); fine structure within β>1 could differ for the real function.

## 7. `[x]` Does C15 survive truncation? — DONE, mostly yes

See `NOTES.md` §T; files `experiment_c15_trunc.py`, `experiment_c15_trunc2.py`.
Predictions were written before measuring.

- **P1 confirmed exactly.** Terminal thresholding preserves the constancy at
  *every* δ (counts identical across n_exp=3,4,5 at all seven δ values, both
  β=1 moduli); the β=3 control diverges. Follows from C24/M2 — identical
  magnitude multisets, so any magnitude threshold keeps identical counts.
- **P2, the practically relevant half: peak cost survives.** Under incremental
  truncation (what PPS actually does) `N_max` is SAME at every δ tested
  (24369 / 24369 / 1028 / 34). That is the quantity bounding memory.
- **C28, the caveat.** `N_final` drifts at aggressive δ (8 → 16 → 32) and ⟨O⟩
  collapses to 0 at δ=1e-1 for the wider circuits while n_exp=3 stays exact.
  The wider circuit is *more fragile at the same δ* despite an identical exact
  spectrum, because it has more gates and so more incremental truncation events.

Net: "period-finding precision is free" holds for memory at all δ, and for
accuracy up to moderate δ. Must not be stated unqualified.

## 8. `[x]` Decode the dominant coefficients — DONE inside step 2

For modexp N=5 the four dominant coefficients are `{x8,e13}` and `{x8}` at
|c| = 0.5 (weights 2 and 1), plus a pair at |c| = 0.1465 (weights 8 and 9). They
sum to exactly ⟨O⟩; the other 3082 sum to exactly 0. See `NOTES.md` §W3–W4.

Still open: the two large coefficients are supported almost entirely on the
x-register bit being measured plus one exponent qubit. Whether that generalises
across N, a and observable is untested and would be cheap to check.

---

## 9. `[ ]` Extend the C15 proof to ANY modexp construction

The proof of C15/C23 uses only that the controlled-multiplier block has the
**multiply–swap–unmultiply** form: `u_a(ctrl,a) = M(a) ; SWAP ; M(a⁻¹)⁻¹`. At
a=1 that collapses to `V = A⁻¹SA`, a conjugate of an involution, and everything
follows. This covers both compilations studied here and any
Vedral/Beauregard-style construction — but *not* automatically a modular
exponentiation built some other way.

Question: is the invariance a property of the algorithm or of this circuit
family? Routes:
- Survey other modexp constructions (Zalka, Takahashi–Kunihiro, Gidney's
  windowed/lookup-based arithmetic, Häner–Roetteler–Svore) and check which have
  a block that is an involution at a=1. Windowed arithmetic in particular does
  *not* obviously have the multiply–swap–unmultiply shape.
- Find the weakest sufficient condition. Conjecture: it suffices that the block
  at a=1 be **any** involution on the full space, since (i)–(iv) never use the
  internal structure of V beyond V²=id. If so, the theorem generalises to every
  construction whose a=1 block is self-inverse, which is a much larger class and
  a cleanly checkable criterion.
- Failing that, find a construction where the invariance genuinely **fails**.
  That would be as valuable as generalising it — it would show the effect is
  compilation-dependent, and would need saying loudly in Paper B, whose current
  framing implies it is a property of the arithmetic.

Note the conjecture above is cheap to test: take an existing block, replace it
with a hand-built non-involutive permutation that is still identity on the valid
subspace, and see whether constancy breaks.

## Housekeeping

- `[x]` Yao.jl backend — closed as far as possible without installing it. Its
  docs list Toffoli under **"Clifford Gates: Two-qubit gates"**, wrong on both
  counts (Toffoli is neither two-qubit nor Clifford), so the doc is unreliable;
  and PauliPropagation.jl, which it most plausibly wraps, has **zero** mentions
  of Toffoli/CCX/CCZ and no permutation gate type. No library documents or
  exploits diagonal closure. A definitive check would need Yao installed.
- `[ ]` Decide whether `.env` should be tracked. It holds only an `LD_PRELOAD`
  path (no secret) but hardcodes an absolute path specific to this machine;
  currently committed.
- `[ ]` `experiment.py` §3, `experiment2.py` §§2–3 and all of `experiment3.py`
  are retracted. They carry warning headers; consider deleting once the papers
  are drafted and nothing references them.
