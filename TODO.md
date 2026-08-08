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

**Scope caveat CLOSED** (`experiment_c15_realfn.py`). Rechecked with the real
table `h[c]=bit_j(a^c mod N)`: P1 holds (constant sparsity for β=1), and α is
confirmed not to organise the data — α=0 (r=3) matches α=1 (r=6) exactly, while
α=1 spans r=6,10,18,22 with wildly different behaviour. Step 6's conclusion
transfers.

**New observation, unpursued:** the real modexp table is measurably non-generic
— sparser than random at r=6 (exactly 0.500 on even t, a clean factor of two)
and denser at r=10. So modexp bit functions carry structure beyond "depends on
e mod r". Random-table densities must not be used as a proxy for real ones.

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

## 9. `[x]` Extend the C15 proof beyond multiply–swap–unmultiply — DONE

**Confirmed: V² = id is the entire condition.** See `NOTES.md` §G; file
`experiment_c15_general2.py`.

Steps (i)–(iv) of the proof never use V's internals beyond V²=id — "identity on
the valid subspace" was context, not an ingredient. Tested with synthetic blocks
unrelated to modular arithmetic:

- controlled swap `x0↔b0` (order 2): support **constant** at 15362 across
  k=1..4 (and 15248 at N=5);
- controlled 3-cycle `x0→b0→b1` (order 3): support **grows** 15362 → 30984 →
  62088 → 123838;
- vacuity check passes (k=0 gives 3116 vs k=1's 15362), so the blocks act.

Generalised statement: *if a circuit contains k blocks each the same permutation
V controlled on its own qubit, V not modifying those controls, and V²=id, then
the Walsh support of any computational-basis pullback is independent of k.*

Paper B's scope caveat is replaced by a **criterion**: the theorem covers any
modexp construction whose a=1 block is an involution — checkable per
construction. Open follow-up: does windowed / table-lookup arithmetic
(Gidney-style) qualify? Now a well-posed question rather than a survey.

**Protocol catch:** pass 1 was vacuous — blocks acted only on b-qubits while the
observable was Z_x0, so both block types came out constant for a trivial reason.
The tell was the must-fail control failing to fail. Without it, a vacuous test
would have "confirmed" the conjecture for the wrong reason.

## 10. `[x]` Does a different modular reduction lift the ceiling? — YES, ~51%

See `NOTES.md` §L2; files `experiment_reduction.py`, `experiment_reduction2.py`.

**First hypothesis refuted.** The msb↔anc CNOT pairing is *not* the cause —
conjugating the ancilla with extra CNOT/Toffoli couplings left defect 1 and the
same w in every variant. Adding XOR couplings cannot break an XOR symmetry.

**Real mechanism (C31), three ingredients:** (a) flipping the msb *is* adding
2^(m−1), which commutes with mod-2^m addition — verified exhaustively, 0/256
violations; (b) anc is coupled to msb only by XOR, so flipping both restores it;
(c) the msb is `b[n]` and the cswaps use `b[:n]`, so it never reaches x.

**Breaking it needs nonlinearity in the msb (C32):**

```
   v0 baseline                 density 0.473  defect 1
   v6 cswap(t0, msb, anc)      density 0.473  defect 1   <- linear, survives
   v4 toffoli(msb, t0, anc)    density 0.716  defect 0   <- BROKEN
   v5 toffoli(msb, t0, t1)     density 0.721  defect 0   <- BROKEN
```

All compute `a^e mod N` correctly. **~51% cost increase from a compilation
choice that changes nothing about the computed function** — the first and only
such effect in the project, everything else having been compilation-invariant.
Constant-factor only; Θ(2ⁿ) is unaffected.

Framing worth keeping: not "compile to X to make PPS cheaper", but "the standard
construction is already cheaper than it needs to be, for a reason nobody
designed".

**Left open:** broken variants sit at 0.716–0.721, not 1.000 (random reaches
1.000), so residual *non-linear* structure remains after the linear structure is
destroyed. Unidentified — a natural next thread.

## 11. `[x]` Identify the residual non-linear structure — SOLVED (C33, C34, C35)

**The residue is a conditional linear structure with an exact ¾ density cap.**
See `NOTES.md` §AF and §RS; files `experiment_affstruct.py`,
`experiment_resid1.py`, `experiment_resid2.py`. Density story now complete:
½ = linear structure (C30), ¾ = conditional structure (C34), no forced cap =
two independent nonlinear monomials, 1.000 = structureless.

How the planned sequence resolved:

1. **Affine structures** — counting lemma: density > ½ forbids linear AND
   affine structures, so the circuit half was settled without measurement.
   At function level the finder explained §I's exact-0.500 rows outright:
   all-ones structures (r=6 AFFINE, r=3 LINEAR, r=10 negative control clean),
   PROVED with a three-ingredient criterion (C33).
2. **Support complement** — it IS a recognisable set: the quadrant
   {z_msb=1, z_anc=0} is exactly empty, 0 violations in 3/3 instances (C34).
   Mechanism: t is restored between reductions, so the wraps are inert on the
   t0=0 half-space (baseline structure survives conditionally, proved) and
   degenerate to CNOTs on the t0=1 half (the structure rotates, verified,
   rather than breaks).
3. **Weight profile** — mooted; the structure was identified without it.
4. **Stacking** — decisive and surprising: v4/v5/v45 inject the same monomial
   msb∧t0 and all sit at the same ¾ cap; a second independent monomial (v3x)
   removes the cap (density 0.98). The residue tracks independent nonlinear
   monomials, not wrap count. Bonus, derived from C29 then verified: the
   broken variant keeps C15 constancy exactly (C35) — the ~51% penalty does
   not forfeit the free exponent register.
5. **0.716 the constant** — closed: instance-dependent and width-dependent;
   it is a finite-size snapshot of the approach to ¾ (98.9% of cap at
   n_exp=6, beside the baseline at 98.7% of ½).

Left open (logged in §RS): the live quadrants sit at 0.91–0.98, below the
generic ~0.995 — a smaller, deeper deficit, presumably further conditional
levels. The general 1 − 2^−(k+1) ladder is a conjecture, not claimed.

## 12. `[x]` Windowed / table-lookup arithmetic — DONE, and it corrected C29

**The criterion was incomplete, and this is the item that found it.** See
`NOTES.md` §WD; files `windowed_arith.py`, `test_windowed.py`,
`experiments/experiment_windowed.py`. Claims C36–C39; C29 regraded to NARROWED.

The literal question has a yes answer that turns out not to decide anything:
**both** windowed designs satisfy V² = id exhaustively, and one of them loses
the invariance completely. So V²=id does not discriminate.

- **C36, the repaired criterion.** C29's *other* hypothesis — each block
  controlled on **its own qubit** — is load-bearing. The general condition is
  that the identity-tail block's dependence on the exponent register be
  **affine**: identity (one fresh control) or constant (uncontrolled) qualify,
  OR does not. Step 9 could not see this because every block it tested had
  exactly one control, where the two hypotheses coincide.
- **C37, lookup (Gidney-style) is covered, and more strongly than expected.**
  The tail table is all-ones, so the QROM permutation is `s ^= 1` regardless of
  j and the block never reads its window. The tail is **dead**, not merely
  constant: support confined to exponent bits {0…α−1}, verified 3/3 with α =
  1,1,2. Directly distinguishable from C24, which says the standard
  construction's support *contains* the all-ones-tail vector — same modulus and
  base, opposite geometry.
- **C38, the honest caveat.** Cost is **bounded and 2-periodic** in the
  tail-window count, not the exact constant C24 gives, because the tail applies
  W unconditionally and W²=id. 28078/62680/28078/62680 at K=0..3; K=0 and K=2
  support **sets identical**. Control β=3 grows.
- **C39, the cause is not windowing.** At w=1 — no window at all — the lookup
  form still has a dead tail. The difference is **emitting the multiply-by-1
  branch instead of optimising it away**. Skipping it (`SelectModExp`) grows
  **exactly ×4.00 = 2^w per window**, the derived rate (OR of w bits has Walsh
  sparsity 2^w), with density pinned at C30's ½ cap — i.e. a β=1 instance is
  put onto the generic β>1 curve. Mirror image of C32: there the standard
  reduction was accidentally cheap, here standard practice would be
  accidentally expensive.

**A prediction was refuted and the refutation was informative.** P5 (restoring
the j=0 branch restores the dead tail) was a conjunction: the flatness held
(support ~32k at every K, no growth) but the dead-bit half failed — the live
window's other bit stays live. Diagnosed, not waved away: the **activation
ancilla** is itself a scratch qubit the pullback ranges over; the lookup design
cancels it (2^w CNOT contributions into s XOR an even number of times), the
select design does not (the ancilla gates a block). Confirmed by checking the
block on the act=0 half-space, where the e1 dependence vanishes.

**VALIDATED AT SOURCE** (arXiv:1905.07682 §3.5, arXiv:1905.09749). C39's
mechanism is Gidney's own stated rationale — *"this also removes the need for
the multiplications to be controlled, because the table lookup can evaluate to
the factor 1 in cases where none of the exponent qubits are set"* — so C37–C39
describe the construction people actually propose. The joint `table[ei, mi]`
indexing we simplified away is confirmed harmless (in a tail window every
ke = 1, so the outer index goes degenerate). His *relabelling* swap, which we
modelled as a physical swap, turns out to be what makes the block an involution
at all, and his `if a is not target: swap(a, b)` line is C38's 2-periodicity
appearing in the real compiler. **Remaining real gap:** the analysis assumes a
**unitary** unlookup; Gidney's measurement-based uncomputation is not unitary,
so C8 does not apply to it as written.

## 12b. `[x]` Is the recurring GF(2) linear/nonlinear pattern one theorem?

**Mostly yes, and it made two of our claims cheaper to state and one of them
wrong.** See `NOTES.md` §GF; file `experiments/experiment_gf2law.py`. Claims
C40, C41; C30/C33 downgraded from "our mechanism" to "cite Carlet Prop. 29".

The unifying *fact* is the defining property of the transform — Walsh
characters of GF(2)ⁿ **are** the affine functions — with two corollaries:
**(a)** affine symmetry confines the support (**classical**: Carlet,
Proposition 29; our C30/C33 are instances, and "partially bent" is the name for
the extremal case); **(b)** affine gating costs one Walsh coefficient, which is
why §WD measured exactly ×2^w per OR-gated window.

- **C40** extends (a) to *conditional* structures: per-cell linear structures
  w_u ⟹ support avoids E = {z : w_u·z = 1 ∀u} ⟹ **density ≤ 1 − 2^−d**,
  d = dim span{w_u}. Verified 4/4 planted with 0 violations, and shown
  affine-invariant (tested, not asserted).
- **C41**: consistency is a **parity** condition — an odd dependency among the
  w_u destroys the cap outright (0.9845 vs 0.8671 at the same cell count).
- **It corrects §RS.** The conjectured 1 − 2^−(k+1) ladder is the wrong shape:
  the parameter is the span dimension, not the conditioning depth.
- **It explains C34 without modexp.** v4's two slice structures (msb⊕anc and
  msb) span d = 2, so E is exactly the quadrant {z_msb=1, z_anc=0} that C34
  found empty, with cap ¾.

**Prior art — pursued as far as open access goes; risk now low.** The one
unread item is Carlet–Tarannikov, DCC 25:263–279 (2002), paywalled at Springer
with no self-archived preprint. Checked instead: **Carlet's own book** (the
comprehensive survey by the same author, which cites that paper on pp. 205,
206, 319 and reproduces its Def. 47 / Prop. 60 — §5.5 read in full), and **the
paper's own abstract**. A covering sequence is a *single global* λ; a *partial*
covering sequence allows two levels; ours has a **different structure vector
per cell**, which is neither. The abstract's stated contributions are
resiliency/correlation-immunity characterisations and constructions — not a
support-confinement law. Also worth citing and distinguishing:
**Maiorana–McFarland** (restrictions *affine* on each coset — same flavour,
stronger hypothesis, different conclusion).

Still state C40/C41 at corollary altitude (Prop. 29 + coset decomposition) and
note they may be folklore. If someone gets institutional access, reading the
DCC body is the last loose end — a ten-minute job.

## 12c. `[ ]` Can the dead 2^−(α+1) fraction be pruned EARLY?

Opened by §BI/C42. For the real Shor input state a fraction 1 − 2^−(α+1) of the
Walsh support contributes **exactly zero** to ⟨O⟩ — 75% at α=1, 87.5% at α=2,
verified 3/3. But PPS propagates backwards and only meets the input state at
the end, so those terms are carried at full cost and *then* discarded: peak
memory (C18) is unchanged. **The question is whether the dead set can be
recognised early.** If a cheap invariant identifies "this branch will end with
exponent support" partway through propagation, that is a real constant-factor
win of 4× (α=1) to 8× (α=2) in peak memory, on top of everything else — and
unlike δ-truncation it is **exact**, not approximate.

> **The naive version is DEAD — answered the same day it was written, by
> derivation, before any effort went into it.** Exponent support is **not
> monotone** under back-propagation, so a term cannot be pruned when it
> acquires it. Through `CCX(a,b,c)` the pullback of `Z_c` is
> `½(Z_c + Z_aZ_c + Z_bZ_c − Z_aZ_bZ_c)`, and the pullback of `Z_aZ_c` is
> `½(Z_aZ_c + Z_c + Z_aZ_bZ_c − Z_bZ_c)` — because `Z_a·Z_a = I`. Verified
> exactly: `Z_aZ_c` produces the term `Z_c` with coefficient +0.50, carrying
> **no** support on the control qubit a. A term with exponent support can
> therefore lose it and go on to contribute, so zeroing it early is simply
> wrong, not merely suboptimal.

What survives the above: is there a *certificate* weaker than "has exponent
support now" that predicts the final exponent support of a whole branch? The
C23/C24 mechanism is known analytically (the tail is entered through a single
parity bit), so the dead set is characterised at the END — the question is
whether that characterisation can be pushed backwards through the propagation.
Unclear, and much less likely to be cheap than it looked. A must-fail control
is easy: the β>1 case should show no such structure.

Related: the same framing applies to any input state via ⟨O⟩ = Σ_z c_z ∏ δ_i,
so a partially-biased register interpolates between "free" and "carried".

## 12d. `[x]` Why is the peak ratio exactly `rot = 2·perm − 2` for modexp?

**SOLVED 2026-08-08 — and it resolved in an afternoon, as the note below
guessed.** See `NOTES.md` §PK; file `experiments/experiment_c17_deficit.py`;
claim **C44**; C17 regraded from "the factor is empirical" to derived; Paper A
§5 rewritten as Proposition 2.

**The identity.** Only Z_c-carrying strings meet the gadget's H and leave the
diagonal, and the four T gates on c are rotations *about* Z_c, which act inside
span{X_c, Y_c} — a closed 2-dimensional space. So they branch **once between
them**, not 2⁴ times, and

```
    N_max^rot = 2·N_max^perm − |B|,    B = { z ∈ S : z_c = 0 }
```

with c the target of the gadget containing the rotation-level peak. **The
deficit is a countable set, not a constant.** Verified 9/9, and sharper than a
count: at the peak every Pauli has X-support ∅ or exactly {c}, and folding the
X_c/Y_c partners recovers the perm peak set *set for set*.

**Why 2.** B = {Z_x0, Z_x0·Z_e0} in 6/6 modexp instances — and those are
exactly the **|coefficient| = ½ Walsh terms**, i.e. the dominant Fourier modes
that §W3/W4 had already singled out from the *final* spectrum by an unrelated
route. They live off the scratch register, never acquire Z_c, never double.

**The correction it forces.** The constant belongs to the **observable**, not
to modexp: same circuit, observable moved to b0 → deficit 4004, anc → 4014,
t0 → 0 (ratio exactly 2, like the adders).

**Left open, deliberately ungraded:** `|B| = 2` itself is measured, not proved.
Proving it needs a characterisation of which peak-time strings avoid the
scratch register — a statement about the peak moment, not about the final
spectrum. The identity above is the theorem; the membership is evidence.

**Method note worth keeping.** Pass 1 had the identity right and the *pairing*
wrong: it compared the peak against the Z-set at the gadget **boundary**,
reconstructed with ~100 lines of gate-index bookkeeping. 0/9. The peak actually
sits mid-gadget, and the gadget target never needed reconstructing — it is the
single bit set in the peak's common x-mask. The checks that only counted B
passed 6/6 through both passes, which localised the failure immediately.

<details><summary>Original note, kept for the record</summary>

Opened 2026-08-08 by the C17 correction. The permutation-native peak reduction
is **2.000000 exactly** for ripple-carry adders (128/64, 512/256, 2048/1024) but
**1.9997** for modular exponentiation — and every modexp instance satisfies

```
  rot = 2·perm − 2       13666 = 2·6834 − 2
                         16386 = 2·8194 − 2
                        128138 = 2·64070 − 2
```

exactly, 3 of 3, across N = 5, 7, 15 and wildly different peak magnitudes. A
constant additive deficit of exactly 2 that is independent of instance size is
not a coincidence and is currently **unexplained**.

**The half that is understood.** At the gadget's Hadamard on target qubit c, a
Z-type string containing Z_c is carried to a mirrored X-sector partner while a
string without Z_c is untouched. So the doubling applies only to the
Z_c-containing subset, which bounds the ratio above by 2 and explains why the
adders attain it (there, at the peak, every string apparently contains Z_c) and
modexp does not.

**The half that is not.** Why the shortfall is exactly 2 rather than
instance-dependent. Two obvious candidates, neither checked: the identity string
Z^0 (which cannot acquire Z_c and so never doubles) plus one partner; or a
parity constraint pinning a second string. Cheapest first move is to dump the
peak-time Pauli set for the smallest modexp instance and simply *look at* which
two strings fail to double — this is a 2-line diagnostic, not a research
programme, and it either resolves in an afternoon or reveals something.

If it resolves, `PAPER_A.md` §5 upgrades from "upper bound of 2, empirical" to a
theorem, which is worth having: it is currently the only quantitative claim in
Paper A that is measured rather than derived.

</details>

*The "cheapest first move" above — dump the peak-time Pauli set and simply look
at it — was exactly right, and is what produced the answer. Both named
candidates (the identity string Z^0; a parity constraint) were wrong.*

## 12e. `[x]` Re-run the reachable sweeps now that the GPU makes them cheap

**DONE 2026-08-08. All three sub-items closed; every prediction confirmed.**
See `NOTES.md` §OS; files `experiments/experiment_c21_onset.py`,
`experiment_c7_scale.py`, `experiment_windowed_scale.py`. C21, C24, C7, C30
and C38 all gained evidence; C43 is new.

- **C21's onset at α = 3 and 4 — MEASURED (§OS1).** Locks at exactly
  n_exp = α+1 in all three instances: N=17 a=2 (r=8) and N=41 a=3 (r=8) at
  n_exp=4, N=17 a=3 (r=16) at n_exp=5, with strict growth at all 7 steps
  below. Two moduli for α=3, so it is not an N=17 artifact. Matched control
  (N=41 a=6, r=40=5·2³ — same modulus, same α, same width, only β differs)
  never locks. **Worth carrying into Paper B: at n_exp=4 the β=1 and β=5 rows
  differ by 66 in 3.4e7, so no single-width cost measurement can see the
  invariant — only the growth separates them.**
- **C24 strengthened to set level (C43, §OS2).** The proof's step (iv) gives
  the surviving coefficient explicitly with no |I| in it, so the supports must
  coincide as *sets* under the (z_rest, tailflag) encoding, not merely in
  size. Derived before measuring, confirmed 3/3 bit-for-bit over supports of
  4.2M / 8.4M / 33.5M.
- **C7 extended to 30 qubits (§OS4).** Two series, slope 1.0054 / 1.0083,
  pooled 1.0062 bits/qubit over q = 15..30; q=30 gives |support| =
  536,271,623 at density 0.49944 in 851 s. All six original rows reproduce
  bit-for-bit on the exact integer path. `PAPER_A.md` §11.1's concession that
  the circuit series is modest can now be softened.
- **C38's 2-periodicity extended to K = 5 (§OS6)** — four full periods.
- **Method (§OS3): the must-fail control caught a vacuous test of mine.** The
  tail-confinement check has only one possible answer at |I| = 1; the control
  passed when it had to fail. Re-scoped, numbers unchanged, conclusions
  narrowed. Third instance of this failure mode in the project.
- **Infrastructure (§OS5):** int32 permutation replay (the *replay*, not the
  FWHT, is what binds memory — int64 needs ~24 GiB at n=30 and fails),
  `accel.pullback_stats`, and exact/cached `lab.measure.support(exact=True)` /
  `lab.measure.stats`. All gated in `test_accel.py`.

**Left open by this item, logged not claimed (§OS4).** The missing fraction of
the C30 hyperplane, 1 − 2·density, falls ~2.2× per bit of modulus in the
n_exp = 2 series and does *not* in the n_exp = 3 series; matched pairs show
N = 33, 35, 77 dropping 8–11× between n_exp = 2 and 3 while N = 21 barely
moves. Cheap follow-up: sweep n_exp at fixed N ≥ 33 with β > 1 and see whether
the deficit really has a cliff there.

## 12f. `[ ]` Two-GPU split FWHT — only if a single run ever needs n = 31+

**Deliberately not done.** Memory doubles per qubit, so the second A4500 buys
exactly **one** more qubit (n = 30 → 31). Measured: float64 FWHT needs 12 GiB
at n = 30 (fits one card) and 24 GiB at n = 31 (does not).

The split itself is easy if it is ever wanted. Halve the array by its **top
bit**; then for every level h < size/2 the butterfly pairs stay inside one half,
so both cards run independently with no communication, and only the **final**
level pairs element i of half 0 with element i of half 1. That is one
peer-to-peer exchange of half the array (~8 GiB each way at n = 31), a couple of
seconds over PCIe.

**Do the cheap thing first.** `accel.wht_exact` already halves memory by using
int32 — legitimate because the FWHT of ±1 data is exactly integer-valued — which
buys the same +1 qubit with no split, and additionally makes the support test
exact rather than thresholded. Only build the split if something genuinely needs
n ≥ 32, and note int32 overflows past n = 30 (bound 2ⁿ vs 2³¹), so a split at
n = 32 needs int64 and therefore both cards anyway.

Meanwhile the second card is not idle-by-design: **run independent sweeps
concurrently**, `CUDA_VISIBLE_DEVICES=0` and `=1`. That is what TODO 12e wants
and it needs no new code.

## 12g. `[ ]` Does the hyperplane deficit fall off a cliff between n_exp 2 and 3?

Opened by §OS4, and cheap. C30 confines the support to a hyperplane, so the
interesting quantity is how much of that hyperplane is *missing*:
1 − 2·density. Across the n_exp = 2 series it shrinks with striking regularity,
~2.1–2.3× per extra bit of modulus (0.0544 → 0.00112 over n = 3..8). Across the
n_exp = 3 series it does not. The matched pairs are the odd part: same N, same
a, only n_exp differing, and N = 21 barely moves (0.01074 → 0.01004) while
N = 33, 35 and 77 each drop **8–11×**.

Either there is a threshold in n relative to n_exp, or it is instance noise
across four points. **Design:** fix N ≥ 33 with β > 1, sweep n_exp by 1 (this
varies exactly one parameter, unlike the two series, which vary N), and look at
the deficit. q = 3n + 4 + n_exp, so N = 33 reaches n_exp = 8 within budget.
Must-fail control: a β = 1 instance, where C21 locks the support and the
deficit must therefore blow up rather than shrink.

Low stakes — it refines the constant in C7, it does not touch Θ(2ⁿ) — but it is
one sweep and the regularity in series A is too clean to leave unexamined.

## 13. `[ ]` Third simulation method on the r = β·2^α invariant

Paper B open problem 4. Two structurally unrelated methods (PPS, MPS) keying
on the same arithmetic invariant is suggestive; a third would make "property
of the algorithm, not the simulator" hard to argue with. Decision diagrams are
the natural candidate (MQT DDSIM): rerun the controlled design (fix N, vary a
so only r changes) and see whether memory keys on β. Mostly integration work,
no new theory.

## 14. `[ ]` Through the inverse QFT — the frontier of the exact model

> **Cîrstoiu checked and it does NOT help here** (2026-08-08, body read). It was
> pulled specifically on the guess that a mixed GF(2)^t × Z/2^t character basis
> would bring it into range. It does not: their group indexes **circuit
> parameters** and their results are about **ensembles** of circuits, and Shor's
> circuit is fixed with no free parameters — there is no ensemble to average
> over. Useful for framing (Pauli-path methods as harmonic analysis on a group,
> which strengthens Paper A's "the analogy becomes an identity" positioning),
> not for the technical problem. Do not re-pull it for this.

Highest risk, highest reach; a Paper C candidate. Everything so far stops
where the identity stops (F10: X/Y pullbacks leave the diagonal and blow up),
but real Shor measures the exponent register *after* the inverse QFT. The
structure is suggestive: the Walsh basis is the character group of (Z/2)^t,
the QFT diagonalises translation on Z/2^t, and the whole cost story is already
2-adic. Is there an exact characterisation of the pulled-back observable in a
mixed character basis (GF(2)^t Walsh × Z/2^t Fourier)? May simply be dense and
structureless — but it is the natural next question the method itself asks.

## Housekeeping

- `[x]` Yao.jl backend — closed as far as possible without installing it. Its
  docs list Toffoli under **"Clifford Gates: Two-qubit gates"**, wrong on both
  counts (Toffoli is neither two-qubit nor Clifford), so the doc is unreliable;
  and PauliPropagation.jl, which it most plausibly wraps, has **zero** mentions
  of Toffoli/CCX/CCZ and no permutation gate type. No library documents or
  exploits diagonal closure. A definitive check would need Yao installed.
- `[x]` ~~Decide whether `.env` should be tracked~~ — **MOOT, deleted
  2026-08-08** along with Julia. It existed only to `LD_PRELOAD` Julia's
  libstdc++.
- `[ ]` `experiment.py` §3, `experiment2.py` §§2–3 and all of `experiment3.py`
  are retracted. They carry warning headers; consider deleting once the papers
  are drafted and nothing references them.
- `[x]` Organize experiments into a proper file structure, with a general
  engine — **DONE 2026-08-08.** `lab/` package (harness with enforced
  predictions/controls, cached measurement, GF(2) structure tools, modarith,
  wrap-registry variants, null models); experiments moved to `experiments/`
  unchanged (records, not code to DRY); `experiment3.py` to `archive/`; two
  new gate suites (`test_lab.py` pins the engine to logged numbers,
  `test_claims.py` re-verifies the headline claims); `experiments/TEMPLATE.py`
  for new work. Deliberately NOT a `src/` package — the unit of value here is
  the readable record, and root-level instruments keep every historical
  import working.
- `[x]` Split claims into a separate file to keep a single source of truth —
  **DONE.** `CLAIMS.md` is now canonical: Paper A claims, Paper B claims, and a
  retracted/dead section, one full row per claim with a pointer to where it is
  established. Both abstracts replaced their tables with a one-line pointer, so
  the duplicate-row drift that produced a stale C7 cannot recur. `NOTES.md` is
  untouched and remains the historical narrative — statuses in its prose may be
  superseded; `CLAIMS.md` wins.
  - Still open: are there other files that should be split out as well? Perhaps
    create some sort of wiki?
- `[x]` ~~Check if julia and dependencies are still needed~~ — **DONE
  2026-08-08. Removed.** No `.py` file imported `juliacall`; the prior-art work
  on PauliPropagation.jl was done by *reading* its source, so nothing
  reproducible depended on it. Dropped `juliacall`, deleted the orphaned 1.1 GB
  depot and `.env`. Venv 1.6 GB → 534 MB, all suites pass, and trap 2 in
  `HANDOFF.md` (the `LD_PRELOAD`/`longdouble` segfault) is retired.

- `[ ]` **Six more dependencies are also unused — decide.** An import scan
  found *only* numpy is imported anywhere: `click`, `matplotlib`, `numba`,
  `qiskit`, `quimb`, `scipy`, `stim` all have **zero** imports. Not removed
  unilaterally because `qiskit` and `stim` were the *source-reading* tools for
  the prior-art sweeps (TODO 1) and may be wanted again; the rest look like
  genuine leftovers. Recommendation: drop `click`, `matplotlib`, `quimb`,
  `scipy`; keep `numba` only if the CPU fallback is ever parallelised; keep or
  drop `qiskit`/`stim` on how likely another source-level sweep is.

- `[x]` **GPU acceleration — DONE 2026-08-08.** `accel.py` + `test_accel.py`
  (suite 9). Both hot paths are memory-bound array passes, so a card with ~10×
  the bandwidth wins by about that: measured **7.1× / 6.9× / 11.5×** on
  permutation replay at q = 17/19/21 and **11.8× / 12.8× / 14.7×** on FWHT at
  2^24/26/28 (RTX A4500). The speedup **grows with n**, which is the useful
  direction. Opt-in via `LAB_GPU=1`; `walsh.py` stays the reference and the GPU
  path is gated against it (exact equality on permutations, identical support
  sets, plus must-fail controls for non-permutation input and oversized
  registers). Capacity: n ≤ 30 on one 20 GiB card.
