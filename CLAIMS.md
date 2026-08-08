# Claims ledger — single source of truth

**This file is the canonical record of claim statuses for both papers.**
`ABSTRACT.md` (Paper A, the Walsh cost model) and `ABSTRACT_SHOR_2ADIC.md`
(Paper B, Shor / 2-adic) no longer carry ledger tables of their own; they point
here. `NOTES.md` is the **historical working narrative** — it records what was
believed at each point in time and therefore contains **superseded statuses**.
Where `NOTES.md` prose disagrees with a row below, **this file wins**.

Each row must be independently reproducible before the abstract goes anywhere.

**All figures below are post-bugfix** (the θ=π propagator bug, `NOTES.md`
STATUS 2). Anything citing pre-fix numbers is void.

Every claim has its full row in exactly **one** section. Claims that both papers
rest on are cross-referenced by ID rather than duplicated — duplicate rows are
what caused the ledgers to drift once already (a stale C7 row, removed in
`3606e9e`).

The "Where established" column cites only pointers stated in the source
material (`ABSTRACT.md`, `ABSTRACT_SHOR_2ADIC.md`, `TODO.md`, `HANDOFF.md`,
`README.md`, `NOTES.md` section headers). A blank cell means no pointer was
recorded, not that none exists.

**Paths note (2026-08-08 restructure).** Experiment scripts moved into
`experiments/` with filenames unchanged; a bare `experiment_*.py` reference
below means `experiments/experiment_*.py`. The headline rows are additionally
re-verified by `test_claims.py`, which runs as part of the correctness gate.

---

## Paper A — Walsh cost model (`ABSTRACT.md`)

| # | Claim | Status | Evidence | Where established |
|---|---|---|---|---|
| **C8** | **PPS term count = Walsh sparsity of the pulled-back bit function, exactly** | **PROVEN (analytic, numerically confirmed)** | 6/6 instances, supports identical, maxerr ≤ 6.7e-16 | `NOTES.md` "F9 — THE RESULT"; `experiment_c8.py` |
| C1 | Toffoli-compiled arithmetic has exact Z-closure | **established** | 0 non-Z terms | `NOTES.md` STATUS 2, "Corrected results" |
| C2 | Fourier-compiled modexp **also** has exact Z-closure | **established** (reverses old C2) | 0 non-Z terms | `NOTES.md` STATUS 2, "Corrected results" |
| C6 | Z-closure follows from the *unitary* being a basis permutation, not from the gate set | **established** | holds for both compilations | `NOTES.md` STATUS 2, "Corrected results" |
| C10 | Adder collapse is affineness, not permutation-ness | **established** | sparsity 1 ⟺ affine; `test_walsh.py` [E] | `NOTES.md` "F9 — THE RESULT" §3; `test_walsh.py` [E] |
| C11 | Walsh gives an exact cost model ~10³× cheaper than the PPS run | **RESTATED** | exact for *final* support; peak is a distinct larger quantity — see C17 | `NOTES.md` "F9 — THE RESULT"; `NOTES.md` "R1 UPHELD" |
| **C17** | Peak ≠ final; permutation-native PPS reduces peak by a factor **≤ 2** and makes it a Walsh quantity | **established; the factor is EMPIRICAL, not proved — "halves exactly" was a rounding artifact, corrected 2026-08-08** | exactly **2.000000** for adders (128→64, 512→256, 2048→1024); **1.9997** for modexp (13666→6834, 16386→8194, 128138→64070), satisfying `rot = 2·perm − 2` in 3/3 — the unexplained part. Upper bound of 2 is argued from the gadget's H exchanging Z- and X-sectors; not yet a proof | `perm_pps.py`, `test_perm_pps.py`; `PAPER_A.md` §5 | `NOTES.md` "R1 FIX — permutation-native PPS"; `perm_pps.py`, `test_perm_pps.py` |
| C12 | PPS-hardness ≡ linear-cryptanalysis resistance | **VALIDATED, with caveat** | endpoints exact; Parseval relation is a bound, tight only for flat spectra. NB the 0.74 figure is stable across N at fixed n_exp, **not** across widths (0.50 at n_exp=1) | `NOTES.md` "C12 VALIDATED" and §X4 (the correction); `experiment_c12.py` |
| **C25** | **S ≥ (1 − NL/2ⁿ⁻¹)⁻²**: published nonlinearity ⟹ PPS cost lower bound, compilation-independent | **PROVED + verified** | tight at affine (1) and bent (2ⁿ); AES NL=112 reproduced; bent circuits give full support end-to-end (64/256/1024) | `NOTES.md` §X (X1, X2); `experiment_crypto.py` (TODO step 5) |
| C26 | The bound is weak away from the extremes | **established, must be stated** | AES 64 vs 239 actual; modexp 4 vs 3086. Uses only max\|c\|, discarding the rest of the spectrum | `NOTES.md` §X3; `experiment_crypto.py` (TODO step 5) |
| C16 | Heavy-tailed spectrum explains the δ non-monotonicity (C14) | **established** | 4 of 15493 coefficients give ⟨O⟩ exactly; the other 15489 sum to zero | `NOTES.md` "F13 — the spectrum shape explains F11" |
| C13 | Identity covers **all** Z-type observables (multi-qubit too), fails for X/Y | **established** | 4/4 Z-type exact; X/Y pullback fully non-diagonal | `NOTES.md` "F10 — scope of the Walsh identity"; `experiment_scope.py` Q1 |
| C14 | Truncation error is **non-monotonic** in δ: δ=1e-1 exact w/ 34 terms, δ=1e-3 off by 0.285 w/ 23482 | **established** | mechanism: small coefficients cancel as a set | `NOTES.md` "F11 — truncation error is NON-MONOTONIC"; `experiment_scope.py` Q2. Sharpened by `NOTES.md` §W1 (incremental truncation also changes what branches) |
| C5 | Truncation can violate \|⟨O⟩\| ≤ 1 | **survives, quantified** | 4/18 runs inadmissible, all at mild δ | `NOTES.md` "F11" and "F5"; `experiment_scope.py` |

**Also relied on by Paper A, full rows in the Paper B section:** C7 (Θ(2ⁿ) for
generic r — the original "results are pre-asymptotic" framing is retracted, see
the Retracted section), C15 and C18 (both moved wholesale to Paper B), C30, C31
and C32 (the density-½ linear structure and the one compilation-dependent cost
effect).

Paper A's ledger stated C30 and C32 in their general, cost-model form; those
statements are preserved verbatim in the trailing notes of the corresponding
Paper B rows.

---

## Paper B — Shor / 2-adic (`ABSTRACT_SHOR_2ADIC.md`)

| # | Claim | Status | Evidence | Where established |
|---|---|---|---|---|
| **C15** | Cost is set by the 2-adic structure of r, at circuit level | **CONFIRMED, controlled** | fix N vary a: r=2 ⟹ support *exactly* 15549 over 256× dim growth; r=6 ⟹ 4.00×/step; reproduced at N=21 and N=5 | `NOTES.md` "C15 CONFIRMED AT CIRCUIT LEVEL", "C15 REFINED"; `experiment_c15.py`, `experiment_c15b.py` |
| **C21** | Invariance has a precise onset: support locks at **n_exp = v₂(r)+1** | **established** | α=1,2 exact (α=2 grows 15493→32143 then freezes); α=3,4 still growing at max reachable width, as predicted | `NOTES.md` "C15 REFINED"; `experiment_c15c.py`, `experiment_c15d.py` (TODO step 4) |
| C22 | λ(N) a power of two ⟹ **every** base free; ⟺ N = 2^a × distinct Fermat primes | **established** | N=15: 7/7 bases free; N=21: 3/11. Odd semiprimes = p·q both Fermat, smallest 15 | `NOTES.md` "C20 STRENGTHENED"; `experiment_c15d.py` (TODO step 4) |
| **C18** | Holds for **peak memory**, not just final support | **established** | N_max exactly 24369 at n_exp = 2,4,6,8 (64× dim growth); control r=6 grows 4.02×/step | `NOTES.md` "R1 consequence — C15 SURVIVES for peak memory"; `perm_pps.py` (see C17) |
| C7 | Generic r ⟹ Θ(2ⁿ); results are not pre-asymptotic | **established** | 24 qubits via Walsh; density 0.473→0.498→½, slope 1.008 bits/qubit | `NOTES.md` "C7 RESOLVED"; `experiment_c7.py` |
| **C30** | The ½ ceiling is a **linear structure** w = b_msb ⊕ anc, not an algorithmic constant | **established; the mechanism is CLASSICAL — cite, do not claim** (Carlet, *Boolean Functions for Cryptography and Coding Theory*, **Proposition 29**: D_e f null ⟺ supp(W_f) ⊆ {0,e}^⊥; such functions are Carlet's **partially bent**). Our contribution is finding it in modexp pullbacks and using it as a PPS cost cap — see §GF | GF(2) rank n−1 in every instance; g(y⊕w)=g(y) pointwise; observable-independent; present in both compilations; random f has full rank | `NOTES.md` §L, §L2; `experiment_linstruct.py` (`gf2_rank_and_kernel`, per `HANDOFF.md`). Paper A stated the general form: "Reversible-arithmetic pullbacks can carry **linear structures**, capping density at 2^−k and bounding them away from bent — **established** — modexp: GF(2) rank n−1, w = b_msb⊕anc, g(y⊕w)=g(y); adder: kernel dim 5 ⟹ density ≤ 2⁻⁵. Bent functions have none, which partly explains modexp's 0.74 vs random's 0.97" |
| **C31** | The structure is forced by (a) msb-flip commuting with mod-2^m addition, (b) XOR-only coupling to anc, (c) msb excluded from the cswaps | **established** | (a) verified exhaustively, 0/256 violations; adding further XOR couplings leaves it intact | `NOTES.md` §L2, "Why the structure is robust"; `experiment_reduction.py` (TODO step 10) |
| **C32** | Breaking it needs **nonlinearity in the msb**, and costs ~51%: density 0.473 → 0.716 | **established** | Toffoli variants break it (full rank); cswap variant does not (linear); all four still compute a^e mod N | `NOTES.md` §L2, "Breaking it requires NONLINEARITY in the msb"; `experiment_reduction.py`, `experiment_reduction2.py` (TODO step 10). Paper A stated it as: "**The only compilation-dependent cost effect found.** Destroying the linear structure costs ~51% (density 0.473 → 0.716) with the computed function unchanged — **established** — needs nonlinear coupling to the msb; linear couplings (CNOT, cswap) leave it intact. Constant-factor only — Θ(2ⁿ) unaffected" |
| **C27** | Invariance survives truncation: exact for terminal thresholding at every δ; peak cost unchanged under incremental | **established** | counts identical across n_exp at all 7 δ values; N_max SAME at every δ; control diverges | `NOTES.md` §T (T1, T2); `experiment_c15_trunc.py`, `experiment_c15_trunc2.py` (TODO step 7) |
| C28 | Accuracy degrades first for the *wider* circuit at aggressive δ | **established, must be stated** | ⟨O⟩ → 0 at δ=1e-1 for n_exp=4,5 while n_exp=3 stays exact; more gates ⟹ more incremental truncation | `NOTES.md` §T2; `experiment_c15_trunc.py`, `experiment_c15_trunc2.py` (TODO step 7) |
| F12 | Function-level dichotomy is absolute | **established** | r=4: sparsity 4 constant to t=24; odd factor: density 1.000000 | `NOTES.md` "F12 — the 2-adic dichotomy at function level"; refined by §I and `experiment_c15_intermediate{,2,3}.py`, `experiment_c15_realfn.py` |
| C19 | Same r = β·2^α invariant governs MPS simulation | **established, cited** | Dang et al. §4: α = trailing zeros, β = odd part "cannot be localised"; §5.2: memory ∝ β² | Dang, Hill & Hollenberg, [arXiv:1712.07311](https://arxiv.org/abs/1712.07311), §4 and §5.2 |
| C20 | N=15 is a degenerate benchmark **for every base**, forced by the modulus | **established, strengthened** | 100% of its bases free; smallest product of two Fermat primes | `NOTES.md` "C20 STRENGTHENED"; `experiment_c15d.py` |
| **C23** | Mechanism: identity blocks apply one involution V, so the circuit depends on the identity tail only through a **parity bit** | **PROVED** | V²=id since V=A⁻¹SA with S disjoint transpositions; blocks commute; character averaging confines support to z_I ∈ {0,1_I} | `NOTES.md` "C15 MECHANISM SOLVED" and "FORMALISED"; `experiment_c15_proof{,2,3,4}.py` (TODO step 3) |
| **C24** | Consequence: support confined to z_I ∈ {0, all-ones}, so size is independent of tail length | **PROVED + verified** | derived before testing; 0 violations in 7/7 β=1 instances, halves individually constant; controls fail with 48189 / 1556046 | `NOTES.md` "FORMALISED — this is now a theorem"; `experiment_c15_proof5.py` (TODO step 3) |
| **C29** | ~~The theorem needs only **V²=id**~~ → The theorem needs V²=id **plus** each block being controlled on its own qubit; it covers any construction with an involutive repeated block *so controlled* | **NARROWED 2026-08-08 (was "PROVED + verified")** | the step-9 evidence stands (synthetic involution constant across k=1..4 at two moduli; order-3 block grows ~2×/block; vacuity passes) but every block tested had exactly one control, so it could not separate the two hypotheses. `SelectModExp`'s tail block has V²=id exhaustively and still loses the invariance (×4.00/window) | `NOTES.md` §G (annotated) and §WD; `experiment_c15_general2.py` (TODO 9), `experiment_windowed.py` (TODO 12) |
| **C36** | The general criterion is **affine control**: V²=id is necessary but not sufficient; the identity-tail block's dependence on the exponent register must be affine — identity (one fresh control qubit) or constant (uncontrolled) qualify, **OR does not** | **established; supersedes C29's phrasing** | LOOKUP tail block is window-independent and keeps the invariance; SELECT tail block is gated on OR(window) (V₀ ≠ V₁ = V₂ = V₃, exhaustive) and loses it, with V²=id holding in both | `NOTES.md` §WD; `windowed_arith.py`, `test_windowed.py` [D]/[E], `experiment_windowed.py` (TODO 12) |
| **C37** | In table-lookup (Gidney-style) modexp the identity tail is **dead, not merely constant**: support confined to z_exp on bits {0…α−1}, containing **no** z with any tail bit set — strictly sharper than C24 | **derived before measuring, then verified 3/3** | live exponent bits exactly [0], [0], [0,1] at (N,a) = (5,4), (7,6), (5,2) with α = 1,1,2, at K = 1 and 2. Same modulus/base contrast: standard contains the all-ones-tail vector 1_I, lookup contains no tail bit at all **Model validated against arXiv:1905.07682 §3.5 / arXiv:1905.09749**: the joint `table[ei, mi]` indexing we simplified away is harmless (every ke=1 in a tail window ⟹ outer index degenerate). | `NOTES.md` §WD; `experiment_windowed.py` P1/P4 (TODO 12) |
| **C38** | For that construction the invariance is **bounded and 2-periodic in the tail-window count**, not exactly constant, because the tail applies W unconditionally and W²=id | **derived before measuring, then verified** | 28078/62680/28078/62680 (N=5 a=4) and 13563/62570/13563/62570 (N=7 a=6) at K=0..3; K=0 and K=2 support **sets identical**, not just equal in size; density → 0. Control β=3 grows 84827 → 1009765 → 4156585 **Corroborated at source**: Gidney's relabelling swap is what makes the block an involution at all, and his `if a is not target: swap(a, b)` emits a physical swap exactly when the block count is odd — C38's parity in the real compiler. | `NOTES.md` §WD; `experiment_windowed.py` P2 + C1 (TODO 12) |
| **C42** | For the **real Shor input state** (exponent register in \|+⟩, so δ=0 there) every z with exponent support contributes exactly 0, giving a useful fraction of exactly **2^−(α+1)** — cost and useful work are both constant in n_exp | **derived from C24, then verified 3/3 across widths** | 75.03% / 74.99% dead at α=1 and 87.51% at α=2; useful counts 3883 / 3879 / 4014 constant in n_exp; the α=2 row jumps exactly at n_exp = α+1, reproducing C21's onset independently. **Not a speedup**: PPS meets the input state only at the end, so peak memory (C18) is unchanged | `NOTES.md` §BI |
| — | *Dictionary:* a **biased product input state** makes ⟨O⟩ = Σ_z c_z ∏_{i∈z} δ_i, i.e. exactly the p-biased WHT of Gangopadhyay et al. (JAMC 2023) with δ_i = ⟨Z⟩_i | **verified numerically** | random product state, N=7 a=6 n_exp=3: direct statevector vs weighted Walsh sum agree to 8.2e-14. Their Theorem 10 (S-spectrum not EA-invariant) explains why §W's weight truncation had to fail while δ-truncation is exact | `NOTES.md` §BI; paper in repo root |
| **C40** | **Conditional-structure law.** If f restricted to each cell of a coset partition has linear structure w_u, then supp(f̂) avoids E = {z : w_u·z = 1 ∀u}, so **density ≤ 1 − 2^−d with d = dim span{w_u}** | **derived before measuring, verified 4/4 planted + affine-invariant; PRIOR ART OPEN** | caps 0.5/0.5/0.75/0.875 at d = 1,1,2,3 with **0 violations**; survives a random GF(2) change of basis (structures verified pointwise, |support| unchanged); random control shows no cap. Subsumes C30 (m=0,d=1) and C34 (d=2) | `NOTES.md` §GF; `experiment_gf2law.py` |
| **C41** | Consistency of that system is a **parity condition**: any dependency among the w_u with **odd** support makes E empty and **destroys the cap entirely** | **derived before measuring, verified as a must-fail control** | same 4 cells: odd dependency → E empty, density 0.9845; even dependency → d=3, cap 0.875, density 0.8671. The cap is set by span dimension + dependency parity, **not** by conditioning depth | `NOTES.md` §GF; `experiment_gf2law.py` C1 |
| — | *Superseded (§RS):* the residual ladder is 1 − 2^−(k+1) in the conditioning depth k | **CORRECTED by C40** | wrong parametrisation: two cells *sharing* a structure give d=1 and a ½ cap (measured 0.5000), not ¾. The ladder is the special case where each level adds one independent vector | `NOTES.md` §GF, §RS |
| **C39** | The cause is **emitting the multiply-by-1 branch**, not windowing: at w = 1 (no window at all) the lookup form still has a dead tail, and skipping j = 0 (`SelectModExp`) removes the β=1 advantage entirely — density pinned at C30's ½ cap, i.e. generic β>1 behaviour | **derived before measuring, then verified** | w=1, K=3: live exponent bits [0]. Select design grows **exactly ×4.00 = 2^w per window** (the derived rate: OR of w bits has Walsh sparsity 2^w), density 0.490–0.492 at every K, at both β=1 moduli **This is Gidney's own stated rationale** (1905.07682 §3.5): *"removes the need for the multiplications to be controlled, because the table lookup can evaluate to the factor 1 in cases where none of the exponent qubits are set"* — so C37–C39 describe the construction people actually propose. | `NOTES.md` §WD; `experiment_windowed.py` P1b/P3 (TODO 12) |
| **C33** | The exact-½ function-level densities (§I's "non-generic" rows) are exact **all-ones structures**: complementing every exponent bit acts as c ↦ (ρ−1−c) mod r with ρ = 2^t mod r, exact when the bit table is antipodal (prime N, even r ⟹ a^(r/2)≡−1) and/or inversion-symmetric (bit_j(x)=bit_j(x⁻¹)) | **PROVED + verified; the linear/affine→coset step is Carlet Prop. 29, cite it** | r=6: all-ones AFFINE at even t; r=3: all-ones LINEAR (inversion symmetry alone — no antipode needed); r=10 negative control clean (mapping right, symmetry broken, density 1.000); random and planted controls behave | `NOTES.md` §AF; `experiment_affstruct.py` (TODO 11.1) |
| **C34** | A single nonlinear monomial msb∧t0 does not destroy the linear structure — it **demotes it to a conditional one**: support exactly avoids the quadrant **{z_msb=1, z_anc=0}**, capping density at **¾**; the 0.716→0.742 trajectory is the approach to that cap. Two independent monomials remove the cap (density 0.98) | **verified 3/3; mechanism proved except one link** | quadrant exactly empty in 3/3 instances; conditional g(y⊕w)=g(y) on the t0=0 half PROVED (t restored ⟹ wraps inert) + 0 violations, must-fail control fails (7584); t0=1 slice structures exact and exhaustive (v4: msb; v5: msb⊕t1⊕anc) but their existence is argued from wrap-linearity, not yet proved in general | `NOTES.md` §RS; `experiment_resid1.py`, `experiment_resid2.py` (TODO 11) |
| **C35** | Breaking the linear structure and C15 constancy are **independent**: the broken variant keeps support exactly constant in n_exp | **PROVED (via C29) + verified** | wrapped a=1 block is A′⁻¹SA′, still a conjugate of the involutive cswap layer; v4 at N=7 a=6: support exactly 23488 at n_exp = 2, 3, 4 | `NOTES.md` §RS; `experiment_resid1.py` (TODO 11) |

**Depends on Paper A, full rows in the Paper A section:** C8 (the Walsh identity,
used as the measurement instrument) and C17 (`perm_pps.py`, which supplies the
peak-memory numbers and makes the n_exp = 8 measurement tractable).

---

## Retracted / dead

Nothing in this section may be cited. Kept so a fresh session does not
re-derive it. Statuses for the C-rows are as they stood in Paper A's ledger; the
F- and H-rows are as recorded in the `NOTES.md` honesty log and STATUS blocks.

| # | Claim | Status | Evidence | Where killed |
|---|---|---|---|---|
| C3 | Matched-instance A/B across compilations | **INVALID as run** | qubit counts differ (10 vs 15) — ancilla confound | `NOTES.md` STATUS 2, "CONFOUND that kills the A/B as designed" |
| C4 | Compilation gap is ~4 orders of magnitude | **REFUTED** | direction reverses; mostly an ancilla artifact | `NOTES.md` STATUS 2 |
| C9 | Truncation reverses the compilation advantage | **REFUTED** | rested on wrong ⟨O⟩ | `NOTES.md` STATUS 2 |
| F1 | Gate-level permutation circuits are FREE for PPS with Z-type observables | **NARROWED; explanation replaced** | right conclusion (adder collapses to 1 term), wrong mechanism (permutation-ness); the real reason is affineness — superseded by C8/C10 | `NOTES.md` STATUS 1 and honesty log |
| F2 | The QFT is the bottleneck, not the arithmetic | **RETRACTED** | bit-reversal bug (trap 4) plus a toy where the observable never met the arithmetic; the old 138 → 872306 scaling table is an artifact — do not cite it | `NOTES.md` STATUS 1, honesty log |
| F3 | δ is not a working dial for this family | **SUPERSEDED by F5/C5** | the slope argument rested on the toy circuits; the norm-violation evidence replaces it | `NOTES.md` STATUS 1 |
| F4 | PPS cost depends on the gate-level arithmetic (Fourier lacks Z-closure) | **RETRACTED** | killed by the θ=π PPS bug; **both** compilations end with zero non-Z terms | `NOTES.md` STATUS 2, honesty log |
| F6 | Toffoli compilation is ~8.5× better than Fourier | **RETRACTED** | direction reverses — Fourier is cheaper (451 vs 15493), and that gap is mostly an ancilla-count effect | `NOTES.md` STATUS 2, honesty log |
| F7 | Walsh sparsity is the driver, compilation sets the ceiling | **SUPERSEDED by F9/C8** | half right: sparsity is exactly the driver, but there is no separate compilation ceiling — the apparent one was an ancilla-count confound | `NOTES.md` STATUS 2, honesty log |
| F8 | Truncation destroys the structural advantage | **RETRACTED** | ⟨O⟩ was wrong (+1 where truth is −1), so the δ-sweep conclusions were meaningless | `NOTES.md` STATUS 2, honesty log |
| H1 | π/4 angles break truncation ("arithmetic is hard for PPS because all T angles are π/4") | **DEAD** | premise holds (adder: 42 branching gates, all \|sin\|=0.7071) but the conclusion is wrong — killed by F1, it all cancels | `NOTES.md` "Falsified hypotheses", honesty log |
| H2 | The π/4-uniform angle is what breaks the power law | **DEAD** | killed by Gharibyan et al. App. C / Fig. 9b — fixed correlated θ_X = π/4 and the power law still held | `NOTES.md` "Falsified hypotheses", honesty log |
| — | *Earlier guess (Paper B):* support confined to low exponent bits | **DISPROVED** | added qubits live in ~half the support terms — superseded by C23 | `NOTES.md` "C15 MECHANISM SOLVED", honesty note in `ABSTRACT_SHOR_2ADIC.md` |
| — | *Earlier guess (Paper B):* identity block is affine over GF(2) | **DISPROVED** | 14336/32768 violations | `NOTES.md` "Failed routes (do not retry)" |
| — | *Earlier guess (Paper B):* the (z, z⊕e) pairing is the C15 mechanism | **DISPROVED** | an n_exp=1 artifact, not r-dependent | `NOTES.md` §W5, "Failed routes (do not retry)" |
| — | *Paper A's original thesis:* "Compilation, Not Algorithm, Determines Pauli-Path Simulability" | **RETRACTED** | the apparent compilation effect was an ancilla-count confound; both compilations have exact Z-closure | `ABSTRACT.md` STATUS block; superseded v1 abstract kept there for the record |
| — | *Paper B:* an intermediate 2-adic law in α | **DISPROVED** | the N=323 reading was sampling aliasing — t = 16, 20, 24 hit residues 4, 2, 0 mod ord₂(9) = 6. The robust structure is the binary β=1 / β>1 split | `NOTES.md` §I; `experiment_c15_intermediate{,2,3}.py` (TODO step 6) |
| — | *Paper B:* the msb↔anc CNOT pairing causes the linear structure | **REFUTED** | conjugating the ancilla with extra CNOT/Toffoli couplings left defect 1 and the same w in every variant; adding XOR couplings cannot break an XOR symmetry. Real mechanism is C31 | `NOTES.md` §L2; `experiment_reduction.py` (TODO step 10) |
| — | *Paper A:* weight truncation as Fourier tail mass | **DISPROVED as stated** | PPS truncates *incrementally*, so a discarded term never branches; the Fourier tail describes **terminal** truncation only. No usable weight cutoff exists | `NOTES.md` §W (W1, W2); `experiment_weight.py`, `experiment_weight2.py` (TODO step 2) |

**C7's original framing is also dead.** Paper A's ledger carried C7 as "Results
are pre-asymptotic — **REFUTED**" (reached 24q via Walsh: density 0.473→0.498 →
½, slope 1.008 bits/qubit ⟹ Θ(2ⁿ)). The ID now carries the opposite, established
statement in the Paper B section; do not cite C7 as a pre-asymptotic caveat.

**Name collision to be aware of.** `NOTES.md` §I uses "H1" locally for a
*periodicity* hypothesis (refuted there); that is unrelated to the honesty-log
H1 above.

---

## F-number to C-number mapping

Early findings were relabelled when they were promoted to claims. The surviving
ones:

| old | now | note |
|---|---|---|
| F5 | C5 | truncation can violate \|⟨O⟩\| ≤ 1 |
| F9 | C8 | PPS term count = Walsh sparsity |
| F10 | C13 | scope of the identity: Z-type yes, X/Y no |
| F11 | C14 | δ non-monotonicity |
| F13 | C16 | heavy-tailed spectrum explains it |
| F12 | — | kept as F12 in the Paper B section |
