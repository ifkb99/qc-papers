# Draft abstract — Paper B (Shor / 2-adic)

**Working title:** *The 2-Adic Structure of the Order Determines Pauli-Path
Simulation Cost for Modular Exponentiation*

Baker, I.

---

> **STATUS — read before circulating.** Split out of the Walsh cost-model paper
> (`ABSTRACT.md`, "Paper A") because the claim, the audience and the evidence
> are separable: Paper A is a methodology result about Pauli propagation on
> permutation circuits; this one is a result about Shor's algorithm specifically,
> and it *uses* Paper A's identity as an instrument. Publishable order is A then
> B, or B alone with A's identity cited as a lemma.
>
> Post-bugfix (θ=π propagator bug, see `NOTES.md` STATUS 2) and post external
> review. Confirmed for both final support and peak memory. The core dichotomy
> is verified in a controlled design at **three** moduli (N = 5, 7, 21); the
> onset rule n_exp = v₂(r)+1 is exact for α = 1, 2 and consistent for α = 3, 4.
> **The circuit-level mechanism is now PROVED** (parity reduction via an
> involution) — see the honesty note, which also records the scope limit.

---

## Abstract (draft)

Shor's algorithm is the canonical target for classical simulation studies of
quantum advantage, yet the difficulty of simulating it is usually discussed in
terms of qubit count and circuit depth. We show that for Pauli-path simulation
(also called sparse Pauli dynamics) the dominant cost parameter is instead an
arithmetic property of the instance: the 2-adic structure of the multiplicative
order r of the chosen base.

Our instrument is the exact correspondence, established separately, between the
Pauli support carried by a Heisenberg-picture propagation through a
basis-permutation circuit and the Walsh–Hadamard spectrum of the corresponding
output-bit Boolean function. Because a^e mod N depends only on e modulo r, and
because the Walsh basis is the character group of (Z/2)^t, the arithmetic of r
relative to powers of two controls the spectrum directly. Writing r = β·2^α with
β odd, the idealised bit function e ↦ bit_j(a^e mod N) is a function of the low α
bits of the exponent alone when β = 1, giving a spectrum of *constant* size
independent of exponent-register width; any odd factor β > 1 makes the period
incommensurate with the Walsh basis and the spectrum becomes maximally spread.
We measure exactly this: at t = 24 the spectrum has 4 nonzero coefficients for
N = 15, a = 7 (r = 4), and is fully dense — density 1.000000 — for every instance
with an odd factor tested.

We then verify the effect on complete, end-to-end verified reversible modular
exponentiation circuits in a controlled design that holds the modulus, ancilla
layout, gate structure and observable fixed and varies only the base, so that r
alone changes. For N = 7 with a = 6 (r = 2) the Pauli support is *exactly
invariant* at 15549 terms, and peak memory *exactly invariant* at 24369 terms,
across a 64- to 256-fold growth in Hilbert-space dimension; for a = 3 (r = 6) the
same quantities grow by a clean factor of 4.00 per two added qubits. The result
reproduces at N = 21 and N = 5. Consequently the exponent register — which sets
the precision of the continued-fractions post-processing, and which one would
expect to be the expensive resource — is free for Pauli-path simulation when r
is a power of two and quadruples cost per two qubits otherwise.

The invariance has a precise onset. Writing α = v₂(r), the controlled-multiplier
block attached to exponent bit i multiplies by a^(2^i) mod N, which is the
identity on the valid subspace exactly when i ≥ α; the support therefore locks as
soon as the first such block appears, at **n_exp = α + 1**, and not before. We
confirm this for α = 1 and α = 2 (where the support grows once, 15493 → 32143,
before freezing) and find α = 3 and α = 4 still growing at the largest width we
can reach, as predicted. The locked value is close to half the Hilbert-space
dimension at the lock point, so the support saturates to half density and then
freezes in absolute terms while density falls fourfold per two added qubits.

The invariance is not an artifact of exact arithmetic. Applying a coefficient
threshold to the finished spectrum preserves it at *every* threshold, since the
surviving magnitudes are identical across widths; and under the incremental
truncation that Pauli-path simulators actually perform, the peak term count —
the quantity that bounds memory — is likewise unchanged at every threshold we
test. Accuracy is independent of width up to moderate thresholds but degrades
first for the wider circuit at aggressive ones, not because its exact spectrum
differs (it is identical) but because incremental truncation has more gates to
act upon. We therefore state the practical claim as: peak cost is free in the
exponent register, and accuracy is free up to moderate truncation.

For instances with an odd factor, which is the generic and cryptographically
relevant case, we find the Walsh density converges monotonically to one half
(0.473 → 0.498 over 15 to 24 qubits, growth 1.008 bits per qubit), so
Pauli-path cost is Θ(2^n): asymptotically no better than state-vector
simulation. The limiting constant of one half is not accidental and not
algorithmic: the support carries a linear structure w = b_msb ⊕ anc, pairing the
accumulator's sign bit with the modular-reduction comparison ancilla, which
confines it to a hyperplane. We verify this is observable-independent and shared
by both compilations, which use the same add/subtract/restore reduction, so the
Θ(2^n) conclusion is robust while the constant should be read as a property of
the reduction discipline rather than of modular exponentiation. The structure is
forced by three features of the standard construction acting together: flipping
the accumulator's most significant bit is the same operation as adding a power
of two and therefore commutes with the modular adder; the comparison ancilla is
coupled to that bit only through exclusive-or; and the swap network transferring
the accumulator into the multiplicand register excludes the most significant bit
entirely, so it never reaches the measured register. We confirm the diagnosis by
substituting reductions that couple the sign bit nonlinearly, which leaves the
computed function unchanged but destroys the linear structure and raises the
density from 0.473 to 0.716 — a cost increase of roughly one half arising purely
from a compilation choice, and the only such effect we observe. Reaching 24 qubits is possible only because the Walsh route computes
the exact cost without running the simulation, which stalls near 17.

Finally, we observe that the same decomposition r = β·2^α governs
matrix-product-state simulation of Shor's algorithm, where α is the number of
trailing zeros of r and the memory reduction is by a factor β² carried by the
odd part. Two structurally unrelated classical methods keying on the same
arithmetic invariant suggests a property of the algorithm rather than of either
simulator.

The practical corollary is stronger than a remark about one instance. Since every
order divides the Carmichael function λ(N), a modulus with λ(N) a power of two
lies in the free branch for *every* base — and λ(N) is a power of two exactly
when N is a power of two times a product of distinct Fermat primes. The odd
semiprimes with this property are precisely p·q with p and q both Fermat primes:
15 = 3×5, 51 = 3×17, 85 = 5×17, and so on. **N = 15 is the smallest, and all
seven of its usable bases lie in the free branch**, against 3 of 11 for N = 21.
The canonical demonstration instance is therefore degenerate for classical
simulation not because of an unlucky choice of base but because of the choice of
modulus, and the property that makes it the natural smallest example is the same
one that makes it uninformative as a benchmark.

---

## Honesty note that must survive into any submitted version

**Updated twice. The circuit-level invariance is now PROVED, not empirical.**

The original argument (r | 2^α ⟹ the bit function depends only on the low α bits
of e) is exact at *function* level but does **not** apply at circuit level: the
added exponent qubits remain live, because `u_a(ctrl, a^{2^i})` with
`a^{2^i} = 1` is the identity only on the valid subspace. That earlier claim is
withdrawn.

The correct mechanism is a **parity reduction**, and it is proved:

> `u_a(·,1) = A⁻¹SA` with S a product of disjoint transpositions, so S² = id and
> hence **V := u_a(·,1) is an involution**. All identity blocks apply this same V
> on distinct controls that V never modifies, so they commute and compose to
> `V^p` with `p = ⊕_{i≥α} e_i`. Averaging the Walsh character over that tail
> confines the support to `z_I ∈ {0, 1_I}` — two values regardless of tail
> length — so the support size carries no dependence on the number of exponent
> qubits.

This explains constancy and liveness together rather than in tension, predicts
the onset n_exp = α+1, predicts the branch relation h(y) = g(y ⊕ e_prev), and
fails for β>1 as required. Its sharpest consequence — support confined to
z_I ∈ {0, 1_I} — was derived before being tested and then confirmed with zero
violations across seven instances, with the two halves individually constant.

**Scope, now a criterion rather than a caveat.** The proof uses nothing about
the block beyond V² = id and the fact that V does not modify the controls. We
confirm this is the operative condition by substituting synthetic blocks
unrelated to modular arithmetic: an involutive one leaves the support invariant
across four widths, while an order-3 one makes it grow by roughly a factor of
two per block, at two moduli, with a vacuity check confirming both genuinely act
on the observed bit. The theorem therefore covers **any** construction whose a=1
block is an involution — a property checkable of a given construction rather
than a family one must belong to. Multiply–swap–unmultiply qualifies because
A⁻¹SA is a conjugate of a product of disjoint transpositions; whether
windowed/table-lookup arithmetic qualifies is a well-posed open question.

## Claims ledger (Paper B)

| # | Claim | Status | Evidence |
|---|---|---|---|
| **C15** | Cost is set by the 2-adic structure of r, at circuit level | **CONFIRMED, controlled** | fix N vary a: r=2 ⟹ support *exactly* 15549 over 256× dim growth; r=6 ⟹ 4.00×/step; reproduced at N=21 and N=5 |
| **C21** | Invariance has a precise onset: support locks at **n_exp = v₂(r)+1** | **established** | α=1,2 exact (α=2 grows 15493→32143 then freezes); α=3,4 still growing at max reachable width, as predicted |
| C22 | λ(N) a power of two ⟹ **every** base free; ⟺ N = 2^a × distinct Fermat primes | **established** | N=15: 7/7 bases free; N=21: 3/11. Odd semiprimes = p·q both Fermat, smallest 15 |
| **C18** | Holds for **peak memory**, not just final support | **established** | N_max exactly 24369 at n_exp = 2,4,6,8 (64× dim growth); control r=6 grows 4.02×/step |
| C7 | Generic r ⟹ Θ(2ⁿ); results are not pre-asymptotic | **established** | 24 qubits via Walsh; density 0.473→0.498→½, slope 1.008 bits/qubit |
| **C30** | The ½ ceiling is a **linear structure** w = b_msb ⊕ anc, not an algorithmic constant | **established** | GF(2) rank n−1 in every instance; g(y⊕w)=g(y) pointwise; observable-independent; present in both compilations; random f has full rank |
| **C31** | The structure is forced by (a) msb-flip commuting with mod-2^m addition, (b) XOR-only coupling to anc, (c) msb excluded from the cswaps | **established** | (a) verified exhaustively, 0/256 violations; adding further XOR couplings leaves it intact |
| **C32** | Breaking it needs **nonlinearity in the msb**, and costs ~51%: density 0.473 → 0.716 | **established** | Toffoli variants break it (full rank); cswap variant does not (linear); all four still compute a^e mod N |
| **C27** | Invariance survives truncation: exact for terminal thresholding at every δ; peak cost unchanged under incremental | **established** | counts identical across n_exp at all 7 δ values; N_max SAME at every δ; control diverges |
| C28 | Accuracy degrades first for the *wider* circuit at aggressive δ | **established, must be stated** | ⟨O⟩ → 0 at δ=1e-1 for n_exp=4,5 while n_exp=3 stays exact; more gates ⟹ more incremental truncation |
| F12 | Function-level dichotomy is absolute | **established** | r=4: sparsity 4 constant to t=24; odd factor: density 1.000000 |
| C19 | Same r = β·2^α invariant governs MPS simulation | **established, cited** | Dang et al. §4: α = trailing zeros, β = odd part "cannot be localised"; §5.2: memory ∝ β² |
| C20 | N=15 is a degenerate benchmark **for every base**, forced by the modulus | **established, strengthened** | 100% of its bases free; smallest product of two Fermat primes |
| **C23** | Mechanism: identity blocks apply one involution V, so the circuit depends on the identity tail only through a **parity bit** | **PROVED** | V²=id since V=A⁻¹SA with S disjoint transpositions; blocks commute; character averaging confines support to z_I ∈ {0,1_I} |
| **C24** | Consequence: support confined to z_I ∈ {0, all-ones}, so size is independent of tail length | **PROVED + verified** | derived before testing; 0 violations in 7/7 β=1 instances, halves individually constant; controls fail with 48189 / 1556046 |
| **C29** | The theorem needs only **V²=id**; it covers any construction with an involutive repeated block | **PROVED + verified** | synthetic involution constant across k=1..4 at two moduli; synthetic order-3 block grows ~2× per block; vacuity check passes |
| — | *Earlier guess:* support confined to low exponent bits | **DISPROVED** | added qubits live in ~half the support terms — superseded by C23 |
| — | *Earlier guess:* identity block is affine over GF(2) | **DISPROVED** | 14336/32768 violations |

## Dependencies on Paper A

- The Walsh identity (Paper A, C8) is used as the measurement instrument.
- `perm_pps.py` (Paper A, C17) supplies the peak-memory numbers; it is what makes
  the n_exp = 8 peak measurement tractable (7.4 s at q = 21 versus rotation-level
  propagation failing to finish q = 17 in 20 minutes).

## Open problems

1. **Prove the circuit-level invariance.** The valid-subspace argument does not
   cover it. This is the paper's main gap.
2. **Intermediate 2-adic structure.** N = 323, r = 144 = 16·9 shows partial
   sparsity at small t (density 0.981 at t = 16) washing out to 1.000 by t = 24.
   Is there a quantitative law in α versus t, rather than the current binary
   split?
3. **Does the criterion survive noise/truncation?** All figures here are δ = 0.
   Paper A's C14 shows truncation behaves non-monotonically on these circuits.
4. Whether the MPS correspondence extends to other methods (decision diagrams,
   tensor-network contraction orders) or is specific to these two.

## Literature anchors

- Dang, Hill & Hollenberg, [arXiv:1712.07311](https://arxiv.org/abs/1712.07311) —
  MPS Shor, 60 qubits. §4 gives α = trailing zeros of r, β = r/2^α odd part;
  §5.2 gives the β² memory factor. **The 2-adic decomposition is explicit in the
  body even though the abstract says only "its factors".**
- Beauregard, [quant-ph/0205095](https://arxiv.org/abs/quant-ph/0205095) —
  Fourier modexp construction.
- Cuccaro et al., [quant-ph/0410184](https://arxiv.org/abs/quant-ph/0410184) —
  ripple-carry adder used for the Toffoli compilation.
- Gharibyan et al., [arXiv:2507.10771](https://arxiv.org/pdf/2507.10771) — the
  PPS resource framework whose extrapolation machinery this bypasses in this
  regime.

## Reproduction

```
uv run python experiment_c15.py     # controlled: fix N, vary a  (N=7 and N=21)
uv run python experiment_c15b.py    # exact constancy sweep, n_exp 2..10
uv run python experiment_c7.py      # scaling to 24 qubits, density -> 1/2
```
