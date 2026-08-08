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
> is verified in a controlled design at two moduli; the *mechanism* at circuit
> level is empirical rather than proved — see the honesty note below, which must
> survive into any submitted version.

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
reproduces at N = 21. Consequently the exponent register — which sets the
precision of the continued-fractions post-processing, and which one would expect
to be the expensive resource — is free for Pauli-path simulation when r is a
power of two and quadruples cost per two qubits otherwise.

For instances with an odd factor, which is the generic and cryptographically
relevant case, we find the Walsh density converges monotonically to one half
(0.473 → 0.498 over 15 to 24 qubits, growth 1.008 bits per qubit), so
Pauli-path cost is Θ(2^n): asymptotically no better than state-vector
simulation. Reaching 24 qubits is possible only because the Walsh route computes
the exact cost without running the simulation, which stalls near 17.

Finally, we observe that the same decomposition r = β·2^α governs
matrix-product-state simulation of Shor's algorithm, where α is the number of
trailing zeros of r and the memory reduction is by a factor β² carried by the
odd part. Two structurally unrelated classical methods keying on the same
arithmetic invariant suggests a property of the algorithm rather than of either
simulator. A practical corollary is that the standard N = 15, a = 7
demonstration instance (r = 4) lies in the degenerate branch of both methods,
and is therefore a poor benchmark for classical-simulation difficulty.

---

## Honesty note that must survive into any submitted version

At **function** level the mechanism is exact and provable: r | 2^α implies the
bit function depends only on the low α bits of e, bounding the Walsh support by
2^α independently of register width.

At **circuit** level the invariance is verified but *not* explained by that
argument. Direct measurement shows the added exponent qubits remain **live** —
for r = 2, n_exp = 6, qubits 15–18 are each set in ≈7779 of the 15549 support
terms. The support does not shrink into a subspace; it is relabelled while its
cardinality is preserved. The cause is that `u_a(ctrl, a^{2^i})` with
a^{2^i} = 1 is the identity only on the *valid* subspace; as a full-space unitary
it acts nontrivially on invalid inputs, and the constancy arises through
cancellation there. **Report the circuit-level result as an empirical regularity
across five widths and two moduli, not as a corollary of the subspace argument.**
Proving it is the main open problem for this paper.

---

## Claims ledger (Paper B)

| # | Claim | Status | Evidence |
|---|---|---|---|
| **C15** | Cost is set by the 2-adic structure of r, at circuit level | **CONFIRMED, controlled** | fix N vary a: r=2 ⟹ support *exactly* 15549 over 256× dim growth; r=6 ⟹ 4.00×/step; reproduced at N=21 |
| **C18** | Holds for **peak memory**, not just final support | **established** | N_max exactly 24369 at n_exp = 2,4,6,8 (64× dim growth); control r=6 grows 4.02×/step |
| C7 | Generic r ⟹ Θ(2ⁿ); results are not pre-asymptotic | **established** | 24 qubits via Walsh; density 0.473→0.498→½, slope 1.008 bits/qubit |
| F12 | Function-level dichotomy is absolute | **established** | r=4: sparsity 4 constant to t=24; odd factor: density 1.000000 |
| C19 | Same r = β·2^α invariant governs MPS simulation | **established, cited** | Dang et al. §4: α = trailing zeros, β = odd part "cannot be localised"; §5.2: memory ∝ β² |
| C20 | N=15 a=7 (r=4) is a degenerate benchmark | **established** | lies in the free branch of both PPS and MPS |
| — | *Circuit-level mechanism* (support confined to low exponent bits) | **DISPROVED as stated** | added qubits are live in ~half the support terms; see honesty note |

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
