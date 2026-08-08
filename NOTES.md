# PPS / Shor working notes

Working notes, not a report. Dense on purpose. Written to restore context fast.

**Question being probed:** does the Pauli Path Simulation (PPS) resource-prediction
framework of Gharibyan et al. ([arXiv:2507.10771](https://arxiv.org/pdf/2507.10771))
transfer from the brickwork/kicked-Ising circuits it was derived on to
*structured arithmetic* circuits (Shor-like)?

**Standing constraint that bounds the whole project:** efficient classical
simulation of Shor's for arbitrary N ⟹ classical poly-time factoring. So nothing
here can be a route to "simulate Shor's at scale." The only live question is
*characterising where PPS breaks*, which is a tooling/limits result, not a
speedup result. Do not lose sight of this.

---

## Conventions (get these wrong and everything silently breaks)

**Pauli representation** (`pauli.py`): string = `(x, z)` bitmask pair,

```
P(x,z) = i^{popcount(x & z)} * X^x Z^z
```

The `i^{|x&z|}` is what makes P Hermitian (turns XZ into Y). Consequence:
**all coefficients in the Heisenberg expansion stay real.** Qubit 0 is the
least significant bit; `to_matrix` builds with `np.kron(m, M)` in that order.

`pauli_mult(p,q) -> (r, k)` with `P(p)P(q) == i^k P(r)`, where
`k = |a&b| + |x&z| - |rx&rz| + 2*|b&x|  (mod 4)`. Verified exhaustively n=3.

**Everything is a Pauli rotation.** `U = exp(-i θ σ / 2)`, one conjugation rule:

```
U† P U = P                            if [P,σ]=0
       = cos(θ) P + sin(θ) (i σ P)     if {P,σ}=0
```

Clifford = θ=±π/2 ⟹ cos=0 ⟹ maps to a *single* Pauli, no branching.
All branching comes from non-Clifford angles. This is the whole cost mechanism.

**Gate decompositions** (`circuits.py`, all verified up to global phase):
- `H = Rz(π/2) Rx(π/2) Rz(π/2)`
- `X = Rx(π)`, `T = Rz(π/4)`, `S = Rz(π/2)`
- `CNOT(c,t) = Rzx(-π/2) · Rz_c(π/2) · Rx_t(π/2)`, σ for the ZX term is `(1<<t, 1<<c)`
- `Toffoli` = standard 15-gate Clifford+T, **exactly 7 T gates** (asserted in tests)
- `CP(θ) = Rz_a(θ/2) Rz_b(θ/2) Rzz(-θ/2)`

**Ordering.** `Circuit.gates[0]` is applied to the state *first*. PPS therefore
iterates `reversed(circuit.gates)` (Heisenberg: innermost conjugation is the
last gate). Getting this backwards gives plausible-looking wrong answers.

**Expectation readout.** `<0|P(x,z)|0>` = 1 if `x == 0`, else 0. So
`<O> = sum of coeffs over Z-type terms`.

---

## Verified infrastructure

`test_core.py` — all pass. Do not trust any result if these regress.

```
[1] pauli_mult / i·σ·P / commutation vs dense, n=3, all 4096 pairs
[2] H, X, T, S, CNOT(0,1), CNOT(1,0), Toffoli decompositions; Toffoli T-count == 7
[3] Cuccaro ripple adder correct on all inputs, nbits=2,3; T-count == 14·nbits
[4] PPS(δ=0) == dense expectation, random circuits (max err ~1e-9)
[5] PPS(δ=0) == dense on a real adder circuit
```
`experiment2.py` §0 additionally verifies `CP`, `QFT` vs DFT, `QFT∘QFT⁻¹ = I`.

---

## STATUS 2 — CRITICAL BUG FOUND; F4/F6/F7/F8 ALL RETRACTED

**A bug in `pps.py` invalidated every result involving an X, Y or Z gate.**

```python
if abs(s) < 1e-12:      # WRONG: "identity up to phase"
    new = terms
```

`U†PU = cos(θ)P + sin(θ)(iσP)` for anticommuting P. At **θ = π**, `sin = 0` but
`cos = −1`, so anticommuting terms must be **negated**. The shortcut treated
every θ=π gate — i.e. every X, Y, Z — as a no-op. Fixed: the fast path now
requires `cos > 0`.

**Why the test suite missed it:** the random circuits in `test_core.py` [4] drew
from `{h, t, cnot, rx, rz}` and never emitted an `x()`. The ripple adder is built
only from CNOT/Toffoli — also no X. The bug could only fire in `ToffoliModExp`
(`_load` uncontrolled, `cc_add_mod`'s `qc.x(msb)`) and in `build()`'s `x(x[0])`.
Regression test added as `test_core.py` [4b], which fails on the old code.

**How it was caught:** Walsh (exact integer arithmetic) disagreed with PPS. The
discriminating test was running PPS in `longdouble` — the error was *identical*
to `float64` (1.000e+00), ruling out precision and proving a logic bug.
**Keep that trick: if extra precision doesn't move the error, it isn't precision.**

### Corrected results (all re-run post-fix)

```
                       exact PPS (delta=0)
 N,a  compilation   q  gates    N_max   N_fin  nonZ      <O>  true
 5,2      Fourier  10   5359     2530     451     0  -1.0000    -1
 5,2      Toffoli  15  21247    40770   15493     0  -1.0000    -1
 7,3      Fourier  10   5359     2494     461     0  -1.0000    -1
 7,3      Toffoli  15  21715    40796   15539     0  -1.0000    -1
```

Every prior claim that flipped:

- **F4 RETRACTED.** "Fourier lacks Z-closure" was pure artifact. **Both**
  compilations end with **zero** non-Z terms. Z-closure is a property of the
  *unitary being a basis permutation*, not of the gate set.
- **F6 RETRACTED.** The direction reverses: Fourier is *cheaper* (451 vs 15493),
  not 8.5x more expensive.
- **F7 SUPERSEDED by F9** (below) — right instinct, wrong framing.
- **F8 RETRACTED.** ⟨O⟩ was wrong (+1 where truth is −1), so the δ-sweep
  conclusions were meaningless.
- **C5 partially survives**: norm violation still occurs under truncation
  (Fourier N=7 δ=1e-3 gives ⟨O⟩ = −1.02086), but it is now the *only*
  truncation claim standing.

### CONFOUND that kills the A/B as designed

The two compilations use **different qubit counts** (10 vs 15) with different
ancilla layouts. They agree on the *valid subspace* but implement **different
permutations of their respective full Hilbert spaces**, so their pullbacks are
different operators over different-sized domains. 451 vs 15493 is therefore
**not** a compilation effect — it is mostly an ancilla-count effect.
**Any future A/B must match total qubit count.**

---

## F9 — THE RESULT: PPS term count = Walsh sparsity, exactly

For a circuit implementing a basis permutation π, `π†Z_jπ` is the diagonal
operator `(−1)^{g(y)}` with `g(y) = bit j of π(y)`. Expanding a diagonal operator
in the Pauli basis is *precisely* the Walsh–Hadamard transform of `(−1)^g`.
Therefore:

> **the number of Pauli terms PPS carries = the Walsh sparsity of g**

Verified to machine precision, 6/6 instances, both compilations, supports
identical (not merely counts):

```
 N,a  compilation   q   walsh     pps  match    maxerr
 5,2      Fourier  10     451     451    YES  6.66e-16
 5,2      Toffoli  15   15493   15493    YES  3.75e-16
 7,3      Fourier  10     461     461    YES  3.61e-16
 7,3      Toffoli  15   15539   15539    YES  4.44e-16
 4-bit adder Z_b0  10       1       1    YES  1.11e-16
 4-bit adder Z_b2  10      10      10    YES  1.11e-16
```

Consequences:

1. **Exact predictive cost model.** Walsh runs in O(2ⁿ·n) and took 0.05s where
   PPS took 49s. No extrapolation, no power-law fitting, no truncation
   heuristics — for permutation circuits the answer is computable outright.
   This sidesteps the paper's Eq. 17 machinery entirely *in this regime*.
2. **Compilation-invariant.** Kills the "compilation, not algorithm" thesis.
3. **Explains F1 correctly at last.** Adder low bit = `a0⊕b0⊕c0`, affine ⟹ Walsh
   sparsity 1. Sparsity 1 ⟺ affine is a standard theorem. The permutation
   property was never the operative fact.
4. **Bridge to cryptanalysis.** Walsh sparsity / linearity is *the* central
   quantity in linear cryptanalysis. "PPS-hard reversible circuit" ≈ "Boolean
   function resistant to linear approximation". That literature is deep and
   directly importable — likely the most valuable thread here.

**Caveat:** Walsh is itself O(2ⁿ), so this predicts cost rather than beating it.
Its value is as an *exact* cost model and an explanation, not a faster simulator.

## REVIEW ROUND (external) — 3 of 4 findings upheld, 1 refuted

An independent review raised four objections. Verified each rather than
accepting; results below. **Net: the core (F9/C8) is unharmed, C15 is
strengthened, one real gap found and fixed, one review claim was wrong.**

### R1 UPHELD (important) — Walsh predicts N_final, not N_max (peak memory)

PPS memory cost is the *peak* term count during propagation, not the final one.
Walsh sparsity equals the **final** count. These differ:

```
        circuit    walsh/final    rot N_max   ratio
 3-bit adder Z_b0            1          128    128x
 4-bit adder Z_b0            1          512    512x
 5-bit adder Z_b0            1         2048   2048x
      modexp N=5         3086        13666    4.4x
     modexp N=15        31176       128138    4.1x
```

**So F1's "the adder is FREE for PPS" is misleading** — it peaks at the full
lightcone (2^(q-1)) before collapsing to 1. And C11's "exact cost model" predicts
*result complexity*, not peak memory. Both need restating. For modexp the ratio
is a bounded ~4x, so **C7's Θ(2ⁿ) asymptotic is unaffected**.

### R1 FIX — permutation-native PPS (`perm_pps.py`), a real improvement

Root cause: standard PPS decomposes Toffoli into Clifford+T and propagates
through the *rotations*, which include Hadamards — so the operator leaves the
diagonal mid-circuit even though the Toffoli as a whole is a permutation.

Treating X/CNOT/Toffoli as **atomic permutation gates** keeps everything Z-type
at every step. Conjugation rules (from `(-1)^{popcount(perm(y) & z)}`):

```
X(q)           : Z^z -> (-1)^{z_q} Z^z                        1 term
CNOT(c,t)      : Z^z -> Z^{z XOR (z_t << c)}                  1 term
Toffoli(a,b,c) : z_c=0 -> Z^z                                 1 term
                 z_c=1 -> 1/2[Z^z + Z^{z^a} + Z^{z^b} - Z^{z^a^b}]   4 terms
```

Measured (verified against Walsh and against rotation-level PPS,
`test_perm_pps.py` all pass):

```
        circuit    final   rot N_max   perm N_max   saving
 3-bit adder Z_b0      1         128           64     2.0x
      modexp N=5    3086       13666         6834     2.0x
     modexp N=15   31176      128138        64070     2.0x
```

> **CORRECTED 2026-08-08 — "exactly 2.0x everywhere" is a ROUNDING ARTIFACT of
> the `%.1f` in the table above.** The true ratios are **2.000000** for the
> adders (128/64, 512/256, 2048/1024) but **1.9997** for modexp:
> 13666/6834, 16386/8194, 128138/64070. Every modexp instance satisfies
> `rot = 2·perm − 2` exactly, 3/3 — a clean regularity we have **not**
> explained. The upper bound of 2 is arguable from the gadget's H exchanging
> the Z- and X-sectors (strings containing Z_c get a mirrored partner, strings
> without Z_c do not), which also explains why the adders attain it and modexp
> does not; the constant deficit of 2 is unexplained. Claim C17 regraded: the
> factor is **empirical, not proved**. Caught by Fable reviewing `PAPER_A.md`
> — the abstract said "halves peak memory exactly" while §5 supported only
> "2.0× on all 6 instances", which are different epistemic states.

At the peak the rotation-level run holds Z and non-Z terms in near-equal
measure, so dropping the non-Z part is close to a factor of two. Modest, but two
things matter more than the constant:

1. The peak becomes a **well-defined Walsh quantity**: after k gates the term
   count *is* the Walsh sparsity of the k-gate suffix's pullback, so
   `N_max = max over suffixes of Walsh sparsity`. Exact, not extrapolated.
2. It is **much faster** — logical ops instead of rotations (1699 vs 21247 for
   N=5) and half the terms. Did q=21 in 7.4s where rotation-level PPS failed to
   finish q=17 in 20 minutes.

Actionable claim for a PPS engine: *implement reversible arithmetic blocks as
atomic permutation primitives; the propagation never leaves the diagonal.*

### R1 consequence — C15 SURVIVES for peak memory, not just final count

The critical test the review did not run. Using `perm_pps`:

```
   a   r  n_exp    q   N_final  N_max(peak)  peak growth
   6   2      2   15     15549        24369            -
   6   2      4   17     15549        24369        1.00x
   6   2      6   19     15549        24369        1.00x
   6   2      8   21     15549        24369        1.00x

   3   6      2   15     15539        24412            -
   3   6      4   17     64353        98018        4.02x
```

Peak is **exactly** constant at 24369 across a 64x growth in dimension for r=2,
against 4.02x/step for r=6. **C15's headline holds for memory cost, not merely
result size.** (Also confirmed at rotation level: N_max 40753 at both n_exp=2
and n_exp=3.)

### R2 UPHELD — the C15 *mechanism* stated in ABSTRACT.md is wrong

Checked directly: for r=2, n_exp=6, the added exponent qubits (15–18) are each
set in ~7779/15549 support terms — **genuinely live**, not confined. So the
abstract's "the Walsh support is confined to the low bits of the exponent
register" is false at circuit level. It is true for the *idealised* function
(F12), where the domain is e alone. Cause: `u_a(ctrl, 1)` is identity only on
the valid subspace; as a full-space unitary it entangles the control with ancilla
garbage. The count is constant; the support is not confined.

NOTES already carried this caveat; **the abstract did not, and has been fixed.**

### R3 REFUTED — the "2-adic vs factors of r" objection was based on the abstract only

The review fetched only Dang et al.'s abstract ("depends on its factors") and
concluded MPS keys on general smoothness while PPS keys on powers of two,
offering r=6 as a counterexample. **The full text says otherwise.** Section 4:

> "α is the number of trailing zeroes in the binary representation of r" …
> "due to the **odd factor β ≡ r/2^α** of r which cannot be localised to
> specific qubits"

and §5.2: MPS matrices for the qubits in A "reduced by a factor of **β²**".

So Dang et al. decompose r = β·2^α *exactly* as 2-adic valuation times odd part,
and it is the **odd part β that costs memory** — β=1 (r a power of two) is the
free case. r=6 has β=3, so it is in the expensive branch for MPS too, and is
**not** a counterexample. The correspondence is real and now citable with a
specific mechanism. Claim strengthened rather than weakened.

### R4 UPHELD in substance, but its numbers are confounded

True and useful: at rotation level the Toffoli decomposition passes through H
gates, so the operator does leave the diagonal mid-circuit — "Z-closure" is a
property of the *endpoint*, not the trajectory. That is exactly R1, and
`perm_pps.py` fixes it.

But the specific figures offered (Toffoli peak non-Z 13664 vs Fourier 848) compare
**q=14 against q=9** — the same ancilla/qubit-count confound already logged as
killing C3. Not a valid compilation comparison.

---

## C15 MECHANISM SOLVED (TODO step 3) — the parity reduction

Paper B's main open problem. Two routes failed first; the third worked and every
link is now verified.

### Failed routes (do not retry)

1. **Affineness of the identity block.** If `u_a(ctrl,1)` were affine over GF(2),
   conjugation would relabel the Walsh spectrum and preserve sparsity. **It is
   not** — 14336 violations of 32768, and conjugating `Z_j` by it gives 2504
   support terms rather than 1. (`experiment_c15_proof.py`)
2. **The (z, z⊕e) pairing.** Exact at n_exp=1 for *both* r=2 and r=6, broken at
   n_exp≥2 for both. An n_exp=1 artifact. (§W5)

### The mechanism

Adding an exponent qubit appends a block controlled on it, so the pulled-back
bit function splits into two branches — g (control off) and h (control on) —
with Walsh coefficients

```
    c_(z,0) = ( ĝ(z) + ĥ(z) ) / 2        c_(z,1) = ( ĝ(z) − ĥ(z) ) / 2
```

Size is preserved *with the new qubit live* precisely when **|ĥ(z)| = |ĝ(z)|**
pointwise: then exactly one of each pair survives. Verified
(`experiment_c15_proof2.py`):

```
                                    BOTH  EXACTLY-ONE   M2 magnitudes
 N=7 a=6 (r=2)  n_exp=3               0        15549   match, err 0.00e+00
 N=7 a=6 (r=2)  n_exp=4               0        15549   match, err 0.00e+00
 N=5 a=2 (r=4)  n_exp=4               0        32143   match, err 0.00e+00
 N=7 a=3 (r=6)  n_exp=3           14934          844   FAILS
 N=21 a=4 (r=3) n_exp=3         1036508         1878   FAILS
```

Bit-for-bit magnitude agreement for β=1; decisive failure for β>1.

The sign pattern `ĥ(z) = ε(z)ĝ(z)` is then a **linear character**
(`experiment_c15_proof3.py`): multiplicative on 5000/5000 sampled triples, and
`ε(z) = (−1)^⟨z,v⟩` holds on **every** support element (15549/15549,
32143/32143, 15509/15509). Moreover **v is a single bit — always the previous
exponent qubit.** So `h(y) = g(y ⊕ e_prev)`: a translation, which preserves
Walsh magnitudes exactly.

A single-bit shift equal to another control is the fingerprint of a **parity
dependence**, and that is what it is (`experiment_c15_proof4.py`):

- **V1.** For i ≥ α, `a^(2^i) mod N = 1`, so every such block applies the *same*
  fixed permutation `V = u_a(·,1)`, controlled on its own qubit. **V is an
  involution**: `V² = id` exactly, 32768/32768 fixed points (2097152/2097152 at
  N=21). Consistent with `V² = u_{a²}`, so `u_a(·,a)` is an involution iff
  `a² ≡ 1` — which the data confirms (true for a=4 mod 5, a=6 mod 7, a=8 mod 21;
  false for a=2 mod 5, a=3 mod 7).
- **V2.** Consequently the composite depends on those controls only through
  `p = ⊕_{i≥α} e_i`. Verified directly: flipping any **two** identity-block
  controls together leaves the pulled-back bit function **pointwise unchanged**.
  True for N=7 a=6, N=5 a=4, and N=5 a=2 (α=2).

### The chain, and what it explains

> For β=1, blocks with i ≥ α all apply the same involution V. The circuit
> therefore depends on the entire identity tail through **one** parity bit,
> however many qubits it spans.

| observation | explained by |
|---|---|
| support size constant in n_exp | one effective variable regardless of tail length |
| **added qubits stay LIVE** | the parity involves every one of them |
| `h(y) = g(y ⊕ e_prev)` | toggling any one control just flips the parity |
| locks at n_exp = α+1 (C21) | the first identity block sits at i = α |
| grows for β>1 | no block is the identity, so V differs per block |

The "added qubits are live" puzzle raised in the external review — which
defeated the naive valid-subspace argument — is resolved: liveness and constancy
are *both* consequences of parity dependence, not in tension.

### FORMALISED — this is now a theorem, not an empirical regularity

Both owed links closed, and they were three lines each once the circuit
structure was looked at properly.

**Setup.** `u_a(ctrl, a)` is built as *multiply, swap, unmultiply*:
```
u_a(ctrl,a) = cmult_mod(ctrl,a) ; cswap layer ; cmult_mod(ctrl,a⁻¹)⁻¹
```
For a = 1 the two multiplies are inverse to each other, so as an operator
**V := u_a(·,1) = A⁻¹ ∘ S ∘ A**, with A the controlled multiply-accumulate and S
the cswap layer.

**(i) V² = id.** S is a product of `cswap(ctrl, x_i, b_i)` over i; the target
pairs (x_i, b_i) are disjoint, so the factors commute and each is an involution,
giving S² = id. Hence
`V² = A⁻¹SA·A⁻¹SA = A⁻¹S²A = A⁻¹A = id`. **V is an involution because it is a
conjugate of one.** (Was numerical; now proved.)

**(ii) The identity blocks commute.** V never modifies an exponent qubit — they
appear only as controls — and all identity blocks apply the same V, so on any
basis state `C_{e_i}(V) C_{e_j}(V)` acts as `V^{e_i+e_j}` either way.

**(iii) Parity reduction.** With (i) and (ii),
`∏_{i≥α} C_{e_i}(V) = V^{Σ e_i} = V^{p}` where `p = ⊕_{i≥α} e_i`. The circuit
depends on the whole identity tail through that one bit. (Was numerical.)

**(iv) Support-size independence.** Write `f(y, e_I) = F(y, p)` for the tail
I = {i ≥ α}. Averaging the Walsh character over e_I, and using that the
annihilator of the even-parity subgroup is exactly {0, 1_I}:

```
   z_I = 0        ->  ( (-1)^F(y,0) + (-1)^F(y,1) ) / 2
   z_I = all-ones ->  ( (-1)^F(y,0) - (-1)^F(y,1) ) / 2
   otherwise      ->  0        (character nontrivial on the even-parity subgroup)
```

So the support is confined to **z_I ∈ {0, 1_I}** — two values, whatever |I| is —
and `|support| = #{z_y : c(z_y,0)≠0} + #{z_y : c(z_y,1_I)≠0}`, which contains no
dependence on |I|. ∎

**Verification of (iv)'s sharpest consequence** (`experiment_c15_proof5.py`),
which had not been tested before deriving it:

```
   N   a   r  al  n_exp  |I|  |support|     z_I=0  z_I=all1   other  holds
   7   6   2   1      3    2      15549      7770      7779       0   True
   7   6   2   1      4    3      15549      7770      7779       0   True
   7   6   2   1      5    4      15549      7770      7779       0   True
   5   2   4   2      4    2      32143     16089     16054       0   True
   5   2   4   2      5    3      32143     16089     16054       0   True
  21   8   2   1      3    2    1037174    518574    518600       0   True
CONTROL (beta>1):
   7   3   6          4          64353      8075      8089    48189  False
  21   4   3          3        2074894    259435    259413 1556046  False
```

`other = 0` in every β=1 case, and the **two halves are individually constant**
across |I| = 2, 3, 4 — not merely the total. Controls fail massively.

**Scope.** The proof uses only that the block has the multiply-swap-unmultiply
form, so it covers **both compilations here** (the Fourier `ModExp` builds `u_a`
identically) and any Vedral/Beauregard-style construction. It does not
automatically transfer to a modular exponentiation built some other way.

**Status: C15 is a theorem for this circuit family.** What remains for Paper B is
presentation, not proof.

## C15 REFINED (TODO step 4) — it locks at n_exp = v2(r)+1, not "always"

Sweeping a third modulus caught an over-claim. N=5 with a=2 has r=4 (α=2), where
every previous sweep used r=2 (α=1). It is **not** constant from n_exp=2:

```
   N   a   r  alpha  n_exp    q  |support|  locked?
   5   4   2      1      1   14       2708        -
   5   4   2      1      2   15      15509     grew
   5   4   2      1      3   16      15509   LOCKED
   5   4   2      1      5   18      15509   LOCKED

   5   2   4      2      1   14       3086        -
   5   2   4      2      2   15      15493     grew
   5   2   4      2      3   16      32143     grew
   5   2   4      2      4   17      32143   LOCKED
   5   2   4      2      5   18      32143   LOCKED

  17   2   8      3      2   21    1037405     grew
  17   2   8      3      3   22    2093137     grew     (lock predicted at 4)
  17   3  16      4      3   22    2093202     grew     (lock predicted at 5)
```

**Rule: the support locks at n_exp = α + 1, where α = v2(r).** Verified exactly
for α=1 and α=2; α=3 and α=4 are still growing at the largest width reachable
(q ≤ 22), consistent with locking at 4 and 5.

**Mechanism (this one *is* clean).** The u_a block for exponent bit i multiplies
by `a^(2^i) mod N`, which equals 1 exactly when r | 2^i, i.e. when i ≥ α. So
blocks i ≥ α are identity on the valid subspace. With n_exp qubits the blocks
are i = 0..n_exp−1, so at least one identity block exists iff n_exp ≥ α + 1.
The support locks the moment the first identity block appears.

**Why this was missed:** every earlier sweep used r=2, i.e. α=1, and started at
n_exp=2 = α+1 — exactly on the threshold. Pure luck. Had the original controlled
design used r=4, the first two data points would have disagreed and the claim
would have looked false.

**Corrected statement for Paper B:**

> support is constant in n_exp **for n_exp ≥ v2(r) + 1** when r is a power of
> two, and Θ(2^q) as soon as r has an odd factor.

Secondary observation: the locked value is close to **half** the Hilbert space
at the lock point (density 0.473 at α=1, 0.490 at α=2, and 0.499 already at
n_exp=3 for α=3). So the support grows to half-density, then freezes in absolute
terms while density falls by 4× per two added qubits.

## C20 STRENGTHENED — N=15 is degenerate for EVERY base, and necessarily so

Every order divides the Carmichael function λ(N), so if λ(N) is a power of two
then **every** base has β=1 and the whole modulus is in the free branch.

```
    N   factors  lambda(N)  pow2?          orders present
    5         5          4    YES               [1, 2, 4]
   15       3x5          4    YES               [1, 2, 4]
   21       3x7          6     no            [1, 2, 3, 6]
   33      3x11         10     no           [1, 2, 5, 10]
   51      3x17         16    YES        [1, 2, 4, 8, 16]
   85      5x17         16    YES        [1, 2, 4, 8, 16]
  143     11x13         60     no   [1,2,3,4,5,6,10,12,15,20,30,60]
```

```
  N=15: 7 usable bases, 7 with beta=1  -> 100%
  N=21: 11 usable bases, 3 with beta=1 ->  27%
```

λ(N) is a power of two exactly when **N = 2^a × (product of distinct Fermat
primes)** — since p−1 must be a power of two for each odd prime p, and p^k needs
k=1. Known Fermat primes: 3, 5, 17, 257, 65537.

**So the odd semiprimes in the free branch are exactly p·q with both p and q
Fermat primes: 15 = 3×5, 51 = 3×17, 85 = 5×17, … and 15 is the smallest.**

That is a sharper version of C20 than "N=15, a=7 happens to have r=4". The
canonical demonstration instance is degenerate **for every base**, and it is
degenerate *because* it is the smallest product of two Fermat primes — the same
property that makes it the natural smallest demo. The degeneracy is forced by
the choice of N, not by the choice of a.

## L — THE DENSITY-½ CEILING EXPLAINED: a linear structure from the reduction ancilla

Two observations had gone unexplained across the whole project:

- **C7:** circuit density converges to 0.498 — approaching ½ **from below** and
  never crossing it — where random Boolean functions give 1.000.
- **Step 6:** the real modexp table at r=6 gives density exactly 0.500 on even t
  where a random table of the same period gives 1.000.

Approaching ½ from below without crossing is the signature of a **linear
constraint on the support**. Confirmed:

```
                                |supp|  density   GF(2) rank  defect
  modexp N=5 a=2 n_exp=1          3086  0.188354      13/14        1
  modexp N=7 a=6 n_exp=2         15549  0.474518      14/15        1
  modexp N=15 a=7 n_exp=1        31176  0.237854      16/17        1
  random f, n=14                 16384  1.000000      14/14        0
```

The support has GF(2) rank **n−1**, so it lies in a hyperplane and the density
is capped at exactly ½. Dually, that means a **linear structure**: a nonzero w
with `g(y ⊕ w) = g(y)` for all y — verified pointwise.

### What w is

**w = b_msb ⊕ anc** — the accumulator's sign bit XOR the modular-reduction
comparison ancilla. `{b3, anc12}` at N=5,7; `{b4, anc15}` at N=15.

The pairing is visible directly in the source. `cc_add_mod` does

```python
qc.cnot(msb, self.anc)                        # set the negative flag
...
qc.x(msb); qc.cnot(msb, self.anc); qc.x(msb)  # uncompute it
```

so flipping msb and anc *together* is a symmetry of the reduction step.

### Scope: circuit property, not observable, and shared by both compilations

```
  Toffoli, N=7:  Z_x0 Z_x1 Z_x2 Z_b0  all give defect 1, w = b3,anc12
  Fourier, N=5:  defect 1, w = b3,anc7      (its own sign bit + ancilla)
  Fourier, N=7:  defect 1, w = b3,anc7
```

**Observable-independent** and present in **both** compilations, which share the
add / subtract-N / conditional-restore reduction. So it follows from the
reduction *discipline*, not from one implementation — but it is not a property
of modular exponentiation as an algorithm. A reduction that computes its
comparison flag differently might not have it.

### Three consequences

1. **C7's ½ is explained exactly** — the support is filling a hyperplane whose
   density is ½ by construction.
2. **C7's Θ(2ⁿ) conclusion is unaffected** (½·2ⁿ is still Θ(2ⁿ)), but **the
   constant ½ must be reported as construction-dependent**, not algorithmic.
   This is a correction to how C7 is currently phrased.
3. **Feeds back into C12.** Bent functions have **no** linear structures at all,
   so possessing one bounds modexp away from the bent bound automatically. That
   partly explains why modexp sits at 0.74 of the bent bound where random
   functions reach 0.97–0.99 — it is not merely "less random", it carries a
   specific, standard cryptanalytic weakness.

The ripple adder shows the same phenomenon more strongly: kernel dimension 5,
so density ≤ 2⁻⁵ = 0.031 (measured 0.0098).

## L2 — STEP 10 RESOLVED: the standard reduction is *accidentally* PPS-friendly

First hypothesis (that the msb↔anc CNOT pairing causes it) was **refuted**:
conjugating the ancilla with extra couplings — CNOT from t0, Toffoli from t0&t1,
CNOT from b0 — left the structure completely unchanged, defect 1 with the *same*
w in every case. Adding XOR couplings cannot break an XOR symmetry.

### Why the structure is robust — three ingredients

**(a) An msb flip commutes with the arithmetic.** Flipping the msb of an m-bit
register *is* adding 2^(m−1), and mod-2^m addition is commutative:
`(b ⊕ 2^(m−1)) + c == (b + c) ⊕ 2^(m−1)`. Verified exhaustively — **0 violations
over all 256 (b,c) pairs at m=4**. So the flip propagates through every
`_add_const` untouched.

**(b) anc is coupled to msb only by XOR** (`cnot(msb,anc)` and the x/cnot/x
uncomputation), so flipping both together restores anc exactly.

**(c) The msb never reaches the observed register.** msb is `b[m−1] = b[n]`, and
the cswaps in `u_a` use `zip(x, b[:n])` — they **exclude** it. So the x register
cannot see the flip.

Together these force `g(y ⊕ w) = g(y)` with w = msb ⊕ anc, hence density ≤ ½.

### Breaking it requires NONLINEARITY in the msb

```
   mode                                 density  defect          w
   v0  baseline                        0.472809       1   b3,anc12
   v6  cswap(t0, msb, anc)             0.472809       1   b3,anc12   <- linear, survives
   v4  toffoli(msb, t0, anc)           0.716064       0   FULL RANK  <- BROKEN
   v5  toffoli(msb, t0, t1)            0.720764       0   FULL RANK  <- BROKEN
```

All four compute `a^e mod N` correctly — the wraps are controlled on scratch
qubits that are 0 on the valid subspace. v6 fails to break the symmetry because
a swap is *linear* over GF(2); only the Toffoli variants, which use the msb
nonlinearly, destroy it.

**Density jumps 0.473 → 0.716, a ~51% cost increase**, from a modification that
changes nothing about what the circuit computes.

### The finding

> **The standard Beauregard/Vedral modular reduction is accidentally friendly to
> Pauli-path simulation.** Its factor-of-two saving is not designed in — it is a
> byproduct of the msb being excluded from the swap network while commuting with
> the adder. A construction that touches the msb nonlinearly forfeits it and
> costs ~51% more.

This is the **first genuine compilation-dependent cost effect** in the project.
Everything else — the Walsh identity, C15, the 2-adic dichotomy — turned out
compilation-invariant, and the one earlier claim to the contrary (F4/F6) was a
bug. Note it is a *constant-factor* effect: 0.716·2ⁿ is still Θ(2ⁿ), so C7's
asymptotic conclusion is untouched.

The framing is also the reverse of the usual one. Not "compile to X to make PPS
cheaper", but "the standard compilation is already cheaper than it needs to be,
for a reason nobody designed."

**Still open:** the broken variants sit at 0.716–0.721, not 1.000, while random
functions reach 1.000. So residual non-linear structure remains after the linear
structure is destroyed. Unidentified.

**→ RESOLVED (TODO 11): see §RS.** It is a conditional linear structure with an
exact ¾ density cap; the quadrant {z_msb=1, z_anc=0} is exactly empty.

## G — THE THEOREM GENERALISES (TODO step 9). V²=id is the whole condition.

> **PARTIALLY SUPERSEDED 2026-08-08 by §WD (TODO 12).** The heading overstates.
> V²=id is the whole condition *among block families of the form tested here* —
> each block controlled on its own fresh qubit. It is **not sufficient in
> general**: `windowed_arith.SelectModExp`'s identity-tail block satisfies
> V²=id exhaustively and still loses the invariance, because it is gated on
> OR(window) rather than on a single qubit. The second hypothesis of the
> theorem below ("controlled on its own qubit") is load-bearing and was
> invisible at step 9 because every block tested there had exactly one control.
> The repaired criterion is in §WD. C29 is regraded accordingly in `CLAIMS.md`.

The C15 proof was stated for the multiply–swap–unmultiply construction. But
re-reading it, steps (i)–(iv) never use V's internals beyond **V² = id** —
"identity on the valid subspace" was context for *why* `a^(2^i)=1` produces such
a block, not a proof ingredient. So the conjecture was that the theorem covers
any construction with an involutive repeated block.

**Confirmed.** Synthetic blocks with nothing to do with modular arithmetic
(`experiment_c15_general2.py`): a controlled swap `x0 ↔ b0` (order 2) and a
controlled 3-cycle `x0 → b0 → b1` (order 3), each controlled on its own fresh
qubit.

```
  N=7 a=6, block = involution        N=7 a=6, block = cycle3
    k  qubits  |support|              k  qubits  |support|
    0      14       3116              0      14       3116
    1      15      15362  SAME        1      15      15362
    2      16      15362  SAME        2      16      30984   +15622
    3      17      15362  SAME        3      17      62088   +46726
    4      18      15362  SAME        4      18     123838  +108476
```

Reproduced at N=5 (constant 15248 vs growth to 123132). Vacuity check passes in
every case: k=0 gives 3116 against k=1's 15362, so the blocks genuinely act on
the observed bit.

### Generalised statement

> **Theorem.** Let a circuit contain k blocks, each the *same* permutation V
> controlled on its own qubit, where V does not modify those control qubits. If
> **V² = id**, the Walsh support of the pullback of any computational-basis
> observable is independent of k (for k ≥ 1).
>
> Proof: steps (i)–(iv) verbatim; they use only V²=id and commutation.

So the scope caveat in Paper B can be **dropped and replaced by a criterion**.
Instead of "this holds for Vedral/Beauregard-style constructions", the claim is
"this holds for any modexp construction whose a=1 block is an involution", which
is a property one can check of a given construction rather than a family one has
to belong to. Multiply–swap–unmultiply satisfies it because `A⁻¹SA` is a
conjugate of a product of disjoint transpositions.

> **Correction (§WD).** The criterion as stated in this paragraph is
> **necessary but not sufficient**. It must also require that the block's
> dependence on its control qubits is *affine* — satisfied trivially by "one
> fresh control qubit" (the identity function) and by "no control at all" (a
> constant function), violated by OR. See §WD.

Worth checking against windowed / table-lookup arithmetic (Gidney-style), where
the a=1 block is not obviously self-inverse — that is now a well-posed question
rather than an open-ended survey.

### Another protocol catch: the control that failed to fail

Pass 1 was **vacuous**. The synthetic blocks acted only on b-register qubits
while the observable was Z_x0, so they never touched the observed bit and the
pullback was unchanged for a trivial reason — both the involution *and* the
3-cycle came out constant. The tell was precisely that **the control did not
fail**. Fixed by making the blocks act on x0 and adding an explicit vacuity
check (k=0 vs k=1 must differ). Without the must-fail control the vacuous test
would have "confirmed" the conjecture for the wrong reason.

## WD — TODO 12: WINDOWED ARITHMETIC. The criterion was incomplete.

`windowed_arith.py`, `test_windowed.py`, `experiments/experiment_windowed.py`.
All predictions declared before the sweeps; two must-fail controls, both failed
as required.

### The question, and why it had a trap in it

Step 9 (§G, C29) reduced Paper B's scope caveat to "the a=1 block is an
involution". TODO 12 asked whether Gidney-style windowed / table-lookup
arithmetic qualifies. **It does — and so does a construction that loses the
invariance completely.** The criterion as stated does not discriminate.

Windowing consumes w exponent bits per block, and there are two natural ways:

- **LOOKUP** (`WindowedModExp`, Gidney-style). Always look up
  `T[j] = a^(j·2^(kw)) mod N` into a scratch register, always multiply x by the
  looked-up *register*, always unlook it up.
- **SELECT** (`SelectModExp`, `skip_zero=True`). Apply the constant multiplier
  `u_a(a^(j·2^(kw)))` controlled on `[window == j]`, **skipping j = 0** because
  multiplying by 1 does nothing. At w = 1 this is exactly `ToffoliModExp`.

Both compute `a^e mod N` with every scratch qubit returned to |0⟩ (checked for
all exponents at N = 5, 7, 9, w = 1, 2, 3; cross-checked against the state
vector and against the unrelated `ToffoliModExp` on the valid subspace).

In the identity tail (`a^(2^(kw)) = 1`, i.e. every window past α when β = 1):

```
                       tail block   V² = id   reads the window?   gated on
  LOOKUP               nontrivial   YES       NO                  (constant)
  SELECT skip_zero=T   nontrivial   YES       YES                 OR(window)
```

Exhaustive over the whole state space at N=5 a=4 (`test_windowed.py` [D]/[E]):
LOOKUP moves 458752 of 524288 states, V²=id, and `V_j` is the *same*
permutation for all four window values including j = 0. SELECT moves 28672 of
65536, V²=id, and `V_0 ≠ V_1 = V_2 = V_3` — i.e. the block is V gated on
OR(window). Control: a **live** window block does read its window, so the
independence result is not vacuous.

### C36 — the repaired criterion

> **V² = id is necessary but NOT sufficient.** C29's other hypothesis — each
> block controlled on **its own qubit** — is load-bearing. The general
> condition is that the identity-tail block's dependence on the exponent
> register be **affine**: the identity function (one fresh control qubit, the
> standard construction) or a constant function (no control at all, the lookup
> construction) both qualify; **OR is nonlinear and does not.**

Step 9 could not see this because every block it tested had exactly one
control, where the two hypotheses coincide. §G is annotated accordingly.

This is the same linear/nonlinear dichotomy that runs through §L2 (a linear
wrap leaves the ½ structure intact, a nonlinear one destroys it) and §AF/§RS,
now appearing one level up — in how a block reads its controls rather than in
how it mixes its data.

### C37 — LOOKUP: the tail is not merely constant, it is DEAD

Because the tail table is all-ones, the QROM permutation is `s ^= 1` regardless
of j, so the tail block never reads its window qubits at all. Predicted before
measuring: the pullback is *constant* in every tail exponent bit, so the Walsh
support is confined to `z_exp` supported on bits `{0 … α−1}` only.

```
  lookup design, live exponent bits in the support   (predicted, then measured)
    N=5 a=4  r=2  α=1   ->  [0]      K=1 |S|=62680    K=2 |S|=28078
    N=7 a=6  r=2  α=1   ->  [0]      K=1 |S|=62570    K=2 |S|=13563
    N=5 a=2  r=4  α=2   ->  [0,1]    K=1 |S|=252111   K=2 |S|=56867
```

**This is a strictly sharper confinement than C24, and the two are directly
distinguishable.** C24 says the standard construction's support *contains* the
all-ones-on-tail vector `1_I` — the tail is alive, through one parity bit. The
lookup design contains **no z with any tail bit set at all**. Measured at the
same modulus and base (N=5, a=4): standard has `1_I` present; lookup has zero
tail bits present. Same logical map, opposite support geometry in the exponent.

### C38 — the invariance survives, but 2-periodic, not exactly constant

The tail applies the fixed involution W **unconditionally**, so the circuit is
`W^K ∘ P_live` and depends on K only through its parity.

```
  |support| vs tail-window count K (w = 2)
    N=5 a=4 (β=1)    28078   62680   28078   62680     2-periodic
    N=7 a=6 (β=1)    13563   62570   13563   62570     2-periodic
    N=7 a=3 (β=3)    84827 1009765 4156585    --       CONTROL: grows
```

Stronger than equal counts: at N=5 a=4 the K=0 and K=2 support **sets are
identical** after dropping the dead tail coordinates — the same function, not
merely the same size. Density falls to zero (0.054 → 0.030 → 0.0033 → 0.0019)
while the β=3 control climbs to the generic ~½ (0.162 → 0.481 → 0.496).

So Paper B's headline survives windowing — cost independent of exponent-register
width when β = 1 — but for this construction it must be stated as **bounded and
2-periodic**, not as the exact constant C24 gives. The two values differ by
2.2× (N=5) and 4.6× (N=7).

### C39 — the cause is the multiply-by-1 branch, not the window

The decisive control is **w = 1**, where there is no window at all, yet the
lookup form still writes `T[0] = T[1] = 1` and multiplies unconditionally.
Predicted before measuring: the tail should still be dead. It is —
N=5 a=4, w=1, K=3: live exponent bits `[0]`, |S| = 62680 (the odd-K value).
So the effect is not "consume several exponent bits at once"; it is
**emitting the multiply-by-1 branch instead of optimising it away.**

### SELECT: skipping j = 0 forfeits the 2-adic advantage entirely

```
  |support| vs K, select design, β = 1
    N=5 a=4    32123  128981  515526  2061438     ×4.02 ×4.00 ×4.00
    N=7 a=6    32098  128956  515388  2060884     ×4.02 ×4.00 ×4.00
```

Exactly **×4.00 = 2^w per window**, which is the derived rate: the tail
dependence runs through `⊕_k OR(window_k)`, and OR of w bits has Walsh sparsity
2^w. In density terms the support is pinned at **0.491–0.492 at every K** —
that is C30's ½ linear-structure cap, i.e. the *generic* β>1 behaviour. So
select-multiply windowing does not merely weaken the β=1 saving, it **removes
it**: a free instance is put onto the same Θ(2^n_exp) curve as a generic one.

Framing worth keeping, and it is the mirror image of §L2/C32: there the
standard reduction was *accidentally cheap*; here the standard practice of
skipping a trivial branch would be *accidentally expensive*. Both are
compilation choices that leave the computed function untouched.

### Putting the j = 0 branch back — the sharpest number in the section

`SelectModExp(skip_zero=False)` emits the multiply-by-1 branch. In the tail all
2^w branches are then the *same* block u_a(·,1) and exactly one fires for any
window value, so the block degenerates to an uncontrolled V. Predicted before
measuring; the flatness half confirmed exactly:

```
  |S| vs K, select design, N=5 a=4, beta=1
    skip_zero=True    32123   128981   515526   2061438      x4.00 per window
    skip_zero=False   32056    32075    32056     32075      2-periodic, flat
    N=7 a=6, False    32024    32072    32024     32072      2-periodic, flat
```

**64× at K=3, and the gap grows as 4^K without bound — from adding a branch
that computes nothing.** Same arithmetic, same qubit count, same everything a
compiler would report.

### An honest split: half of P5 was refuted, and the refutation was useful

P5 was a conjunction — "2-periodic **and** tail bits dead". The first held; the
second **failed**: the live window's other bit (e1) stays in the support, where
α = 1 predicts only e0. Not waved away:

> The **activation ancilla is itself a scratch qubit**, and the pullback ranges
> over all its values. In the LOOKUP design every branch contributes a CNOT
> into s, so the 2^w contributions XOR the ancilla an even number of times and
> it cancels. In SELECT the ancilla *gates* a block, so it does not.

Checked directly on window 0 at N=5 a=4, where the multiplier depends only on
e0: the select block is **not** independent of e1 in general but **is**
independent of e1 on the act = 0 half-space; the lookup block is independent
everywhere. So the refuted half is an ancilla artifact of this particular
activation scheme, not a failure of the mechanism — and it is a reminder that
in this project the pullback sees scratch qubits in every state, which is the
same fact that made the LOOKUP/SELECT distinction exist at all.

### Consistency worth noting: |S| does not depend on w

At matched tail-block parity the lookup design gives **62680 at w = 1, 2 and
3** (21, 21 and 24 qubits). Not a coincidence and not a bug: the arithmetic
skeleton is the same n modular additions per multiply whatever w is, the
lookups reduce to the same net permutation, and the extra qubits w = 3 needs
are untouched — so it is literally the same function in a bigger space.

### VALIDATED AGAINST THE SOURCE (2026-08-08, after the fact)

`WindowedModExp` was built from a reconstruction of Gidney's construction, not
from his papers — the single load-bearing unverified premise in this section.
Now checked against **arXiv:1905.07682 §3.5** (pseudocode `times_equal_exp_mod`)
and **arXiv:1905.09749 §"windowed arithmetic"**. The model holds, and one claim
is stated outright in the source.

**C39's mechanism is Gidney's own design rationale, verbatim:**

> *"We can reduce the number of multiplications that are needed by iterating
> over small windows of the exponent and looking up the corresponding factor to
> multiply by for each one. **This also removes the need for the
> multiplications to be controlled, because the table lookup can evaluate to
> the factor 1 in cases where none of the exponent qubits are set.**"*
> — 1905.07682 §3.5

and 09749: *"the n_e controlled multiplications we needed to perform become
n_e/c_exp **uncontrolled** multiplications"*. So the real construction is the
unconditional-multiply one; `SelectModExp` is precisely what Gidney's design
avoids. **C37–C39 are about the construction people actually propose to run.**

This is §L2/C32's framing a second time, and stronger: Gidney made the
multiply uncontrolled to save Toffolis, and that choice — made for a completely
unrelated reason — is exactly what preserves the β=1 PPS advantage. Nobody
designed it for that.

**The joint-table simplification is confirmed harmless.** Gidney's lookup is
indexed jointly, `table[ei, mi]` with entries `(ke·f·2^j) mod N` over
`ke in kes`. In a tail window every `ke = 1`, so the entries stop depending on
`ei`: the outer index goes degenerate and the lookup does not read the exponent
window — regardless of the multiplication window size. That is the argument
§WD gave from reasoning; it is now checked against the real construction. Our
`m_window = 1` costs nothing here.

**A detail we did NOT model, and it corroborates C38.** Gidney uses a
*relabelling* swap (`a, b = b, a`, free at compile time) where we emit a
physical swap network. It matters for the involution: his block without the
relabel is `(a,b) → (−b, a+b)`, which squares to `(−(a+b), a)` — **not** an
involution. Counting the relabel gives `(a,b) → (a+b, −b)`, which does square
to the identity. So W²=id needs the relabelling included, and our physical-swap
model captures the right permutation. Better: because the relabelling alternates
with block count, his code ends with

> `if a is not target: swap(a, b)`

— i.e. the real compiler emits a physical swap **exactly when the block count
is odd**. C38's 2-periodicity is not an artifact of our modelling; it is
visible in Gidney's own source.

**The unlookup caveat is correct and now precisely sourced.** Gidney uncomputes
a lookup by measuring the output qubits in the X basis and repairing the
resulting phase negations with a smaller fixup lookup (1905.07682, "Uncomputing
a table lookup"; 1905.09749 cites appendix C of [8]), at √L instead of L
Toffolis. That is **not unitary**, so C8's Walsh identity does not apply to it
as written. The logical map is unchanged (deferred measurement), but the object
PPS would propagate is not the same object. Keep the caveat.

### Scope limits to carry into Paper B

- The construction windows the **exponent** and multiplies bit-by-bit over x
  (m_window = 1). Gidney also windows the multiplication — **checked at source
  and it does not affect the argument**, see the validation subsection above.
- The analysis assumes a **unitary** unlookup. Gidney's measurement-based
  uncomputation is genuinely outside it, and this is the one real gap.
- N = 5, 7 at w = 1, 2, 3 and K ≤ 3. The mechanism arguments are exact and
  exhaustive over the state space; the sweeps are two moduli.
  **K extended to 5 — see §OS.**

## OS — TODO 12e: THE ONSET MEASURED AT α = 3 AND 4 (and C7 to 30 qubits)

TODO 12e was not a research question, it was a *compute* question: several
sharp predictions had been made and left unmeasured because the runs were too
slow. With the CUDA backend they are minutes. Files:
`experiments/experiment_c21_onset.py`, `experiment_c7_scale.py`,
`experiment_windowed_scale.py`. All three declare predictions before measuring
and carry must-fail controls.

### OS1 — C21's onset is now MEASURED at α = 3 and α = 4, not merely predicted

This was the single most valuable rerun: C21 ("the support locks at
n_exp = v₂(r) + 1") was exact at α = 1, 2 and, at α = 3 and 4, only
"still growing at the largest width we can reach, as predicted". Now:

```
  N=17 a=2  r=8  α=3      N=17 a=3  r=16 α=4      N=41 a=3  r=8  α=3
 n_exp  q   |support|    n_exp  q   |support|    n_exp  q   |support|
   1   20     255,104      1   20     255,356      1   23   2,070,878
   2   21   1,037,405      2   21   1,037,374      2   24   8,346,567
   3   22   2,093,137      3   22   2,093,202      3   25  16,766,478
   4   23   4,188,525 <    4   23   4,188,537      4   26  33,539,711 <
   5   24   4,188,525 =    5   24   8,379,626 <    5   27  33,539,711 =
                           6   25   8,379,626 =
```

**Locks exactly at n_exp = α + 1 in all three, with strict growth at every
step below it (7/7).** α = 4 is the deepest onset ever measured here, and the
two α = 3 rows are different moduli, so the rule is not a property of N = 17.

**Matched must-fail control, and it is the good kind: same modulus, same α,
same width, only β differs.** N = 41 has elements of order 8 (β = 1) *and* of
order 40 = 5·2³ (β = 5). The β = 5 row grows at every step —
2,072,174 / 8,346,761 / 16,766,484 / 33,539,777 / 67,086,624 — and never
locks. Note how close it is to its β = 1 twin at n_exp = 4 (33,539,777 vs
33,539,711, a gap of 66 in 3.4e7): **at any single width the two are
indistinguishable; only the growth separates them.** That is worth carrying
into the papers — the invariant is not visible in a cost measurement at one
size.

### OS2 — C24 strengthens from equal COUNTS to identical SETS (C43)

C23/C24 (iv) computes the surviving coefficient explicitly, and the expression
does not mention |I| at all — so the proof gives more than the constant-size
statement that was tested. Write each support element as (z_rest, tailflag)
with tailflag ∈ {0,1} saying whether z_I = 0 or 1_I. Then the *set* of such
pairs must be identical at consecutive widths, not merely the same size.

**Derived before measuring, then confirmed 3/3** (N=17 a=2 at n_exp 4 vs 5;
N=17 a=3 at 5 vs 6; N=41 a=3 at 4 vs 5) — bit-for-bit identical sorted key
arrays over supports of 4.2M, 8.4M and 33.5M elements. Logged as **C43**.

The two halves are individually constant as C24 says, and unequal to each
other (e.g. 2,094,285 vs 2,094,240) — which is expected: C24 predicts each
half is width-independent, not that they match.

### OS3 — METHOD: the must-fail control caught a vacuous test of my own

C2 **passed on the first run**, i.e. failed to fail. Diagnosis: the tail
confinement test "z_I ∈ {0, 1_I}" is vacuous at |I| = 1, because a single bit
*is* either all-zeros or all-ones. The β = 5 control was being asked a
question with only one possible answer — and so, silently, were three of P4's
six rows.

Nothing measured was wrong; the scope of what it could testify to was. Both
checks were re-scoped to |I| ≥ 2, which leaves P4 with 3 genuinely
non-vacuous confirmations (still 0 violations) and gives the control 33.5M
violations out of 67M, as required.

**Second time in the project a must-fail control has caught a vacuous
measurement** — the first was step 9 (§G), where synthetic blocks acted only
on b-qubits while the observable was `Z_x0`. Not the third: the project's
other two self-caught errors (a sweep drawing a new random table per t, and
degenerate random tables at small r) came from "vary exactly one parameter"
and from null-model hygiene, not from this rule. The rules are not
interchangeable and the tally should not be inflated.

### OS4 — C7's circuit series extended to 30 qubits

See `experiment_c7_scale.py`. Two series, each varying only N at fixed n_exp,
measured on the exact integer path.

**P0 first, because it gates everything else:** all six logged C7 rows
reproduce bit-for-bit on the exact (`!= 0`) support test — 15,493 / 15,539 /
127,936 / 1,037,322 / 8,347,241 / 8,346,759. So the thresholding artifact of
§GF does **not** touch the circuit-level numbers, which is what the arithmetic
says it should do: coefficients are integers/2^q, and at q ≤ 30 the smallest
nonzero magnitude is 2⁻³⁰ ≈ 9.3e-10, three orders above the 1e-12 threshold.
(It bites at function level, where n is much larger.)

```
  A: n_exp = 2, vary N              B: n_exp = 3, vary N (β>1 only)
   N   n   r    q   |support|  dens    N   n   r    q   |support|  dens
   5   3   4   15      15,493  .4728    7   3   6   16      30,712  .4686
   7   3   6   15      15,539  .4742   21   5   6   22   2,076,089  .4950
  15   4   4   18     127,936  .4880   33   6  10   25  16,766,640  .4997
  21   5   6   21   1,037,322  .4946   35   6  12   25  16,766,067  .4997
  33   6  10   24   8,347,241  .4975   77   7  30   28 134,188,859  .4999
  35   6  12   24   8,346,759  .4975
  77   7  30   27  66,952,166  .4988   slope  A 1.0054   B 1.0083
 143   8  20   30 536,271,623  .4994          pooled 1.0062  (was 1.008)
```

**q = 30, N = 143, a = 5: |support| = 536,271,623, density 0.49944.** Half a
billion Pauli terms, exactly counted, in 851 s on one A4500 — 64× the Hilbert
space C7 settled on. The slope stays at 1.006 bits/qubit over q = 15..30, and
the C30 linear structure w = b_msb ⊕ anc is still exactly a structure (0
support elements with ⟨z,w⟩ = 1) in 13/13 circuits.

Both controls behave: growing a β=1 instance in n_exp instead of N collapses
the density by exactly 4× per two qubits (0.4733 → 0.1183 → 0.0296, ratios
0.2500/0.2500), and an affine observable reads sparsity exactly 1 at every
width. So the ~½ reading is a property of these circuits, not of the
instrument or the register size.

**Open observation, logged not claimed.** Write the *missing fraction* of the
C30 hyperplane, 1 − 2·density. In series A it shrinks with almost perfect
regularity, ~2.1–2.3× per extra bit of modulus:
0.0544 / 0.0239 / 0.0107 / 0.00494 / 0.00234 / 0.00112 at n = 3..8. In series
B it does **not**: 0.0627 / 0.0100 / 0.00063 / 0.00022 at n = 3, 5, 6, 7.
Matched pairs (same N and a, only n_exp differing) show the same split —
N = 21 barely moves (0.01074 → 0.01004) while N = 33, 35, 77 each drop by
8–11×. So there is *something* about how completely the support fills the
hyperplane that depends on n_exp for larger n and not for smaller, and two
series is far too little to name it. **Do not quote the "halves per qubit"
reading; it is series A only.** A future session wanting a cheap question:
sweep n_exp at fixed N ≥ 33 with β > 1 and see whether the deficit really
falls off a cliff between n_exp = 2 and 3.

### OS6 — C38's 2-periodicity holds to K = 5, i.e. three full periods

§WD stopped at K = 3, which is two periods — the bare minimum that can be
called a period at all, and exactly the shape that produced the retracted
"intermediate 2-adic law" (§I: three points of an oscillation read as a trend).
So this was worth two more points rather than none.

```
  K              0        1        2        3        4         5    q at K=5
  N=5 a=4   28,078   62,680   28,078   62,680   28,078    62,680       29
  N=7 a=6   13,563   62,570   13,563   62,570   13,563    62,570       29
  control   84,827 1,009,765 4,156,585 1.66e7 6.66e7  266,296,901      29
  (N=7 a=3, r=6, beta=3 — no tail window, so nothing is dead; ×4 per window)
```

Exact repeats, not near-misses. Two strengthenings beyond the counts, both
predicted first:

- **Sets, not sizes.** The K = 4 support restricted to the K = 0 coordinates
  equals the K = 0 support as a set, and likewise K = 5 against K = 1 — 4/4.
  If W² = id really makes the longer circuit the same permutation on a wider
  register, that is what has to happen, and it does.
- **C37's dead tail survives arbitrary tail length.** Live exponent bits are
  exactly {0…α−1} at K = 4 and 5 in both moduli; no tail bit ever enters the
  support, at any length.

### OS5 — infrastructure that came out of it (all gated in `test_accel.py`)

- **int32 permutation replay.** The replay, not the transform, is the binding
  memory constraint: `idx ^= ((idx >> c) & 1) << t` keeps the array plus two
  temporaries live, so int64 needs ~24 GiB at n = 30 and does not fit on a
  20 GiB card. Images are < 2ⁿ, so int32 is exact for n ≤ 30 and needs ~12
  GiB. **This is what makes q = 30 reachable at all**, and it is why the
  handoff's "n ≤ 30" figure — which was measured on the FWHT alone — was
  optimistic for circuit-level work until now.
- **`accel.pullback_stats`** returns count / density / per-mask GF(2) parity
  counts without moving the support to the host. At q = 30 the support is
  4 GiB as int64; a density sweep wants three scalars.
- **`lab.measure.support(..., exact=True)` and `lab.measure.stats`**, both
  cached under their own keys so exact and thresholded results can be
  compared rather than silently substituted.

## PK — TODO 12d SOLVED: the peak deficit is a COUNTABLE SET, and it is the
## two dominant Fourier modes

`experiments/experiment_c17_deficit.py`. Claim **C44**; C17 regraded from
"the factor is empirical" to derived. Paper A §5 rewritten around it.

**The question.** Permutation-native propagation halves peak memory —
2.000000 exactly for ripple-carry adders, but **1.9997** for modexp, where all
three logged instances satisfied `rot = 2·perm − 2` exactly. A constant
additive deficit of 2, independent of instance size, was unexplained.

### The derivation (read off the gadget, no measurement needed)

Two facts about the standard Clifford+T Toffoli gadget settle it.

1. A Z-type string commutes with every Z-rotation, so it **cannot branch until
   an H turns a Z into an X**. The gadget's only H acts on the target c, so
   only strings carrying Z_c ever leave the diagonal.
2. Once a string is X_c-type, the gadget's four T gates on c are rotations
   about Z_c, and those act *within* span{X_c, Y_c}:
   `X_c → cos·X_c − sin·Y_c`, `Y_c → cos·Y_c + sin·X_c`. **That space is
   closed**, so the four T gates branch ONCE between them, not 2⁴ times.

Hence every Z_c-carrying string contributes exactly 2 Paulis and every other
string exactly 1:

```
    N_max^rot = 2·N_max^perm − |B|,     B = { z ∈ S : z_c = 0 }
```

with S the permutation-native peak set and c the target of the gadget in which
the rotation-level peak falls. **The deficit is not a constant; it is the
number of live strings that miss that gadget's own target qubit.**

### Verified 9/9, and sharper than a count

At the peak, every Pauli has X-support either empty or exactly {c} — the
x-masks are literally `{0, 128}` for the q=14 instances and `{0, 512}` for
q=17 (and `{8}`/`{16}`/`{32}` for the adders, where the identity sector is
absent precisely because |B| = 0). Folding the X_c/Y_c partners back onto
their parent Z-string recovers S **set for set**, with multiplicity 2 on
z_c = 1 and 1 on z_c = 0, in all nine instances.

### Why 2 — and this is the pretty part

```
  instance     perm      rot   deficit   B
  N=5  a=2    6,834   13,666      2      {x0}, {x0,e0}
  N=7  a=3    8,194   16,386      2      {x0}, {x0,e0}
  N=15 a=7   64,070  128,138      2      {x0}, {x0,e0}
  N=5  a=4    8,194   16,386      2      {x0}, {x0,e0}
  N=7  a=6    6,910   13,818      2      {x0}, {x0,e0}
  N=15 a=2   65,538  131,074      2      {x0}, {x0,e0}
  adders          -        -      0      empty
```

B is the same pair in every instance — and **those two strings are exactly the
Walsh coefficients of magnitude ½**, checked against the full final spectrum,
6/6. §W3/W4 found that same pair from a completely unrelated direction
(coefficient-magnitude analysis of the final spectrum: "the two large
coefficients are supported almost entirely on the x-register bit being measured
plus one exponent qubit"). So:

> **The deficit is 2 because the modexp bit function has exactly two dominant
> Fourier modes, and those modes live off the scratch register — they never
> acquire Z_c, so they never double.**

### The correction this forces: the constant belongs to the OBSERVABLE

Must-fail control C2, and it fires hard. Same circuit (N=5, a=2), observable
moved:

```
  observable    perm      rot   deficit
  x0  (standard) 6,834   13,666        2
  b0            12,196   20,388     4004
  anc           12,206   20,398     4014
  t0             1,024    2,048        0     <- ratio exactly 2, like the adders
```

"Modexp has deficit 2" is a statement about the **standard computational-basis
observable**, not about modular exponentiation. Paper A now says so.

### Grading, and what is still open

- `N_max^rot = 2·N_max^perm − |B|` — **derived + verified 9/9.** This is what
  upgrades Paper A §5 from "upper bound of 2, empirical".
- `|B| = 2` for the standard observable, with B the two |c| = ½ modes —
  **verified 6/6, mechanism identified, NOT proved.** Proving it needs a
  characterisation of which peak-time strings avoid the scratch register, which
  is a statement about the peak moment rather than about the final spectrum.
  Do not write it as a theorem.

### Method: the first pass was wrong, and the way it was wrong is worth keeping

Pass 1 (P1/P2) had the identity right but paired the peak against the Z-set at
the **gadget boundary**, reconstructed by re-emitting the circuit and mapping
gate indices to logical ops. 0/9. Diagnosis from dumping the peak set: the peak
sits deep *inside* the gadget, after the Toffoli's own 4-way expansion has run
at rotation level, so no boundary carries the right set — the partner is the
perm-level **peak**. And the gadget target c never needed reconstructing at
all: it is the single bit set in the peak's common x-mask. **A hundred lines of
index bookkeeping were solving a problem that did not exist.** Meanwhile the
checks that only *counted* B (P4/P5) passed 6/6 through both passes, which is
what made the diagnosis quick — the failure was localised to the boundary
machinery, not to the idea.

## BI — BIASED INPUTS (Gangopadhyay–Kumar–Stănică–Gangopadhyay, JAMC 2023)

Paper dropped in the repo root:
`Stability-of-the-Walsh-Hadamard-spectrum-of-cryptographic-Boolean-functions-with-biased-inputs.pdf`
(J. Appl. Math. Comput. 69:3337–3357, 2023). Read; two things in it are
directly ours, one is a warning, and one earlier negative result of ours now
has a principled explanation.

### The dictionary — their biased inputs ARE our biased input STATES

They study `W_f^{(1+ε)}(u) = Σ_x (1+ε)^{wt(x)} (−1)^{f(x)⊕u·x}`, the WHT when
inputs are i.i.d. Bernoulli(½+δ), and define the **stability transform**
`S_f(u) = Σ_x wt(x)(−1)^{f(x)⊕u·x}` as its first-order term in ε.

In our setting the same object appears for a completely different reason. For a
**product input state** `⊗_i (cos θ_i|0⟩ + sin θ_i|1⟩)` we have
`⟨Z^z⟩ = ∏_{i∈z} δ_i` with `δ_i = cos 2θ_i`, so

> **⟨O⟩ = Σ_z c_z ∏_{i∈z} δ_i**

— the PPS expectation is exactly a δ-biased evaluation of the Walsh spectrum,
graded by wt(z). **Verified numerically**: random product state, N=7 a=6
n_exp=3, 16 qubits — direct statevector 0.133418110567 vs the weighted Walsh
sum 0.133418110567, err 8.2e-14. So their whole framework is "PPS with a biased
product input", and their stability spectrum is the sensitivity of a PPS
estimate to input-state bias. Their δ = 0 is our maximally-mixed direction and
δ = 1 is the computational basis.

### C42 — for the REAL Shor input, most of the support is dead weight

Shor's initial state puts the exponent register in |+⟩, i.e. **δ = 0 on those
qubits**, so *every* z with any exponent-register support contributes **exactly
zero** to ⟨O⟩. Combined with C24 (`z_I ∈ {0, 1_I}` on the tail, exp bits below α
free) this predicts a useful fraction of exactly **2^−(α+1)**. Derived, then
measured:

```
  N=7 a=6  α=1   n_exp=2,3,4:  |S|=15549   dead 75.03%   useful 3883
  N=5 a=4  α=1   n_exp=2,3,4:  |S|=15509   dead 74.99%   useful 3879
  N=5 a=2  α=2   n_exp=2:      |S|=15493   dead 75.00%   useful 3873
  N=5 a=2  α=2   n_exp=3,4:    |S|=32143   dead 87.51%   useful 4014
```

75% at α=1 and 87.5% at α=2 = 1 − 2^−(α+1), and the α=2 row **jumps exactly at
n_exp = 3 = α+1**, reproducing C21's onset from an independent direction. The
useful count is itself constant in n_exp, so C15 holds for useful work as well
as for total cost.

**Honest limit — this is not a free speedup.** PPS propagates the observable
backwards and only meets the input state at the end, so a term is only known to
be dead once propagation is finished. It does **not** reduce peak memory (C18)
as stated. What it does is separate *cost* from *useful work*: at α=2 seven of
every eight Pauli terms carried are irrelevant to the answer. Whether the dead
set can be predicted early enough to prune is open and worth asking.

### Reference for the p-biased machinery

O'Donnell, *Analysis of Boolean Functions*, is free on arXiv
([arXiv:2105.10386](https://arxiv.org/abs/2105.10386)) — **Chapter 8** is the
p-biased Fourier expansion, which is the proper setting for the dictionary
above and for anything further in this direction. Use it rather than
re-deriving.

### Why §W's negative result was structural, not bad luck

§W concluded that weight truncation is the wrong knob while δ (magnitude)
truncation is exactly right. Their **Theorem 10** gives the reason: the
S-spectrum is covariant only under **weight-preserving orthogonal** A, and they
exhibit a counterexample showing it is **not invariant under extended affine
equivalence**. Hamming weight is not an affine invariant; the Walsh support and
|c_z| are. So a weight-graded truncation is basis-dependent by construction,
whereas everything else in this project (C8, C40, the PPS cost model) is affine
covariant. §W's failure was forced.

Their grading is on wt(x) and ours on wt(z), but the transform swaps the two
sides and the argument is the same either way. Ties directly to §GF: affine
structure is what the Walsh basis respects, and weight is not part of it.

### What it does NOT give us

No cost model, no permutation circuits, no connection to the 2-adic story.
Their results are about bent and symmetric functions (Theorem 4's bound
|S_f(u)| ≤ (n/2)2^{n/2}, Maiorana–McFarland stability, symmetric
classifications) — none of which our pullbacks are. Cite it for the biased-input
dictionary and for Theorem 10; do not lean on it for anything else.

## GF — THE RECURRING GF(2) PATTERN: one fact, two corollaries, one folklore gap

`experiments/experiment_gf2law.py`. Prior art checked at source first; the
derivation was done before any measurement; 7/7 checks pass with two must-fail
controls failing as required.

### The question

Three results here have the same shape — affine ingredients are free, nonlinear
ones are expensive: §L2/C30–C32 (linear structure ⟹ ½ cap; a linear wrap is
harmless, a Toffoli destroys it), §RS/C34 (one nonlinear monomial demotes it to
a *conditional* structure, ¾ cap), §WD/C36 (affine block control keeps C15,
OR does not). Is that one theorem?

### The answer: one FACT, two corollaries — and one of them was already textbook

The fact is the defining property of the transform: **the Walsh characters of
GF(2)ⁿ are exactly the affine functions**, so `(−1)^affine` is a *single*
coefficient and anything nonlinear is spread. Two distinct corollaries:

**(a) Affine symmetry confines the support.** This is **classical and we must
cite it, not claim it.** Carlet, *Boolean Functions for Cryptography and Coding
Theory*, **Proposition 29**: *"The derivative D_e f equals the null function
(resp. function 1) if and only if supp(W_f) is included in {0ⁿ, e}^⊥ (resp. in
its complement)."* That is C30 (linear structure ⟹ hyperplane ⟹ density ≤ ½)
and C33 (affine structure ⟹ the opposite coset) verbatim. Functions whose Walsh
support is an affine subspace are Carlet's **partially bent** functions.

**(b) Affine gating costs one Walsh coefficient.** If a circuit's dependence on
a register runs through a gating function φ, that register contributes a factor
`1 + |supp(φ̂)|` to the support. φ affine ⟹ φ̂ is a single character ⟹ factor 2
(C24's `z_I ∈ {0, 1_I}`); φ = OR of w bits ⟹ 2^w coefficients ⟹ factor 2^w per
window. §WD measured **exactly ×4.00 at w = 2**. Same fact, different place.

### C40 — the conditional law (the piece with no citation found)

Derived before measuring. Let a coset partition split F₂ⁿ into cells H_u, and
suppose f restricted to H_u has linear structure w_u. Decomposing the transform
over the partition, `c_z = 2⁻ⁿ Σ_u (−1)^{u·z_C} A_u(z')` and Proposition 29
kills `A_u(z')` whenever `w_u·z' = 1`. Hence

> **supp(f̂) ∩ E = ∅ for E = {z : w_u·z = 1 for every u}**, so
> **density ≤ 1 − 2^−d with d = dim span{w_u}** — when that system is
> consistent.

```
  planted, n=14                cells   d  |E|    cap     density  violations
    no conditioning (C30 form)     1   1  8192  0.5000   0.4919      0
    2 cells SHARING one w          2   1  8192  0.5000   0.5000      0
    2 cells, independent w         2   2  4096  0.7500   0.7466      0
    4 cells, independent w         4   3  2048  0.8750   0.8649      0
```

**This corrects §RS.** The conjectured "1 − 2^−(k+1) ladder" is the wrong
parametrisation: the cap is set by the **span dimension of the per-cell
structures**, not by the conditioning depth. Two cells *sharing* a structure
give d = 1 and a ½ cap (measured 0.5000), not ¾. The ladder is only the special
case where each new level contributes one new independent vector.

### C41 — consistency is a PARITY condition, and it can destroy the cap

The system `{w_u·z = 1}` is solvable iff every linear dependency among the w_u
has **even** support. An odd dependency (`w_0 ⊕ w_1 ⊕ w_2 = 0` forces
`1⊕1⊕1 = 1 ≠ 0`) makes E empty and the cap vanishes entirely.

```
  4 cells, ODD dependency    d=2  consistent=False  |E|=0     density 0.9845
  4 cells, EVEN dependency   d=3  consistent=True   |E|=2048  density 0.8671  (cap 0.8750)
```

Same cell count, same construction, caps 1.0 vs 0.875 **purely from the parity
of a dependency**. This was the must-fail control and it failed as required.

### P4 — the project's own numbers are instances

Read off the real v4 modexp (N=5, a=2, wrap `msb_t0_anc`): the t0=0 slice has
structure msb⊕anc, the t0=1 slice has msb. They span d = 2, so
`E = {z : z_msb⊕z_anc = 1 and z_msb = 1} = {z_msb=1, z_anc=0}` — **exactly the
quadrant C34 found empty** (counts 7454 / 8032 / **0** / 7978), cap ¾, measured
density 0.7161. C30 is the m = 0, d = 1 case. So C34's mechanism is not special
to modexp; it is this law with d = 2.

### Scope, and the honest altitude

- **Affine invariance tested, not asserted** (P5): under a random GF(2) change
  of basis the cells become cosets of a generic subspace, the structures
  transform to L⁻¹w_u, the per-cell structures still hold pointwise, E is still
  exactly avoided, and |support| is unchanged. So the law covers arbitrary
  subspace partitions, not just coordinate-aligned ones.
- **This is a COROLLARY of textbook material**, not a deep theorem: Proposition
  29 plus the standard coset decomposition of the Walsh transform. We did not
  find the combination stated and it may be folklore. The contribution is
  recognising that this project's density caps are all instances of it, plus
  the parity condition.
- **PRIOR ART — checked as far as open access allows; residual risk low but
  not zero.** The paywalled item is Carlet–Tarannikov, *Covering sequences of
  Boolean functions and their cryptographic significance*, DCC 25:263–279
  (2002) = ref [326] of Carlet's book. Its **body was not read**. What was
  checked instead, all legitimately:
  - **Carlet's own book** — the comprehensive modern survey *by the same
    author*, which cites [326] on pp. 205, 206, 319 and reproduces its
    definitions and its Walsh characterisation (Def. 47, Prop. 60). §5.5 was
    read in full. A covering sequence is a **single global** sequence λ with
    `Σ_a λ_a D_a f(x)` constant; **partial** covering sequences (ref [231])
    relax that to two levels on a set and its complement. Ours has a
    **different structure vector per cell**, which is neither — no single λ
    reproduces it. The book does not state C40.
  - **The paper's abstract**, from Carlet's own publications page: its
    contributions are characterisations of balancedness / correlation immunity
    / resiliency, subclasses of resilient functions, and degree and
    nonlinearity bounds. None is a per-cell support-confinement law.
  - **No self-archived preprint** exists on the author's page; there is no
    legitimate free copy of the body that was found.
  - **Adjacent and worth citing anyway: Maiorana–McFarland** (book §5.1.1) —
    functions whose restrictions to each coset of a subspace are **affine**,
    whose support is confined by the image of φ. Same *flavour* (structure on
    cells constrains the support), different hypothesis (affine restriction,
    far stronger than a linear structure) and different conclusion.

  Net: C40/C41 may still be folklore, and remain an easy corollary of
  Proposition 29 plus the coset decomposition — that is the altitude to claim
  them at. But "covering sequences already contains this" now looks unlikely.

### Why this is worth having

It converts density measurement into density *prediction*. Given a construction,
find the per-cell structures, take the span, check the dependency parities, and
the cap follows without simulating anything — and it says which modifications
can possibly help: only ones that add an independent structure vector, or that
introduce an odd dependency.

## AF — TODO 11.1: AFFINE structures, and the function-level ½ rows EXPLAINED

`experiment_affstruct.py`; predictions written before measuring, all confirmed.

**Counting lemma (no measurement needed).** A linear structure g(y⊕w)=g(y)
confines the Walsh support to the hyperplane <z,w>=0; an affine structure
g(y⊕w)=g(y)⊕1 confines it to the coset <z,w>=1. Either caps density at ½.
So the broken variants (0.716–0.722 > ½) can carry NO structure of either
kind — TODO 11.1's circuit-level half is settled by counting, before any
measurement. The rank test is blind to the affine kind (a coset off the
origin can still span the full space); the finder translates the support by
one member and takes the GF(2) kernel of the difference set.

**The exact-0.500 rows of §I are exact all-ones structures.** Mechanism,
derived before measuring, confirmed 6/6 with controls:

- complementing every exponent bit is e ↦ 2^t−1−e, which acts on c = e mod r
  as c ↦ (ρ−1−c) mod r, with ρ = 2^t mod r;
- for prime N and even r, a^(r/2) ≡ −1 (the only square root of 1 besides 1),
  and for odd N negation flips bit 0 — the table is antipodal,
  h[c+r/2] = 1⊕h[c];
- when additionally h[−c] = h[c] (i.e. bit_j(x) = bit_j(x⁻¹ mod N) on ⟨a⟩ —
  an instance property), the composition is an exact structure at w = all-ones.

```
 N=7  a=3 r=6 : even t: w = all-ones AFFINE (eps=1), pointwise True; odd t: none
 N=7  a=2 r=3 : even t: w = all-ones LINEAR (eps=0), pointwise True; odd t: none
 N=11 a=2 r=10: nothing, density 1.000  <- negative control: mapping is right at
                t≡0 mod 4 (ρ=6=r/2+1) but bit0 is NOT inversion-symmetric
                mod 11 (2⁻¹=6 keeps it, 4⁻¹=3 flips it)
```

Controls: random f → nothing; planted g = f(y′)⊕y_top → w = top bit AFFINE,
found. Circuit level (all forced by counting, run as tool sanity): v0 →
exactly the known linear w and nothing affine; v4/v5 → empty kernel.

So §I's "measurably non-generic" observation is explained: the 0.500-on-even-t
rows ARE exact structures, with a checkable criterion (right residue class of
t; antipodality; inversion symmetry of the bit table). Note r=3 — odd r, no
antipode — gets its structure from inversion symmetry alone, so the ½ deficit
does not require an even order.

**Grade: PROVED + verified** (the derivation is the proof; the negative
control isolates the inversion-symmetry ingredient). Claim C33.

## RS — TODO 11 RESOLVED: the 0.716 residue is a CONDITIONAL linear structure

`experiment_resid1.py`, `experiment_resid2.py`; predictions in the file
headers, written before measuring.

### The finding

> Injecting a single nonlinear monomial msb∧t0 does not destroy the linear
> structure — it **demotes it to a conditional one**. The support of every
> broken variant avoids the quarter **{z_msb = 1, z_anc = 0}** exactly
> (0 violations, 3/3 instances), capping density at **¾**. The measured
> 0.716 → 0.742 trajectory is the approach to that cap: at n_exp=6 the
> variant sits at 98.9% of ¾ where the baseline sits at 98.7% of its ½ cap.

Physical reading: **a Z on the accumulator's sign bit never appears without a
Z on the comparison ancilla.**

### Mechanism, read from the construction

**(i) Inert half-space.** The t register is restored between reductions
(Cuccaro MAJ/UMA restores the operand; loads are unloaded), so every wrap
sees the INITIAL t. On {y : t0(y)=0} the v4/v5/v45 wraps are literally the
identity — the circuit IS the baseline there — and L2's structure
w = msb⊕anc survives conditionally: g(y⊕w) = g(y) for all y with t0=0.
Proved by construction; verified pointwise, 0 violations (7584 on the
must-fail t0=1 control; the t0=0 slice equals the v0 slice bit-for-bit).

**(ii) Linear firing branch.** On t0=1 the wrap toffoli(msb,t0,·) degenerates
to cnot(msb,·) — a LINEAR gate — and XOR couplings cannot break an XOR
symmetry (L2 pass 1); they *rotate* it. The t0=1 slice carries its own exact
linear structure w′: **v4: msb alone; v5: msb⊕t1⊕anc**. Verified pointwise,
exhaustive, at N=5 and N=7.

**(iii) Walsh consequence.** For <z,w>=1 the t0=0 half cancels pairwise, so
Ĝ(z) is carried by the t0=1 slice alone, and the slice structure kills its
odd-coset spectrum. Chained: support ⊆ {<z,w>=0} ∪ {<z,w>=1 ∧ <z,w′>=0} —
one quadrant exactly empty, density ≤ ¾.

```
 v4 support density by (z_msb, z_anc) quadrant:
              (0,0)    (0,1)    (1,0)    (1,1)
  N=5 a=2    0.9099   0.9805   0.0000   0.9739     0 support elements in (1,0)
  N=7 a=3    0.9159   0.9802   0.0000   0.9822     0
  N=7 a=6    0.9136   0.9771   0.0000   0.9766     0
```

### Why one monomial ≠ two wraps ≠ two monomials

v4, v5 and v45 all inject the **same** monomial msb∧t0 (v45's two wraps share
it): same inert half-space, same cap — all sit at 0.716–0.722. v3x adds a
second independent monomial msb∧t1: the wraps are then inert only on a
QUARTER space, a quarter-space invariance forces no exact zeros at all, and
density jumps to 0.980–0.983. The residue was never "intrinsic to the
arithmetic" — it tracks the number of independent nonlinear monomials
coupling the msb, not the number of wraps.

### Independence from C15 (derived from C29, then verified)

The wrapped a=1 block is A′⁻¹SA′ — still a conjugate of the involutive cswap
layer — so C29 applies verbatim and constancy survives the breaking: v4 at
N=7 a=6 has support exactly **23488 at n_exp = 2, 3, 4**. The ~51% penalty
does not forfeit the free exponent register. Claim C35.

### The ladder, and corrections to C32's framing

> **CORRECTED 2026-08-08 by §GF (C40/C41).** The ladder below is the wrong
> parametrisation. The cap is **1 − 2^−d with d = dim span{per-cell
> structures}**, plus a **parity condition** on dependencies among them — not
> a function of the conditioning depth k. Two cells *sharing* a structure give
> d = 1 and a ½ cap, not ¾; and an odd dependency removes the cap entirely at
> unchanged d. The k = 0, 1 rows below are right for the reason §GF gives
> (d = 1 and d = 2), and the "two independent monomials → no forced cap" row is
> most likely the parity case, not a larger d. Also: the linear→hyperplane step
> is Carlet's **Proposition 29**, classical — cite it.

baseline: cap ½ (proved, L2) → one monomial: cap ¾ (this section) → two
independent monomials: no forced cap (0.98 measured). The pattern
1 − 2^−(k+1) for codim-k inert subspaces is a natural conjecture but is NOT
claimed; only k = 0, 1 are established.

C32's verb was wrong: a single nonlinear coupling **halves the forced-zero
set** rather than removing it. The "~51%" is a finite-size snapshot of a
¾-cap function beside a ½-cap function (asymptotic ratio 3/2 exactly). And
"0.716" is not a constant: it is instance-dependent (23464 / 23488 / 23579
per 32768) and width-dependent (rising toward ¾). TODO 11.5 closed.

### Open residue, honestly

The three live quadrants sit at 0.91–0.98 — below the generic ~0.995 zero
rate — with the (0,0) quadrant lowest. A further, smaller deficit exists,
presumably deeper conditional levels inherited from the baseline slice. Not
pursued.

**Grades.** Conditional survival on the inert half: **PROVED + verified**.
Quadrant law and ¾ cap: **verified 3/3 instances; mechanism proved except one
link** — the existence of the rotated slice structure w′ is exhaustively
verified per instance and argued from wrap-linearity, not yet proved in
general. Claim C34.

## I — IS THERE AN INTERMEDIATE 2-ADIC LAW? (TODO step 6). No. Answer is negative.

N=323 (r = 144 = 16·9) had shown density 0.981 at t=16 rising to 1.000 at t=24,
which looked like it might be a quantitative law in α rather than C15's binary
split. It is not.

### Two methodology bugs of my own, both instructive

1. **Degenerate tables.** Testing `g(e) = h[e mod r]` with random h: for r=3
   there are only 8 possible tables and 2 are constant, so small-r rows were
   contaminated and appeared *sparse*. Fixed by rejecting constants.
2. **Varied two things at once.** I drew a *new* random table for each t, so
   t-dependence was confounded with table variance — a direct violation of the
   "vary exactly one parameter" rule in `METHOD.md`. Fixed by fixing h and
   sweeping t.

Both produced confident-looking numbers before being caught. Worth remembering
that the protocol's rules catch *my* errors, not just other people's.

### What is solid

**β = 1: sparsity exactly constant in t.** Consistent with the C15 proof.

```
  r= 4 alpha=2 beta=1: sparsity = 4 4 4 4 4 4 4 4 4 4 4 4   (t = 11..22)
  r= 8 alpha=3 beta=1: sparsity = 4 4 4 4 4 4 4 4 4 4 4 4
  r=16 alpha=4 beta=1: sparsity = 16 16 16 16 16 16 ...
```

**β > 1: sparsity is Θ(2^t)** — density bounded well away from 0 at every t.

### What is NOT a law

Density within β>1 does **not** organise by α. With a fixed table and
consecutive t it shows a strong period-`ord₂(β)` oscillation *plus* a slow
upward drift, and neither component is a function of α:

```
  r=12 (beta=3, ord2=2): 1.000 0.625 1.000 0.625 ...   clean period 2
  r= 5 (beta=5, ord2=4): t≡1 mod 4 gives 0.791, 0.815, 0.832 ...  DRIFTS
  r= 7 (beta=7, ord2=3): t≡2 mod 3 gives 0.834, 0.831, 0.850, 0.875 ... DRIFTS
```

So the periodicity hypothesis (H1) is **refuted**: residue classes are not
constant, they drift. Density is periodic-plus-drifting, and function-dependent.

**Conclusion: there is no clean intermediate law. The robust structure is the
binary β=1 / β>1 dichotomy that C15 already asserts.** The N=323 observation is
best explained as sampling aliasing — t = 16, 20, 24 hit residues 4, 2, 0 mod
ord₂(9)=6, so three different points of the oscillation were read as a trend.

This is a negative result that *protects* the paper: it rules out a complication
rather than adding one, and it means C15 should be stated as the binary split
without hedging about intermediate regimes.

### Scope caveat CLOSED — rechecked with the real modexp bit function

`experiment_c15_realfn.py` reruns the same fixed-function, consecutive-t design
using the genuine table `h[c] = bit_j(a^c mod N)` instead of a random one.

**P1 holds.** β=1 gives constant sparsity across t=11..22 (4 for r=4; 8 for
r=16, so ≤ 2^α without being tight).

**Step 6's conclusion is confirmed: α is not the controlling parameter.**

```
  r= 3 alpha=0 beta= 3:  1.000 0.500 1.000 0.500 ...   exact period 2
  r= 6 alpha=1 beta= 3:  1.000 0.500 1.000 0.500 ...   identical to alpha=0
  r=12 alpha=2 beta= 3:  all 1.000
  r=10 alpha=1 beta= 5:  0.770 1.000 0.791 1.000 ...   drift 0.054
  r=18 alpha=1 beta= 9:  ~0.93-1.000                    drift 0.053
  r=22 alpha=1 beta=11:  ~0.92-1.000                    drift 0.015
```

α=1 produces wildly different behaviour across r = 6, 10, 18, 22, while α=0
(r=3) matches α=1 (r=6) exactly. **α does not organise the data.** The period is
set by ord₂(β), the values by β. Same qualitative structure as the random-table
study, so that study's conclusion transfers.

### NEW: the real modexp table is measurably NON-generic

Same r, real table vs random table:

```
  r= 6 (beta=3):  real  1.000 0.500 1.000 0.500 ...
                  random 1.000 1.000 1.000 1.000 ...
  r=10 (beta=5):  real  0.770 1.000 0.791 1.000 ...
                  random 0.500 0.500 0.394 0.197 ...
```

The real function is *sparser* than random at r=6 (exactly 0.500 on even t — a
clean factor of two, suggesting a dead variable or parity constraint) and
*denser* at r=10. So modexp bit functions carry structure beyond "depends on
e mod r". Not pursued further; noted as an open observation, and it means
random-table results should not be used as a proxy for real densities — only
for the qualitative α question they were built to answer.

**→ Now explained: see §AF.** The exact-0.500 rows are all-ones linear/affine
structures with a checkable three-ingredient criterion.

## T — DOES C15 SURVIVE TRUNCATION? (TODO step 7). Peak cost: yes. Accuracy: to a point.

Every C15 figure was δ=0, and nobody runs PPS at δ=0 — so the practical claim
("period-finding precision is free when r is a power of two") was unsupported in
the regime where it would actually be used.

**Predictions were written before measuring** (`experiment_c15_trunc.py` header):

- **P1** Terminal truncation *must* preserve constancy at every δ, since C24/M2
  give identical magnitude multisets across n_exp and a magnitude threshold on
  identical multisets keeps identical counts. Zero exceptions expected; failure
  would mean C24 or M2 is wrong.
- **P2** Incremental truncation (what PPS does) *may* break it — more n_exp means
  more gates, hence more truncation events, and W1 already showed incremental and
  terminal are different operations.
- **P3** Control: β>1 must fail at every δ.

### T1 — P1 confirmed exactly

```
 N=7 a=6 (r=2, beta=1)          delta:  0    1e-6   1e-4   1e-3   1e-2  3e-2  1e-1
   n_exp=3                            15549  15549  15549  13194  1712    75     4
   n_exp=4                            15549  15549  15549  13194  1712    75     4
   n_exp=5                            15549  15549  15549  13194  1712    75     4
 N=5 a=2 (r=4, beta=1)
   n_exp=3,4,5                        32143  32143  32143  26257  1704    36     8   (identical)
 N=7 a=3 (r=6, beta=3)  CONTROL
   n_exp=3                            30712  30712  30712  23849  1442    40     4
   n_exp=4                            64353  64353  62982  44480  1225    51     4
   n_exp=5                           129012 129012 123137  69849   994    43     4
```

Identical across every width at every δ for β=1; the control diverges. **The
invariance is not an artifact of working exactly — it holds at every truncation
level**, because the underlying magnitudes are identical.

### T2 — P2: the practically relevant cost survives; accuracy has a limit

Incremental (real PPS), N=7 a=6:

```
  n_exp   delta     N_max   N_final         <O>
      3   0e+00     24369     15549   -1.000000
      3   1e-02      1028         8   -1.000000
      3   1e-01        34         8   -1.000000
      4   1e-02      1028        16   -1.000000
      4   1e-01        34         0   +0.000000     <- accuracy lost
      5   1e-02      1028        32   -1.000000
      5   1e-01        34         0   +0.000000

  constancy across n_exp:  delta=0     N_max SAME  N_final SAME  <O> SAME
                           delta=1e-4  N_max SAME  N_final SAME  <O> SAME
                           delta=1e-2  N_max SAME  N_final differs  <O> SAME
                           delta=1e-1  N_max SAME  N_final differs  <O> differs
```

**N_max — peak memory, the quantity that actually bounds cost — is SAME at every
δ tested** (24369, 24369, 1028, 34). That is the practically relevant statement
and it survives intact.

What does *not* survive: `N_final` drifts at aggressive δ (8 → 16 → 32), and at
δ=1e-1 the estimate collapses to ⟨O⟩ = 0 for the larger widths. Mechanism is
W1's: more exponent qubits means more gates, hence more incremental truncation
events, so terms that a terminal threshold would have kept are destroyed en
route. The larger circuit is *more* fragile at the same δ despite having an
identical exact spectrum.

### T3 — control confirmed under incremental truncation

The β>1 control was cut short by the CPython 3.14 crash (trap 7) but got far
enough to settle the question:

```
  N=7 a=3 (r=6, beta=3) CONTROL
   n_exp   delta     N_max   N_final         <O>
       3   0e+00     48855     30712   -1.000000
       4   0e+00     98018     64353   -1.000000
       5   0e+00    196060    129012   -1.000000     <- doubling per step
       3   1e-04     48036     30779   -1.001772     <- norm violation
       4   1e-04     94827     62819   -0.997549
```

N_max grows ~2× per added qubit, against SAME at every δ for β=1. And δ=1e-4
produces ⟨O⟩ = −1.001772, i.e. |⟨O⟩| > 1 — another instance of C5, now on the
control. N=5 a=2 (β=1) likewise gave N_max 48972 at both n_exp=3 and 4.

### Verdict for Paper B

The practical claim stands, with a stated boundary:

> Peak Pauli-path cost is independent of exponent-register width for β=1, at
> every truncation level tested. Accuracy is likewise independent up to
> moderate δ, but at aggressive δ the wider circuit degrades first — not because
> its exact spectrum differs (it is identical) but because incremental
> truncation has more gates to act on.

So "period-finding precision is free" is true for *memory* and true for
*accuracy up to moderate δ*, and must not be stated unqualified.

## X — CRYPTANALYSIS IMPORT (TODO step 5). A usable bound, tight at the extremes.

C12 gave a qualitative bridge. This makes it an inequality and tests it against
**published** nonlinearity values rather than self-measured ones.

### X1 — the bound

With normalised coefficients (Σc²=1 by Parseval) and NL = 2ⁿ⁻¹(1 − max|c|):

```
    S · max|c|² ≥ Σ c_z² = 1      ⟹      S ≥ (1 − NL/2ⁿ⁻¹)⁻²
```

Since S is exactly the Pauli support (Paper A), **any published nonlinearity
gives a lower bound on PPS cost, for every circuit computing that function,
compilation-independent, without running anything.**

```
                  function   n  S (actual)      bound  holds  tight       NL  NL/NLbent
         affine (x0^x1^x3)  10           1        1.0   True   True        0     0.0000
  inner product (bent) m=4   8         256      256.0   True   True      120     1.0000
  inner product (bent) m=6  12        4096     4096.0   True   True     2016     1.0000
           AES S-box bit 0   8         239       64.0   True  False      112     0.9333
       AES min over combos   8           -          -      -      -      112   (published: 112)
         modexp N=5 (Z_x0)  14        3086        4.0   True  False     4096     0.5039
                  random f  12        4000      264.2   True  False     1922     0.9534
```

**Independent validation:** the AES S-box nonlinearity comes out at exactly
**112**, the published value, over all 255 nonzero linear combinations. That is
a check of the whole Walsh pipeline against an external constant.

### X2 — closing the loop end-to-end

`Toffoli(x_i, y_i, out)` over i computes `out ^= <x,y>`, the inner product,
which is **bent**. Bent ⟹ flat spectrum ⟹ PPS must carry every Z-string:

```
    m  qubits  perm_pps |supp|  predicted   match
    3       7               64         64    True
    4       9              256        256    True
    5      11             1024       1024    True
```

Support is exactly 2^(2m) — the output qubit forced into every term, the bent
structure filling the rest. **No truncation can help**: every coefficient has
identical magnitude, so there is no "large" subset to keep. This is the clean
worst-case statement, and it is derived from a published property.

### X3 — HONEST LIMITATION: the bound is weak in the middle

Tight only at the two extremes, where Parseval is saturated by a single
magnitude:
- affine: bound 1, actual 1 ✓
- bent: bound 2ⁿ, actual 2ⁿ ✓
- AES: bound 64, actual **239** — off by 3.7×
- modexp: bound 4, actual **3086** — off by 770×

So as a *quantitative* predictor it is poor except at the endpoints, because NL
depends only on `max|c|` and discards the rest of the spectrum. The valuable
statements are the qualitative direction and the exact bent worst case, not the
numeric bound in between.

Sharper bounds would need more than NL — e.g. the full Walsh value/multiplicity
distribution, which *is* published for several crypto families (the AES inverse
has a known spectrum). Worth a follow-up: for functions with a fully published
spectrum, S is known **exactly**, not bounded.

### X4 — correction to a C12 claim

Earlier: "modexp sits at a stable 0.74 of the bent bound across 15→21 qubits."
That stability is across **N at fixed n_exp=2**, not across widths. At n_exp=1
the same instance gives 0.5039, because max|c| is 0.5 there versus 0.268 at
n_exp=2. **Do not state it as width-independent.**

## W — WEIGHT TRUNCATION (TODO step 2). Naive identity fails; knob is unusable.

PPS has two truncation knobs: coefficient threshold δ (all prior work here) and
Pauli **weight** (drop terms acting on > k qubits). This closes the gap.

Setup: for a permutation circuit the surviving Paulis are Z-strings, and the
Pauli weight of `Z^z` is `popcount(z)` — which *is* the Fourier **degree** of the
Walsh coefficient. So weight-k truncation is literally low-degree Fourier
truncation. The hope was that its error would be the Fourier tail, extending the
exact model to the second knob.

### W1 — the naive identity is FALSE (and the reason generalises)

**PPS truncates incrementally, at every gate.** A term discarded early never
branches, so the final result is not the truncated final operator. The Fourier
tail describes **terminal** truncation (truncate the finished operator) only.

```
modexp N=5, 14 qubits, exact <O> = -1.000000
  k   terminal <O>   term err   incremental <O>   incr err   tail mass
  2      -1.015625   1.56e-02         -0.500000   5.00e-01    0.499298
  3      -0.898438   1.02e-01         -1.000230   2.30e-04    0.494888
  4      -0.820312   1.80e-01         -1.000648   6.48e-04    0.470779
  8      -1.242188   2.42e-01         -0.644562   3.55e-01    0.188095
 14      -1.000000   0.00e+00         -1.000000   0.00e+00    0.000000
```

Both are **non-monotonic in k**, and incremental lands on either side of
terminal — at k=3 it is 400x *more* accurate, at k=8 far worse. This same
incremental-vs-terminal distinction applies to δ, and retroactively explains why
δ behaves non-monotonically too (C14): it is not just cancellation, it is that
truncation changes what subsequently branches.

### W2 — weight truncation is the WRONG KNOB for reversible arithmetic

Signed sum of coefficients per weight level (modexp N=5):

```
 weight   #terms      mass   signed sum
      1        3  0.250031    -0.500000
      2       12  0.250671    -0.515625
      3       54  0.004410    +0.117188
      4      163  0.024109    +0.078125
      5      316  0.052185    -0.140625
      7      601  0.071869    +0.312500
      8      574  0.084717    -0.515625
      9      451  0.089630    +0.296875
     13        8  0.002197    -0.031250
```

The levels carry **large, alternating** signed sums that cancel only when summed
over *all* of them. Any weight cutoff slices through that cancellation. There is
no k that both reduces cost and preserves ⟨O⟩.

### W3 — coefficient truncation is the right knob, now *exactly* explained

For the same circuit, the `|c| > 0.1` set (**4 terms** of 3086):

```
  they sum to             : -1.00000000    (exact <O> = -1.000000)
  everything else sums to : +0.00e+00      (3082 terms)
```

An **exact** split, not an approximation. δ-truncation works because the large
coefficients form a set that carries ⟨O⟩ exactly while the remainder sums to
exactly zero. Weight truncation fails because the same split is not aligned with
degree.

### W4 — the dominant coefficients are NOT all low-degree

Weights of the four: **1, 2, 8, 9**. Two of them sit high. So low-degree
truncation cannot substitute for δ — it would discard half the terms that carry
the answer. This is the concrete reason W2 holds.

Corrects an earlier misreading: the apparent "292x low-weight enrichment"
(mass 0.50 at weight ≤ 2 vs 0.026 for a random function) is **just those two
0.5-magnitude coefficients** at weights 1 and 2 (0.5² + 0.5² = 0.5). It is not a
broad low-degree structure — `k(99% mass) = 11` confirms the bulk is high-weight.

### W5 — the (z, z⊕e) pairing lead: REFUTED

Top-12 listing showed every coefficient appearing twice at equal magnitude,
differing by an exponent qubit. That looked like a hard group-theoretic
constraint and the best route to proving C15 (step 3). It is not:

```
   a   r  n_exp  |support|   all exp bits paired?
   6   2      1       3116                   True
   6   2      2      15549                  False
   3   6      1       3206                   True
   3   6      2      15539                  False
```

Exact at n_exp=1 for **both** r=2 and r=6, broken at n_exp≥2 for both
(≈0.4% of partners missing). So it is an n_exp=1 artifact, not r-dependent, and
**not the C15 mechanism. Step 3 remains open.**

Worth keeping: at n_exp=1 the pairs are almost all *opposite*-sign (3084
opposite vs 2 same), so nearly everything cancels pairwise and the two
same-sign pairs carry ⟨O⟩ — consistent with F13/W3.

**Bottom line for step 2:** the exact model does *not* extend to weight
truncation, and weight truncation should not be used on reversible arithmetic.
That is a negative result, but a sharp and actionable one, and W1 gives a
mechanism that also improves the understanding of the δ results.

## C7 RESOLVED — the results are NOT pre-asymptotic

F9 removed the need for a faster simulator: since PPS term count *equals* Walsh
sparsity, exact PPS cost is computable in O(2ⁿn) without running PPS. That
reached **24 qubits** (dim 1.7e7) where PPS itself stalls around 15–17.

`experiment_c7.py` / `/tmp/c7a.log`, Toffoli modexp, observable Z_x0:

```
   N   n   r  qubits        2^q   sparsity  density  r pow2    time
   5   3   4      15      32768      15493  0.47281     yes    0.1s
   7   3   6      15      32768      15539  0.47421      no    0.0s
  15   4   4      18     262144     127936  0.48804     yes    0.8s
  21   5   6      21    2097152    1037322  0.49463      no   16.5s
  33   6  10      24   16777216    8347241  0.49753      no  258.2s
  35   6  12      24   16777216    8346759  0.49751      no  262.1s

  log2(sparsity) grows 1.008 bits/qubit  (1.000 = exactly Theta(2^n))
  density: 0.473 -> 0.474 -> 0.488 -> 0.495 -> 0.498 -> 0.498
```

**Density converges monotonically to 1/2 and the growth slope is 1.008
bits/qubit.** So PPS on modular exponentiation costs Θ(2ⁿ) — asymptotically no
better than a state-vector simulation. The small-n numbers were already in the
asymptotic regime; C7 is answered, and negatively for the threat.

For contrast, PPS at 17 qubits ran >20 min without finishing while the Walsh
computation of the same quantity took **0.01s**.

## C15 CONFIRMED AT CIRCUIT LEVEL — and it is the sharpest result so far

The tension in F12 (function-level dichotomy vs circuit-level density stuck at
~0.5) is resolved: **the earlier runs all used n_exp=2, where a 2-bit exponent
register leaves no room for periodicity to show.** Sweep n_exp and the dichotomy
appears immediately.

Controlled design: fix N (identical circuit, ancillas, gate structure), vary `a`
so that only r changes. `experiment_c15.py`, `experiment_c15b.py`.

**N=7, a=6 (r=2, a power of two), observable Z_x0:**

```
 n_exp  qubits        2^q   sparsity   density  vs n_exp=2    time
     2      15      32768      15549  0.474518        SAME    0.0s
     4      17     131072      15549  0.118629        SAME    0.4s
     6      19     524288      15549  0.029657        SAME    2.9s
     8      21    2097152      15549  0.007414        SAME   27.3s
    10      23    8388608      15549  0.001854        SAME  176.2s
```

**CONTROL — N=7, a=3 (r=6, has an odd factor), same circuit size:**

```
 n_exp  qubits   sparsity   density    growth
     2      15      15539  0.474213         -
     4      17      64353  0.490974     4.14x
     6      19     258691  0.493414     4.02x
     8      21    1037728  0.494827     4.01x
    10      23    4155634  0.495390     4.00x
```

Sparsity is **exactly** 15549 across a 256x growth in Hilbert-space dimension
for r=2, against a clean 4.00x-per-step growth for r=6. Same modulus, same
circuit, only the base differs.

**Independent confirmation at N=21** (`experiment_c15.py` §2): a=8 (r=2) gives
1037174 at both n_exp=2 and n_exp=4 while the dimension quadruples; a=2 (r=6)
goes 1037322 → 4186980 and a=4 (r=3) goes 512784 → 4152181.

**Mechanism.** aᵉ mod N depends only on e mod r. If r | 2ᵏ then only the low k
bits of e matter, so every additional exponent qubit adds a variable the Walsh
support cannot touch — the support is pinned to a fixed subspace. Any odd factor
in r makes the period incommensurate with the GF(2) basis and the support
spreads over everything.

**Statement of the result.**

> For reversible modular exponentiation with computational-basis observables,
> Pauli-path simulation cost is **independent of the exponent-register size when
> the order r is a power of two**, and **Θ(2^q) as soon as r has an odd factor.**

Consequences:

1. **Period-finding precision is free, or fatal, depending on r.** The exponent
   register is what sets the accuracy of the continued-fractions step; here
   enlarging it costs PPS *nothing* when r is a power of two and quadruples cost
   per two qubits otherwise.
2. **The textbook demo is the degenerate case, quantitatively.** N=15 a=7 has
   r=4. Every "we simulated Shor on N=15" result sits in the corner where PPS
   cost does not grow at all. Cryptographic N has r with odd factors
   generically, i.e. the Θ(2ⁿ) branch.
3. **Same 2-adic dependence as MPS.** Dang, Hill & Hollenberg report memory
   depending on *the factors of r*
   ([arXiv:1712.07311](https://arxiv.org/abs/1712.07311)). Two structurally
   unrelated classical methods keying on the same arithmetic property is worth
   stating as a shared fact about the algorithm, not a coincidence of either.

**Caveat.** `live_variables` reports all q bits live even in the r=2 case, so the
constancy is *not* simply "the function ignores the extra qubits". The extra
exponent qubits do influence the function (through behaviour on invalid inputs)
without enlarging the Walsh support. The clean subspace argument above explains
the valid-input structure; the full-space statement is empirical over 5 sizes
and 2 moduli, not proved.

## F12 — the 2-adic dichotomy at function level (superseded by C15 above)

Function-level probe (`/tmp/c7b.log`), g(e) = bit 0 of aᵉ mod N over e ∈ [0,2ᵗ),
which isolates the algorithm from ancilla layout:

```
    N   a    r   t        2^t   sparsity  density  r pow2
   15   7    4  12       4096          4  0.00098     yes
   15   7    4  24   16777216          4  0.00000     yes    <- CONSTANT in t
   21   2    6  24   16777216   16777216  1.00000      no    <- FULLY DENSE
   35   3   12  24   16777216   16777216  1.00000      no
  143   5   20  24   16777216   16777216  1.00000      no
  323   5  144  16      65536      64293  0.98103      no
  323   5  144  24   16777216   16777216  1.00000      no
```

Sharp dichotomy, and the reason is elementary: if r | 2ᵏ then aᵉ mod N depends
only on the low k bits of e, so g is a function of k variables and its Walsh
support is ≤ 2ᵏ — **constant in t**. If r has any odd factor, the periodicity is
incommensurate with the GF(2)ᵗ Walsh basis and the spectrum is **maximally
spread**. N=323 (r=144=16·9) shows the intermediate case: the factor 16 buys
partial sparsity at small t, washed out by t=24.

Two consequences:

1. **This is the same phenomenon Dang, Hill & Hollenberg report for MPS** — that
   memory depends on *the factors of r* rather than r itself
   ([arXiv:1712.07311](https://arxiv.org/abs/1712.07311)). Two unrelated
   classical methods, same 2-adic dependence. Worth stating as a shared
   structural fact rather than a coincidence of either method.
2. **The textbook demo is the degenerate case.** N=15, a=7 has r=4, a power of
   two — sparsity 4, constant forever. Every "we simulated Shor" result on N=15
   sits in the trivially-simulable corner. For cryptographic N, r is generically
   not a power of two, so PPS is at the fully-dense worst case.

## F13 — the spectrum shape explains F11 (aggressive truncation wins)

F11 was empirical; C12's spectrum data explains it. modexp N=5, Z_x0:

```
   keep |c| >    #kept     sum kept      |err|
        0e+00    15493    -1.000000   0.00e+00
        1e-03    13005    -0.987061   1.29e-02
        1e-02     1579    -0.955322   4.47e-02
        3e-02       75    -0.771973   2.28e-01   <- worst
        1e-01        4    -1.000000   0.00e+00   <- exact, 4 terms
```

**4 of 15493 coefficients reproduce ⟨O⟩ exactly; the other 15489 sum to zero.**
The spectrum is heavy-tailed by a factor of 33 (max|c| = 0.2676 against the
flat-spectrum value 1/√S = 0.0080). Keep the spikes → exact. Keep the spikes and
*some* of the cancelling sea → residue. That is the non-monotonicity, now
predicted from spectrum shape rather than observed.

Practical rule, counterintuitive and actionable: **for these circuits truncate
aggressively, not mildly.** Mild δ is the worst regime.

## C12 VALIDATED — the cryptanalysis bridge is quantitative

With c_z normalised, L = 2ⁿ·max|c|, NL = 2ⁿ⁻¹(1 − max|c|), bent bound
2ⁿ⁻¹ − 2^(n/2−1). Parseval gives Σc² = 1, so a flat spectrum over S terms has
each |c| = 1/√S and NL ≈ 2ⁿ⁻¹(1 − 1/√S) — tying PPS cost S to nonlinearity.

```
                  function   n  sparsity   dens   max|c| 1/sqrt(S)  NL/NLbent
         adder b0 (affine)  12         1  0.000  1.00000   1.00000    0.00000
                  adder b1  12         4  0.001  0.50000   0.50000    0.50794
                  adder b4  12        46  0.011  0.50000   0.14744    0.50794
            modexp N=5 a=2  15     15493  0.473  0.26758   0.00803    0.73649
           modexp N=15 a=7  18    127936  0.488  0.25763   0.00280    0.74382
           modexp N=21 a=2  21   1037322  0.495  0.25588   0.00098    0.74463
            random f, n=14  14     16384  1.000  0.03357   0.00781    0.97404
            random f, n=18  18    262144  1.000  0.00965   0.00195    0.99229
```

Both endpoints land exactly: affine ⟹ sparsity 1, NL 0. Random ⟹ density 1.000,
NL/NLbent 0.97–0.99 (near-bent, as expected). **modexp sits at a stable ~0.74 of
the bent bound** — strongly nonlinear, but measurably *not* random, and the
figure is flat across 15→21 qubits.

**Caveat that must be stated:** the Parseval prediction over-estimates NL by
~35% for modexp (measured 12000 vs predicted 16252) because the spectrum is
heavy-tailed, while for random functions it is accurate to 1–3%. So S ↔ NL is a
**bound that is tight only for flat spectra**, exact at the affine and bent
endpoints, not an identity. Do not overclaim it.

## F10 — scope of the Walsh identity (exactly where it stops)

`experiment_scope.py` Q1, Toffoli modexp N=5, 14 qubits:

```
         observable   type   walsh     pps   nonZ  match
               Z_x0      Z    3086    3086      0    YES
               Z_x1      Z    2926    2926      0    YES
          Z_x0 Z_x1      Z    2848    2848      0    YES
     Z_x0 Z_x1 Z_x2      Z    2514    2514      0    YES
               X_x0      X       -   60358  60358    n/a (hit cap)
               Y_x0      Y       -   60358  60358    n/a (hit cap)
```

Holds for **every Z-type observable**, single- or multi-qubit (generalise via
`(-1)^{popcount(π(y) & zmask)}`). Fails completely for X/Y-type: the pullback
leaves the diagonal, every surviving term is non-Z, and it blows past the cap.

Honest scope statement: **computational-basis observables on permutation
circuits.** That is exactly what one measures in Shor (bits of the output
register), so the restriction is natural rather than convenient — but it must be
stated, not glossed.

## F11 — truncation error is NON-MONOTONIC in δ (aggressive beats mild)

`experiment_scope.py` Q2, Toffoli modexp N=5, truth = −1:

```
   delta    N_max        <O>     |err|  admissible
   0e+00    40770   -1.00000  6.66e-16         yes
   1e-04    39992   -1.00348  3.48e-03    VIOLATED
   1e-03    23482   -0.71532  2.85e-01         yes
   1e-02     2052   -1.00000  0.00e+00         yes
   3e-02      514   -1.00000  0.00e+00         yes
   1e-01       34   -1.00000  0.00e+00         yes
```

**δ=1e-1 gives the exact answer with 34 terms; δ=1e-3 is off by 0.285 with
23482 terms.** Three orders of magnitude more work, far worse answer.

Mechanism: the coefficient spectrum is a handful of large terms plus ~10⁴ small
ones (each ≈2⁻ⁿ) that **cancel among themselves**. Discard all of them and the
cancellation is preserved exactly. Discard *some* and the residue survives.
Mild truncation is the worst regime — it breaks cancellations without removing
the terms that would have completed them.

This independently reproduces the paper's own counterintuitive observation that
"reducing δ does not always improve accuracy" (Gharibyan et al., abstract) on a
completely different circuit family, and gives a concrete mechanism for it.

**C5 survives, quantified: 4/18 truncated runs were inadmissible (|⟨O⟩| > 1).**
All four were at mild δ (1e-4, 1e-3) — never at aggressive δ. The check is free
and one-sided: it can prove a run invalid, never valid.

---

## STATUS 1: F2 RETRACTED, F1 NARROWED (after building real modexp)

Building the real Beauregard modexp (`modexp.py`) invalidated two of the three
findings below. **Read this block before trusting F1–F3.**

- **F2 is RETRACTED.** The "QFT is the bottleneck" numbers came from a QFT with
  a bit-reversal bug (trap 4) *and* a toy where the adder acted on the b-register
  while the QFT acted on the a-register — so the observable barely met the
  arithmetic. With the QFT fixed, the toy sandwich collapses to a constant
  `N_max = 8` at every size. The old scaling table (138 → 872306, 1.57 bits/qubit)
  is an artifact. Do not cite it.
- **F1 is NARROWED.** It holds for *gate-level* permutation circuits only.
- **F3 survives in corrected form**, see F5.

Real-circuit numbers (N=5, a=2, n_exp=3, 11 qubits; modexp = 8038 rotations /
3339 non-Clifford; δ=1e-3):

```
   observable   circuit     N_max         <O>
       Z_exp0    modexp         1   +1.000000     <- trivial: Z on a control qubit
       Z_exp0      shor     46756   +0.072884
         Z_x0    modexp     13313   +1.022676     <- |<O>| > 1 : INVALID estimate
         Z_x0      shor     13313   +0.222552
         Z_x1    modexp     17555   +1.021963     <- also invalid
         Z_x1      shor     17555   +0.162320
         Z_b0    modexp     17126   +0.455090
         Z_b0      shor     17126   +0.025288
```

### F4 — PPS cost depends on the *gate-level* arithmetic, not the logical map

Toffoli-based adder (`circuits.ripple_adder`) is a permutation **gate by gate**,
so Z-strings stay Z-type and everything cancels (F1). Beauregard modexp computes
the *same kind of logical map* but in Fourier space, where the individual gates
are CP/Rz rotations and are **not** permutations. Result: Z_x0 through the modexp
gives `N_max = 13313`, not 1.

**⟹ Same logical arithmetic, opposite PPS profile, depending only on whether you
compile to Toffoli or to Fourier rotations.** This was thread 3 (a curiosity);
it is now the main result and the main live question.

Also note: adding the inverse QFT changes `N_max` **not at all** for x/b-register
observables (13313 both, 17555 both, 17126 both) — it only shifts `<O>`. The
modexp dominates. This is the direct refutation of F2.

### F5 — Truncation produces provably invalid estimates (corrected F3)

`<Z_x0> = +1.022676` at δ=1e-3. For any Pauli observable |<O>| ≤ 1, so this is
not a slightly-wrong answer, it is an *impossible* one. Truncation error here is
not small and not bounded. This is much harder evidence than the old F3 slope
argument, and it does not depend on the toy circuits.

**Cheap validity check to keep using: assert |<O>| ≤ 1.** Free, and it caught
this immediately.

---

## Findings (F1 narrowed, F2 retracted, F3 superseded by F5 — see above)

### F1 — Gate-level permutation circuits are FREE for PPS with Z-type observables

**Scope: Toffoli/CNOT/X circuits only. Does NOT extend to Fourier arithmetic.**

3-bit Cuccaro adder, 42 T gates, observable `Z_b0`:

```
peak Pauli terms 128  →  final terms 1  →  <O> exact at every δ tested (incl 1e-2)
```

**Mechanism:** a classical reversible circuit is a permutation matrix.
Conjugating a *diagonal* operator by a permutation stays diagonal. Z-type Pauli
strings span the diagonals ⟹ Z-strings map to Z-strings, and every intermediate
branch into an X/Y string **cancels exactly**. Branching is real but transient.

**⟹ Toffoli-density does NOT imply PPS-hardness.** This is elementary once seen
and is almost certainly known. Treat as corrected foundation, not a result.

### F2 — The QFT is the bottleneck, not the arithmetic

Adding an inverse QFT on the a-register drives the observable off-diagonal, and
*then* the T gates branch for real:

```
                 Z-type weight fraction:  start   min    end    mean
arithmetic only                           1.000  0.000  1.000  0.484
+ inverse QFT                             1.000  0.000  0.000  0.004
```

Scaling (δ=0, exact PPS; "sandwich" = H layer → adder → inverse QFT):

```
 nbits qubits Tgates | arith N_max fin | sand N_max    fin  ratio
     2      6     31 |          32   1 |        138     52    4.3
     3      8     51 |         128   1 |       1158    384    9.0
     4     10     74 |         512   1 |      10306   2916   20.1
     5     12    100 |        2048   1 |      94006  22592   45.9
     6     14    129 |        8192   1 |     872306      -  106.5 (CAP)

sandwich: log2(N_max) ~ 1.57 bits per qubit   (1.0 = 2^n, 2.0 = 4^n worst case)
```

arith column *always* ends at 1 term (F1). Sandwich never collapses, and the
QFT's multiplier is itself growing exponentially.

**⟹ For Shor-like circuits the PPS bottleneck is the QFT, not the modular
exponentiation — the opposite of where the gate count sits.** This is the most
interesting thing found so far.

### F3 — δ is not a working dial for this family

```
nbits=5 (12 qubits):  δ=1e-2  N_max= 2974   slope 2.09
                      δ=1e-3  N_max=41398   slope 1.14
                      δ=1e-4  N_max=91446   slope 0.34
                      δ=1e-6  N_max=94006   slope 0.01
                      δ=0     N_max=94006   (exact)
```

Below ~1e-4 δ does nothing: coefficient spectrum has a **floor**, not a
power-law tail, so truncation has nothing left to discard. Accuracy fails as a
**cliff**, not gracefully — on a non-degenerate observable: exact through
δ=3e-2, then error jumps to 1.0 at δ=1e-1.

Matters because the paper's practical contribution (extrapolate N_max from cheap
test runs, their Eq. 17 `N_max ~ δ^-m`) needs a power law to extrapolate along.

---

## Falsified hypotheses — do not re-derive

**H1 (DEAD): "arithmetic is hard for PPS because all T angles are π/4, so
cos=sin=1/√2 and nothing is safe to truncate."**
Branch-weight data *supports the premise* (adder: 42 branching gates, all
|sin|=0.7071 exactly, 1 distinct value; brickwork: 96 gates, 96 distinct values,
mean 0.377). But the conclusion is wrong — see F1, it all cancels.

**H2 (DEAD): "the π/4-uniform angle is what breaks the power law."**
Paper's **Appendix C / Fig. 9b** already tested fixed correlated θ_X = π/4 and
the power law *still held*. The angle is not the differentiator. If there's a
real effect it's the arithmetic+QFT **structure**, not the angle value.

---

## Traps already hit (do not repeat)

1. **Floating-point residue counted as real Pauli terms.** Exact cancellations
   leave ~1e-16 junk. `pps.py` now floors at `abs(v) > 1e-13` even when δ=0.
   Before the fix, `experiment.py` §3 was measuring pure fp noise (coefficients
   reported at 2^-53, 2^-105 — garbage). **`experiment.py` §3 output is stale/
   invalid**; §1 and §2 are still fine.
2. **QFT bit-reversal acts on the *input* index (columns), not rows.**
   `U ≈ phase * F[:, bitrev]`. Cost an hour. `QFT∘QFT⁻¹=I` passing does *not*
   validate the convention.
3. **Degenerate observables.** In the H→adder→QFT toy, every a-register
   observable has `<O> = 0` by symmetry (an adder has no period structure). The
   only non-zero ones found were ancillas that return to |0>. **Cannot do a real
   accuracy-vs-δ study without genuine modular exponentiation.** This is the
   main reason the modexp generator was built — and building it immediately
   killed F2, so the instinct was right.
4. **QFT swap network goes BEFORE the body, not after.** The H/CP body puts a
   bit reversal on the *input* index (`U = phase * F[:, bitrev]`). Forward QFT is
   therefore `Body . Swaps` — swaps applied to the state first. Getting this
   backwards silently produces a unitary that still satisfies `QFT∘QFT⁻¹ = I`
   **and still passes a casual eyeball test**, but breaks Fourier arithmetic
   (bit reversal does not commute with addition) and *silently changed all of
   F2/F3*. Costliest bug of the session. `circuits.qft(..., swaps=True)` is now
   verified against the plain DFT with sign +1.
5. **(SUPERSEDED — environment changed.)** Run everything with `uv run python`
   from `research/` (Python 3.14 + `.venv`, numpy 2.4.6, scipy, stim, qiskit,
   quimb, numba, juliacall). Only `source .env` when Julia/PauliPropagation.jl is
   needed — its `LD_PRELOAD` of Julia's libstdc++ **segfaults numpy longdouble**
   (exit 139), which cost a debugging cycle. The old advice to use
   `/usr/bin/python3` predates the venv and no longer applies.
6. **`pkill` before a heredoc in the same command kills the write.** Two scripts
   vanished this way. Write the file first, kill second — or use the Write tool.
7. **CPython 3.14 crashes on long perm_pps runs.**
   `Fatal Python error: _TAIL_CALL_CACHE: Executing a cache.` — a bug in 3.14's
   new tail-calling interpreter, not in this code. It killed two background jobs
   mid-run and looked like a hang or a stray `pkill`. If a long run dies with no
   traceback, suspect this first. Workaround: rerun (it is intermittent), or use
   a 3.12/3.13 interpreter for long jobs.

---

## Caveats that decide whether F3 is real

- **Scale gap is severe.** Mine: 6–14 qubits, 31–129 T gates. Paper: 127 qubits,
  5000–8000 gates. Paper says the power law only emerges after ~1/3 of the
  circuit. My circuits may simply be **pre-asymptotic**, which would make F3 an
  artifact. This is the single biggest threat to the result.
- Python dict-based PPS dies around ~15 qubits. Paper's **Appendix B** gives the
  bit-packed representation (ν_P vectors in uint64 arrays) — that's the fix if
  scale is needed.
- `<O>=0` degeneracy above means F3's cliff was measured on a trivial observable.

---

## F6 — C3 RESOLVED: Z-closure is exact, but the benefit is ~8x, not 10^4

`toffoli_arith.py` built and verified (`test_toffoli_arith.py` A–G pass,
including a cross-check that both compilations give identical `a^e mod N`).
The matched A/B, N=5 a=2 n_exp=2, observable Z_x0, **δ=0 (exact)**:

```
compilation  qubits  rots  nonCliff   N_max  N_final  Z-type  non-Z  max|non-Z|
   Fourier       10  5359      2262  131064   131064     508 130556    4.71e-02
   Toffoli       15 21247      4074   40774    15476   15476      0    0
```

**Z-closure is real and exact.** The Toffoli compilation ends with *zero*
non-Z-type terms — not "small", exactly zero. The Fourier compilation ends with
130556 non-Z terms carrying real weight (largest 4.7e-2, so not fp noise).

**Why Fourier lacks closure even though it computes a permutation:** Beauregard's
circuit is a permutation only *on the valid subspace* (b=0, anc=0, x<N). As a
full 2^n unitary it is not a permutation matrix — the phase rotations only
conspire on that subspace. Toffoli/CNOT/X is a permutation matrix on the whole
Hilbert space, unconditionally.

**But the payoff is modest:** 15476 vs 131064 final terms (~8.5x), N_max 40774 vs
131064 (~3.2x). Not the "four orders of magnitude" the draft abstract claimed.

## F7 — the real cost driver is Walsh sparsity, which is compilation-invariant

Z-closure bounds support to 2^n Z-strings instead of 4^n Paulis — a genuine
quadratic ceiling. It does **not** make the problem cheap. The Toffoli modexp
still needs 15476 Z-strings.

Reason: for a permutation π and observable Z_j, `π† Z_j π` is the diagonal
operator `(-1)^{f(x)}` where f is the Boolean function giving bit j of the
output. Its Pauli expansion is exactly the **Walsh–Hadamard expansion of f**, so
the term count is the Walsh sparsity of f.

- ripple adder: output bit = `a0 XOR b0 XOR c0`, **linear** ⟹ 1 Walsh coefficient
  ⟹ N_final = 1. That, not "permutation", is why F1 collapsed.
- modexp: bit of `a^e mod N` is **highly nonlinear** ⟹ dense Walsh spectrum
  ⟹ 15476 terms.

**⟹ Compilation sets the ceiling (2^n vs 4^n); the algorithm's Boolean structure
sets where you sit under it.** The draft abstract's thesis ("compilation, not
algorithm") is therefore wrong as stated and needs rewriting — it is *both*, with
distinct roles. F1's original explanation was also wrong (right conclusion,
wrong mechanism).

## F8 — truncation destroys the structural advantage

At δ=1e-3 the ordering **reverses**: Toffoli N_max 24692, Fourier 11581. The
compilation that is better exactly is worse under truncation, because δ-truncation
discards terms that were going to cancel, so the exact Z-closure cancellation
never completes. Norm violations appear on both (⟨O⟩ up to +1.109).

Practical reading: telling someone "compile to Toffoli before running PPS" is
**only** sound advice at δ=0, which is not a regime anyone runs in. This
substantially weakens the applied claim.

---

## Open threads

1. ~~Real modular exponentiation~~ **DONE** — `modexp.py`, verified end to end
   (`test_modexp.py` A–G all pass; N=5,7,15; period peaks land exactly).
2. ~~Re-run F2/F3 on real modexp~~ **DONE** — F2 did not survive. See STATUS block.
3. **[NOW THE MAIN THREAD] Toffoli-arithmetic vs Fourier-arithmetic A/B.**
   F4 says the compilation choice, not the logical function, sets PPS cost.
   Needs: a Toffoli-based modular multiplier to put beside the Beauregard one,
   so the comparison is like-for-like on the *same* N, a, and observable.
   Currently only have a Toffoli *adder* (non-modular) — that is the next build.
   Prediction to test: Toffoli modexp gives small `N_max` (permutation ⟹ Z-type
   preserved), Fourier modexp gives large `N_max`, identical logical circuit.
   If that holds cleanly it is a genuine, checkable statement about when PPS is
   applicable, and it is actionable (compile arithmetic to Toffoli before PPS).
4. Accuracy study now actually possible: non-degenerate `<O>` exist on the real
   circuit. Sweep δ, plot error vs `N_max`, use the |<O>| ≤ 1 check as a
   validity gate. Expect a cliff, but measure it properly this time.
5. Scale test for the pre-asymptotic worry (needs bit-packed PPS, App. B).
   Current Python dict PPS: ~15s for 8000 gates at 11 qubits, N_max ~50k.

6. **[HIGHEST VALUE, not started] Weight-truncation as Fourier tail mass.**
   The Pauli weight of `Z^z` is `popcount(z)`, which is exactly the Fourier
   degree of that coefficient. So for permutation circuits, **weight-truncation
   error is literally the Fourier tail mass** of the pulled-back bit function.
   The current work covers only coefficient-truncation (δ); weight-truncation is
   the other standard PPS knob and is completely untouched. Measuring
   mass-by-weight profiles (adder vs modexp vs the binomial profile of a random
   function) would extend the exact cost model to it. Logged in
   `experiment_review.py` R3.
7. **Prove the circuit-level C15 invariance.** The valid-subspace argument does
   not cover it (added qubits are live). Most promising route: is the support
   confined to an affine subspace? That would also explain density → ½ exactly.
   Paper B open problem 1.
8. Decode *which* qubits carry the 4 dominant Walsh coefficients that reproduce
   ⟨O⟩ exactly (F13) — would turn the aggressive-truncation rule from empirical
   to structural.

## Honesty log

Things believed and then killed, in order. Keep adding to this.

- H1 "π/4 angles break truncation" — killed by F1 (it all cancels).
- H2 "uniform angle breaks the power law" — killed by paper App. C.
- F2 "the QFT is the PPS bottleneck in Shor" — killed by the real modexp
  (bit-reversal bug + a toy where the observable never met the arithmetic).
- F4/F6/F8 "Fourier lacks Z-closure; Toffoli is 8.5x better; truncation
  reverses it" — **all killed by the θ=π PPS bug.** Both compilations have
  exact Z-closure; Fourier is cheaper, not dearer; the ⟨O⟩ values were wrong.
- F1's *explanation* — right conclusion (adder collapses to 1 term), wrong
  mechanism (permutation-ness). Real reason: the output bit is affine.
- F7 "Walsh sparsity is the driver, compilation sets the ceiling" — half right.
  Walsh sparsity is exactly the driver (F9), but there is no separate
  compilation ceiling; the apparent one was an ancilla-count confound.
- C17 "permutation-native PPS halves peak memory **exactly**" — **rounding
  artifact**, corrected 2026-08-08. The test table printed `%.1f`, so 1.9997
  displayed as "2.0x" and the phrase "exactly" was written from the display
  rather than the numbers. True: 2.000000 for adders, 1.9997 for modexp with
  `rot = 2·perm − 2` in 3/3. **The lesson is narrow and worth keeping: never
  write a precision claim from a rounded display.** Caught by external review
  noticing the abstract claimed more than the section delivered.
- C29 "V²=id is the entire condition" (§G, TODO 9) — **over-claimed, narrowed
  by §WD.** Not wrong about anything it tested; wrong about what it had tested.
  Every synthetic block in step 9 was controlled on exactly one fresh qubit, so
  the experiment could not separate "V is an involution" from "the control is a
  single qubit". `SelectModExp` has the first and not the second, and loses the
  invariance. **The generalisation was stated at the altitude of the proof
  sketch rather than of the evidence** — the proof genuinely uses only V²=id
  *given* the block form, and that conditional got dropped in the write-up.
- P5 (TODO 12) "restoring the j=0 branch makes the tail dead" — **half
  refuted, and the half that failed was informative.** Support flatness held
  exactly; the dead-bit half failed because the activation ancilla is itself a
  scratch qubit the pullback ranges over. Logged rather than re-scoped; the
  diagnosis is in §WD and is now a check in `experiment_windowed.py`.
- C24's tail-confinement test at |I| = 1 (TODO 12e) — **not wrong, but
  vacuous, and it had been counted as evidence.** "z_I ∈ {0, 1_I}" has one
  possible answer when the tail is a single bit, so the n_exp = α+1 rows say
  nothing and a β>1 control evaluated there *passes*. Caught only because the
  control passed. Re-scoped to |I| ≥ 2; no number changed, three of six rows
  stopped counting. See §OS3. **The general shape: a test can be sound,
  correctly implemented, and still carry zero information at some parameter
  values — and those are exactly the values a sweep starts at.**

**Lessons, in order of how much they cost:**
1. The toy sandwich (H → adder → QFT) is not a proxy for Shor. Use `modexp.py`.
2. **Randomised tests only cover the gates they sample.** The θ=π bug survived
   because no random circuit ever emitted an X. Enumerate the gate *set*, don't
   sample it.
3. **An independent exact reference is worth more than more tests.** Walsh
   caught what six suites missed, because it computes the same quantity by a
   completely different route.
4. **Precision-sweep to separate bugs from float error.** Same error in
   `longdouble` as in `float64` ⟹ it is a bug, full stop.
5. Matched-instance A/B needs matched *qubit counts*, not just matched logical
   function.
6. **State a generalisation at the altitude of the evidence, not of the proof
   sketch.** C29 said "V²=id is the entire condition" on the strength of blocks
   that all had exactly one control. When every instance you tested shares an
   incidental property, that property is a *hypothesis you did not vary* — list
   it before generalising. The cheap check is to ask what the must-fail control
   would have to look like to break the clause you are about to drop; here it
   was "a block with V²=id and a non-affine control", which took an afternoon
   to build and refuted the claim immediately.
7. **The pullback sees scratch qubits in every state, not just |0⟩.** Two
   constructions identical on the valid subspace can differ completely in Walsh
   support — that is the whole content of §WD, and it is also what refuted half
   of P5. Reason about the permutation on the *full* space, always.

---

## Literature anchors

- Gharibyan et al., PPS practical guide — [arXiv:2507.10771](https://arxiv.org/pdf/2507.10771).
  Eq. 8/9 = branching rule. Eq. 11 = power law. Eq. 17 = N_max. App. B =
  bit-packed rep. App. C = correlated angles (kills H2). App. E = power-law
  deviations. App. F = trouble estimating m.
- Dang, Hill & Hollenberg, MPS Shor, 60 qubits — [arXiv:1712.07311](https://arxiv.org/abs/1712.07311).
  Memory depends on the *factors of r*, not r.
- Begušić & Chan, sparse Pauli dynamics vs IBM 127q — [arXiv:2306.16372](https://arxiv.org/pdf/2306.16372).
- Orús & Latorre, entanglement in Shor scales **linearly** in n — [quant-ph/0311017](https://arxiv.org/pdf/quant-ph/0311017).
  (Efficient MPS needs O(log n). That gap is the whole story.)

## File map

Run everything with `uv run python` from `research/` (see trap 5).

**2026-08-08 restructure:** experiment scripts moved to `experiments/`
(filenames unchanged; run as `uv run python -m experiments.<name>`), fully
retracted ones to `archive/`, and the shared machinery was extracted into the
`lab/` package (`lab.harness`, `lab.measure`, `lab.gf2`, `lab.modarith`,
`lab.variants`, `lab.nulls`) with `test_lab.py` pinning it to the numbers
logged here and `test_claims.py` re-verifying the headline claims. Bare
`experiment_*.py` names below predate the move.

```
perm_pps.py      permutation-native PPS: X/CNOT/Toffoli as atomic gates,
                 stays Z-type throughout, 2x lower peak, ~10x faster
walsh.py         classical_permutation, FWHT, pullback_coefficients,
                 walsh_sparsity, is_affine, permutation_via_statevector
toffoli_arith.py Toffoli-compiled modexp (Cuccaro + add/sub-N reduction)
windowed_arith.py windowed modexp (§WD): WindowedModExp = table lookup +
                 register multiply; SelectModExp = select-multiply with the
                 skip_zero knob; replay()/verify_modexp() gate correctness on
                 circuits far too wide for a permutation array
pauli.py         symplectic algebra, rotation rule
circuits.py      Circuit (σ,θ list), gate decomps, Cuccaro adder, QFT+swaps,
                 cswap/ccphase/inverse, brickwork
statevec.py      O(2^n) per-gate state-vector sim (to_unitary dies past ~12q)
pps.py           propagate(), exact_expectation(), fit_power_law()
modexp.py        Beauregard modexp: phi_add -> cc_phi_add_mod -> cmult_mod
                 -> u_a -> build()/build_shor().  VERIFIED.
test_core.py     correctness gate — run first, always
test_modexp.py   modexp layers A-G — run second
experiment.py    branch weights + coeff distribution  (§3 STALE, trap 1)
experiment2.py   F1 mechanism + gate verification     (F2 section RETRACTED)
experiment3.py   old toy scaling                      (RETRACTED, see STATUS)
```

Order-finding sanity anchor: `ModExp(15, 7, n_exp=4).build_shor()` must give
exponent-register peaks at y = 0, 4, 8, 12 with p = 0.25 each (r = 4).
