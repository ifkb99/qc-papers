---
code: c15-mechanism-solved-todo-step-the-parity
title: "C15 MECHANISM SOLVED (TODO step 3) — the parity reduction"
outcome: solved
claims: [C15, C21]
todo: [3]
---
# C15 MECHANISM SOLVED (TODO step 3) — the parity reduction

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
