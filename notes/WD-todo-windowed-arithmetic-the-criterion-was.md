---
code: WD
title: "TODO 12: WINDOWED ARITHMETIC. The criterion was incomplete."
date: 2026-08-08
outcome: record
claims: [C8, C24, C29, C30, C32, C36, C37, C38, C39]
todo: [12]
---
# WD — TODO 12: WINDOWED ARITHMETIC. The criterion was incomplete.

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
