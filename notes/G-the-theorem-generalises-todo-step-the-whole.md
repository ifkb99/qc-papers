---
code: G
title: "THE THEOREM GENERALISES (TODO step 9). V²=id is the whole condition."
date: 2026-08-08
outcome: superseded
claims: [C15, C29]
todo: [12, 9]
---
# G — THE THEOREM GENERALISES (TODO step 9). V²=id is the whole condition.

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
