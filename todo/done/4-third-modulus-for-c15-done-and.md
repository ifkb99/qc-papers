---
id: 4
state: done
title: Third modulus for C15
outcome: DONE, and it caught an over-claim
claims: [C15, C20, C21, C22]
---
# Third modulus for C15

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
