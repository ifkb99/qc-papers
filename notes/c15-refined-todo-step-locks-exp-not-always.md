---
code: c15-refined-todo-step-locks-exp-not-always
title: "C15 REFINED (TODO step 4) — it locks at n_exp = v2(r)+1, not \"always\""
outcome: record
claims: [C15]
todo: [4]
---
# C15 REFINED (TODO step 4) — it locks at n_exp = v2(r)+1, not "always"

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
