---
code: AF
title: "TODO 11.1: AFFINE structures, and the function-level ½ rows EXPLAINED"
outcome: solved
claims: [C33]
todo: [11]
---
# AF — TODO 11.1: AFFINE structures, and the function-level ½ rows EXPLAINED

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
