---
id: 6
state: done
title: The intermediate 2-adic law
outcome: DONE, negative result
claims: [C15]
---
# The intermediate 2-adic law

**There is no intermediate law.** See `NOTES.md` §I; files
`experiment_c15_intermediate{,2,3}.py`.

- β=1: sparsity exactly constant in t (4, 4, 16 across t=11..22). Matches the
  proof.
- β>1: sparsity is Θ(2^t), density bounded away from 0 — but it does **not**
  organise by α. With a fixed table and consecutive t it shows a strong
  period-`ord₂(β)` oscillation plus a slow upward drift, and neither component
  is a function of α.
- The periodicity hypothesis was refuted too: residue classes drift rather than
  staying constant (r=7, t≡2 mod 3: 0.834, 0.831, 0.850, 0.875).
- **The N=323 "intermediate" observation was sampling aliasing** — t=16,20,24
  hit residues 4,2,0 mod ord₂(9)=6, so three points of an oscillation were read
  as a trend.

Net: the binary β=1/β>1 dichotomy is the robust structure and C15 should be
stated without hedging about intermediate regimes. A negative result that
protects the paper rather than complicating it.

**Two methodology bugs of my own, caught by the protocol:** degenerate random
tables at small r (rejected constants), and drawing a new table per t, which
confounded t-dependence with table variance — a direct violation of "vary
exactly one parameter". Both produced confident-looking numbers first.

**Scope caveat CLOSED** (`experiment_c15_realfn.py`). Rechecked with the real
table `h[c]=bit_j(a^c mod N)`: P1 holds (constant sparsity for β=1), and α is
confirmed not to organise the data — α=0 (r=3) matches α=1 (r=6) exactly, while
α=1 spans r=6,10,18,22 with wildly different behaviour. Step 6's conclusion
transfers.

**New observation, unpursued:** the real modexp table is measurably non-generic
— sparser than random at r=6 (exactly 0.500 on even t, a clean factor of two)
and denser at r=10. So modexp bit functions carry structure beyond "depends on
e mod r". Random-table densities must not be used as a proxy for real ones.
