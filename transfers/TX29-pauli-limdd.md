---
id: TX29
field: "Local-invertible-map decision diagrams (Pauli-LIMDDs)"
status: open
effect: lower-bound
one_line: "Post-modexp state: exponent-first LIMDDs need > r/(2 sqrt(2N)) nodes at one level under a hypothesis on a (C110); x-first need ceil(beta/2) when t - alpha >= 3 bitlen(beta) - 1 and the set-bit positions of beta have gcd 1 (C111); both derived, interleaved orders open"
source: "Vinkhuijzen et al., Quantum 2023, Def. 3 and Thm 1 p.12 (surveyor slate 4); Vinkhuijzen-Coopmans-Laarman arXiv:2401.01322 §2 and appendix Lemma 4 (level-class counting, read 2026-09-19)"
claims: [C110, C111, C112, C108]
notes: [EF, LM]
todo: [69]
---
# TX29 — Local-invertible-map decision diagrams (Pauli-LIMDDs)

## Dictionary

Node sharing up to local Pauli ↔ equivalence classes of sub-states under XOR
translations, sign flips and scalars. For the post-modexp state at the last exponent
qubit, a class is a value of c ⊕ (Ac mod N) (C110 Lemma 1).

## Hypotheses

Exponent-first orders (C110): A = a^(2^j) ≢ ±1 mod every prime factor of N, and 2^t ≥ r.
This fails at every j when 3 | N, and for large j when N has a Fermat-prime factor.
x-first orders (C111): the set-bit positions of the odd part β of r have gcd 1, and
t − α ≥ 3·bitlen(β) − 1; for odd β < 64 a persistence lemma plus an exhaustive check
lower that threshold to T(β) ≤ 2·bitlen(β), g > 1 included. C112 gives the exact class
count at the last exponent level of any order. All three are derived, and each claim's own
file carries its review history. The row stays open because the transfer is bounded: the
bounds hold in the two extreme order classes only.
Interleaved orders are unchecked (TODO 69). Non-Pauli LIM groups are outside the row.

## Consequence for the goal

In exponent-first orders, when no prime p | N has ord_p(a) a power of 2 (C110's hypothesis
for every choice of last exponent qubit), a Pauli-LIMDD does not remove the exponential:
it has at least N^δ/(2√2) nodes once r ≥ N^(1/2+δ). Measured post-hoc at 0.79–1.00 of D_j
(⌈r/2⌉ ≤ D_j ≤ r) on 255 random-semiprime cases meeting the hypothesis, n ∈ {10, 12, 14,
16}, j ∈ {0, t − 1}. Sweeping every base of odd N < 2^9, still under the hypothesis, the
ratio falls to 0.22 (C110, Evidence). Moduli and bases outside the hypothesis (e.g. 3 | N) are
not covered. In x-first orders, under C111's hypotheses, a reduced LIMDD has exactly
⌈β/2⌉ nodes at the first exponent level, and every LIMDD has at least that many. Moving
the row to `obstruction` needs interleaved orders settled (TODO 69), or an explicit
order-class scope written into the row. These are lower bounds,
consistent with TX15.
