---
id: 74
state: open
title: "Attack memory scaling through Grover precision, streamed counting and sampled symmetries"
outcome: "Grover coin solved for known M, given H_elem (C122); TX50 decided as an obstruction (C124); poly-memory factoring methods are the bar"
claims: [C52, C60, C117, C118, C119, C120, C121, C122, C123, C124, C125]
---
# 74: Memory scaling after the cache pilot

Note FG owns the phase-0 slate, and note GR owns the round that decided it.
C122 reduces the Grover precision question to predicate interfaces, given H_elem
at binary t. TODO72/73 stay preserved. C120 is still an engineering side result,
and C121 closes only its own audit. Moderate slowdown is acceptable, with no
numeric budget given. This item lists derivations; it does not authorize runs.

**The bar (TX51).** For the factor output, ECM (heuristically subexponential)
and Pollard rho (heuristic N^(1/4)) already factor with polynomial memory. Any
proposal here must say which exponential it removes, in which quantity, for
which output, relative to them. A memory saving at unchanged exponential time
does not meet the bar. For C52's exponent-output law, the bar is the charged
exact-order route (TX51 (ii)).

**Decided: TX50 is an obstruction** (C124). Under membership plus phase access
(the interface a sector method needs to run), subgroup sectors save nothing over
the charged quotient walk. The known N-only subgroups with cheap membership and
phases save at most 2 or a divisor of the smooth part of r. J_N, QR_N, <b>, G^k
and factor-base subgroups fall outside those hypotheses and are not ruled out.
C124's conjecture is **an assumption, not a research target** (C125, note CJ).
It implies factoring is not in BPP under weak phases (under point phases, that
factoring and F_p^* discrete logs are not both in BPP), so it cannot be proved
unconditionally. For
factored-multiple phases it reduces to the named assumption REA_fact. Lemma P
(whether poly-time phases yield a factored multiple of |H|_R) is the open gap;
pursue it only with a concrete implementation class in hand.

Even a counterexample would, by dominance, only shorten the charged quotient
walk. d <= poly(log N) on a density-1/poly set of a suffices for randomized
polynomial-time factoring (C124's Circularity paragraph); with exponential d the
operative statement is dominance, not factoring.

**Grover line.** TX48's coin is settled for known M (TX48's family has M=F_n),
given H_elem at binary t. Two steps are **dropped**:
* "If the exact expected-cost lemma holds, propose a bounded implementation
  design": the coin is standard, and implementing it would not touch the goal,
  because the remaining cost is the family-specific interface.
* The TX48 short-rational TV-surrogate derivation: the exact coin (given H_elem
  at binary t) makes it unnecessary for the coin, and it matters only for
  approximate interfaces, which no current item needs.

TX49 stays open and deprioritized; it asks about a classical counting algorithm.

For factor predicates:
* Contract A with t in the input is, given H_elem, polynomially equivalent to
  factoring (C122, an instance of TX42).
* Contract B on the two-orientation predicate, at t in the input with
  sin^2(2t theta)>=1/poly, is, given H_elem, equivalent to computing p xor q
  from N. By C123 that factors N heuristically, on average over random balanced
  semiprimes.

**Deprioritized, not dropped.** TX51's end-to-end audit of a published
low-space order routine. It matters only for C52's exponent-output law, because
for the factor output ECM and rho already set a stronger bar.

**Leads recorded, not scheduled.**
* The branch-and-prune tail in C123. Worst-case tracked size is unknown and
  plausibly driven by max(v2(p+q), v2(p-q)). Large-valuation instances fix many
  low bits of p up to a small ambiguity, which is the setting of Coppersmith-type
  known-bits methods. This matters only for the side question of factoring from
  (N, p xor q).
* Contract B over the *restricted* factor predicates (note GR) needs no p xor q
  step.
