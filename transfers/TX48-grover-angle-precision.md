---
id: TX48
field: "Exact random-bit sampling and certified numerical arithmetic"
status: imported
effect: removes-exponential-in-subfamily
one_line: "The large-iteration Grover angle is one exact lazy coin (expected poly(n+log t) for known M, given H_elem); output hardness sits in the predicate interface"
source: "Brassard-Hoyer-Mosca-Tapp quant-ph/0005055 section 2 equations (5)-(8) and Theorem 2; Stoudenmire-Waintal arXiv:2303.11317 sections III.B-C; Brent JACM 23(2) 1976 and Brent-Zimmermann ch. 4 (H_elem, not body-audited); FLINT ball semantics (inclusion only)"
claims: [C60, C117, C119, C122]
notes: [FG, GR]
todo: [74]
---
# TX48 — Amplification-scale structured Grover without expanded fractions

## Dictionary

Extend TX35/TX38's explicit path predicate and final local Hadamards to an
iteration count t supplied in binary. The new question is precision/storage
at large t, not construction of another oracle or a new MPS principle.
TX38's imported fixed-t exact result remains unchanged. This extension was open
when proposed; the hypothesis check below imports it.

For D=2^n, p=M/D and W(y)=sum_x f(x)(-1)^(x dot y), the deriver proposes the
factorization P(0)=cos^2(2t theta), theta=asin(sqrt(p)), and, when the nonzero
event has positive probability, P(y|y!=0)=W(y)^2/[M(D-M)] for y!=0.
It follows by substituting the two-dimensional rotation into C117's amplitude
formula and applying Parseval. This identifies one iteration-dependent coin;
the fixed nonzero law can use the marked-state Fourier sampler with zero
rejection. Its construction must still be charged.

## Hypotheses

[Brassard et al.](https://arxiv.org/pdf/quant-ph/0005055), section 2,
equations (5)-(8), supplies the rotation in a known good/bad subspace. Here
the explicit predicate supplies M and C117's suffix environment. The source
does not supply cheap classical counts or samplers for arbitrary predicates.

Two different proposed contracts remain separate:

* Lazy exact coin: certified intervals compared against progressively revealed
  unbiased random bits. Prove refinement and bit-cost bounds. The target is
  expected polynomial workspace/time in n+log(t+1); worst-case precision and
  time are unbounded. There is no implemented large-t exact sampler yet.
* Certified TV approximation: compute short rational coefficients defining one
  approximate state, then sample its normalized law with exact integer branch
  weights. Prove a global error and precision budget before choosing a cap.
  C60's raw-state error principle is relevant; its separate sparse-block
  implementation is not automatically a certificate for this sampler.

[FLINT's ball-arithmetic documentation](https://github.com/flintlib/flint/blob/main/doc/source/using.rst),
“Ball semantics” and “Quality of enclosures”, supplies inclusion guarantees,
not a polynomial-cost theorem for this proposal. Charge square roots, powering,
rounding, interval replay, random bits, suffix environments and retained outputs.
Exact rational rejection draws have finite expected cost and unbounded
worst-case draw time. Unconditional TV does not guarantee relative accuracy
after rare postselection.

## Hypothesis check (2026-09-22)

C122 carries the result. The angle/fixed-law split holds for every predicate, and
the exact lazy coin has expected polynomial cost. That cost depends on two named
conditions: M must be known, and H_elem must hold at binary t. Polynomial t needs
no H_elem: the coin is then an exact Chebyshev rational in M/D, so M must still be
known (at t=0 the coin is exactly M/D). The TV-surrogate
contract is not needed for the coin; it stays relevant only for approximate
interfaces. FLINT supplies inclusion guarantees, not the cost theorem.

## Consequence for the goal

This removes C117's materialized coefficient growth with amplification iterations,
for known M, given H_elem at binary t. It covers coin cost only: the subfamily is
the explicit path predicate, whose Fourier-mass interface C117 constructs. No
generic oracle or factoring shortcut follows.

For factor predicates:
* Contract A (computational basis) with t in the input is, given H_elem,
  equivalent to factoring (C122, an instance of TX42).
* This row's own contract, the final H (Contract B), carries no factor
  information under canonical ordering.
* With both orientations, at t in the input with non-negligible
  sin^2(2t theta), it is, given H_elem, equivalent to computing p xor q (C122). That factors
  N heuristically and on average (C123).

No implementation is scheduled; TODO74 records why.
