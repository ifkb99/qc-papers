---
id: TX45
field: "Computational number theory: collision aggregation and shift recurrences"
status: open
effect: unknown
one_line: "C121 completes the actual giant-step ratio audit; a short constructible cover with charged cleanup and seeds remains open"
source: "Hittmeir (2020), Harvey (2020), Harvey--Hittmeir (2021); RSA73 structural slate"
claims: [C121]
notes: [FD, FE]
todo: [73]
---
# TX45 — Aggregate candidates with their construction charged

## Dictionary and hypotheses

For residue lists U,V modulo N, product polynomials, multipoint evaluation and
gcds aggregate collisions modulo an unknown factor. Complete factoring also
requires N-only generators with a coverage theorem. Exact coincidences modulo
N and nonunit inputs need handling; a product gcd=N needs recovery/splitting.

## Hypotheses

Primary bodies inspected in the independent slates:

* [Hittmeir](https://arxiv.org/pdf/2006.16729), Definition 3.1,
  Algorithm 3.2/Lemma 3.3: revealing sets and disjointness modulo N.
* [Harvey](https://arxiv.org/pdf/2010.05450), Proposition 2.5 and
  Algorithm 4.1: factorial/progression-evaluation baselines.
* [Harvey--Hittmeir](https://arxiv.org/pdf/2105.11105), Algorithms 4.1/4.3:
  prime/semiprime input, order hypotheses, interval/residue generators and
  exact-match cleanup.

These are existing methods. Deterministic power-of-N costs remain exponential
in bitlen(N); no practical superiority over GNFS is inferred. Later leads are
the balanced-semiprime refinement in
[arXiv:2512.19076](https://arxiv.org/html/2512.19076v1), Theorem 1.1, and
[Harvey--Hittmeir 2026](https://arxiv.org/pdf/2601.11131), Theorem 1.1.
Their full proofs were not audited here; no new memory guarantee is imported.

## Proposed transfer

Blocking lists and regenerating selected tree levels trades storage for
repeated evaluation/generation. Compare structured evaluation, not pairwise
gcds; charge retained trees, integer precision, convolution and recovery.

The sharper proposal asks whether actual giant steps have a short cheaply
constructed geometric-run cover. For a unit alpha, let H(z) be the product of
(z-gamma*alpha^j), 0<=j<t. Elementary telescoping gives

    (z-gamma*alpha^(t-1))*H(alpha*z)
      = alpha^t*(z-gamma*alpha^(-1))*H(z).

Using this identity needs seeds and safe division: take a gcd before inversion,
extract a proper factor when possible, and handle gcd=N separately. Direct
initialization already visits every run element. Few runs and cheap seeds are
unproved for the full N-only generator; discovering a cover through unknown
discrete logarithms is forbidden. Arithmetic progressions are not geometric
runs.

## Consequence for the goal

C121 owns the completed symbolic ratio/seed audit, its scoped counterexample
and the missing constructive lemma. The proposed general short cover remains
open; the audit supplies no demonstrated factoring-memory improvement or generic
factoring obstruction. TODO73 owns subsequent lemma decisions, and FE records
the audit provenance. No numerical recurrence sweep follows from this result.
