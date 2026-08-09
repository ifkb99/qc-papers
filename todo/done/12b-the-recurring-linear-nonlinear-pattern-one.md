---
id: 12b
state: done
title: Is the recurring GF(2) linear/nonlinear pattern one theorem?
claims: [C30, C33, C34, C40, C41]
---
# Is the recurring GF(2) linear/nonlinear pattern one theorem?

**Mostly yes, and it made two of our claims cheaper to state and one of them
wrong.** See `NOTES.md` §GF; file `experiments/experiment_gf2law.py`. Claims
C40, C41; C30/C33 downgraded from "our mechanism" to "cite Carlet Prop. 29".

The unifying *fact* is the defining property of the transform — Walsh
characters of GF(2)ⁿ **are** the affine functions — with two corollaries:
**(a)** affine symmetry confines the support (**classical**: Carlet,
Proposition 29; our C30/C33 are instances, and "partially bent" is the name for
the extremal case); **(b)** affine gating costs one Walsh coefficient, which is
why §WD measured exactly ×2^w per OR-gated window.

- **C40** extends (a) to *conditional* structures: per-cell linear structures
  w_u ⟹ support avoids E = {z : w_u·z = 1 ∀u} ⟹ **density ≤ 1 − 2^−d**,
  d = dim span{w_u}. Verified 4/4 planted with 0 violations, and shown
  affine-invariant (tested, not asserted).
- **C41**: consistency is a **parity** condition — an odd dependency among the
  w_u destroys the cap outright (0.9845 vs 0.8671 at the same cell count).
- **It corrects §RS.** The conjectured 1 − 2^−(k+1) ladder is the wrong shape:
  the parameter is the span dimension, not the conditioning depth.
- **It explains C34 without modexp.** v4's two slice structures (msb⊕anc and
  msb) span d = 2, so E is exactly the quadrant {z_msb=1, z_anc=0} that C34
  found empty, with cap ¾.

**Prior art — pursued as far as open access goes; risk now low.** The one
unread item is Carlet–Tarannikov, DCC 25:263–279 (2002), paywalled at Springer
with no self-archived preprint. Checked instead: **Carlet's own book** (the
comprehensive survey by the same author, which cites that paper on pp. 205,
206, 319 and reproduces its Def. 47 / Prop. 60 — §5.5 read in full), and **the
paper's own abstract**. A covering sequence is a *single global* λ; a *partial*
covering sequence allows two levels; ours has a **different structure vector
per cell**, which is neither. The abstract's stated contributions are
resiliency/correlation-immunity characterisations and constructions — not a
support-confinement law. Also worth citing and distinguishing:
**Maiorana–McFarland** (restrictions *affine* on each coset — same flavour,
stronger hypothesis, different conclusion).

Still state C40/C41 at corollary altitude (Prop. 29 + coset decomposition) and
note they may be folklore. If someone gets institutional access, reading the
DCC body is the last loose end — a ten-minute job.
