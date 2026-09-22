---
id: TX34
field: "Quantum entanglement: Renyi-2 purity, Schmidt tails and approximate tensor trains"
status: obstruction
effect: lower-bound
one_line: "C113 gives a fidelity-robust MPS bond obstruction at the specified interleaved clean-state cut; elementary degree and Fourier bounds suffice"
source: "Schuch-Wolf-Verstraete-Cirac arXiv:0705.0292v2 p.2 eq. (4), p.3; caveat for printed eq. (3); Oseledets 2011 Theorem 2.2 pp.2299-2300; primary-source audit 2026-09-21"
claims: [C19, C107, C108, C109, C110, C111, C112, C113]
notes: [QB]
todo: [70]
---
# TX34 — Purity and Schmidt weights obstruct approximate MPS compression

## Dictionary

The normalized clean post-modexp state, split across physical qubits, gives a
Schmidt spectrum. Its second moment is subsystem purity. Equal positive
amplitudes on the modular-power graph express that purity as an ordered
rectangle count, including degenerate rectangles. C113 owns the exact state,
fixed interleaving, fidelity convention, arithmetic estimates and proof.

This differs from TX9/TX17's exact rank of an observable and TX29's exact
Pauli-LIMDD equivalence classes. Neither supplies a state-approximation bound
without a separate argument.

## Hypotheses

[Schuch et al.](https://arxiv.org/pdf/0705.0292v2), p.2 eq. (4) and p.3's
Renyi-alpha>1 discussion, relate approximation obstructions to discarded
Schmidt weight. Their conventions differ from C113's normalized squared
fidelity, so C113 derives its necessary-bond inequality directly. Exact rank
alone does not control retained weight.

**Printed-source caveat:** p.2 eq. (3) prints an unsquared Euclidean norm
bounded by twice the sum of discarded weights. Taken literally this fails:
for sqrt(99/100)|00>+(1/10)|11>, the best unnormalized rank-one Euclidean error
is 1/10, exceeding twice the discarded weight, 2/100. The source audit checked
the page image. This row does not import that literal inequality.

[Oseledets](https://doi.org/10.1137/090752286), Theorem 2.2, pp.2299–2300,
bounds TT-SVD's global Frobenius error by the square root of the sum of squared
local unfolding errors. It supplies correct error accounting across all cuts,
not polynomial access to an implicitly specified coefficient tensor. The body
was checked online; the audit discloses that it could not archive that PDF.

The arithmetic hypothesis is resolved for the specified cut in C113. Its
order-only bound works for every odd modulus with the order kept explicitly.
A separate finite-Fourier argument also makes C107's full-frequency small-kappa
hypothesis sufficient; C109 supplies applicable semiprime conditions. This uses
interval occupancy of the complete subgroup and does not transplant the scalar
OBDD theorem. No additional fourth-order character estimate is needed here.

[Dang–Hill–Hollenberg](https://arxiv.org/pdf/1712.07311v4), sections 4 and
5.1–5.2, provides a relevant Shor sampling baseline with a different ordering,
a work-register qudit, and measurement collapse. Its representation factors
are not a comparison at C113's physical-bit cut and full-state contract.

## Consequence for the goal

The transfer is an obstruction for approximation at the specified cut under
C113's conditions, including a named unbounded prime family. The theorem also
covers composite moduli when its stated order or character-sum bound is useful;
it makes no claim about every semiprime base. All orders, output-only sampling,
generic gate-by-gate simulation, and efficient construction of an approximant
remain outside this result. TX15's order-finding barrier is unchanged.

TODO 70 closes with a proof and source audit, without a numerical diagnostic.
Note QB records the independent authorship, fresh reviews and integration.
The audit is bounded; no novelty or exhaustive prior-art claim is made.
