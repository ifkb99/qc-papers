---
id: TX50
field: "Abelian symmetry decomposition and output-only simulation"
status: obstruction
effect: barrier
one_line: "Sector sampling needs only d=r/|<a> cap H| and chi(a^d); any interface supplying them is quotient order finding, so known N-only H with cheap membership and phases save at most 2 or the smooth part of r"
source: "memory74 deriver D2 projector argument; C52/C116; extensions of TX18/TX36/TX37"
claims: [C52, C115, C116, C124, C125]
notes: [FG, GR]
todo: [74]
---
# TX50 — Retain phases while sampling symmetry sectors

## Dictionary

For clean modular multiplication on units modulo N, a supplied abelian
subgroup H acts by multiplication and commutes with every arithmetic block.
Exponent-only output observables also commute with that action. Decompose
work |1> into character sectors, sample a sector with its appropriate weight,
and propagate only that sector while retaining character phases on coset
transitions. This targets an output marginal, not a compact full pure state.

## Hypotheses

The proposed elementary projector argument uses the equal sector weights of
|1> under the regular subgroup action and orthogonality of sectors. It applies
only at clean block boundaries (C115/TX37), with commuting exponent-output
instruments. Arbitrary dirty ancilla gates or later work measurements are not
covered. C52 is the limiting fully diagonalized case; it requires the order.
This extends C116's clean modular structure, not its existing theorem.

An actual algorithm needs an N-only subgroup constructor, uniform character
sampler, coset representatives and efficient phase evaluation. Order discovery,
discrete logarithms and cyclotomic precision are not free. Compare the reachable
quotient orbit with the original orbit r, not with the whole ambient unit group.
For example, H={+1,-1} at N=7,a=2 does not reduce its three-state reachable
orbit. Large ambient symmetry alone proves no useful saving.

The target is reduced state storage plus subgroup/interface and precision
workspace at comparable time. Charge a trajectory per requested output, not
just the retained output bits. A finite-precision version needs a global
instrument-error budget for the exponent marginal. No full cost bound or
RSA-applicable family is established. In particular, choosing H=<a> already
requires the order-aware information that TX15 says to charge.

## Consequence for the goal

A more speculative representation route: remove redundant sectors from the
requested observable while retaining the phases that produce interference.
Unlabelled coset merging is invalid. First supply one explicit family with a
cheap subgroup interface and a smaller reachable quotient; otherwise stop at
a constant-factor or circularity obstruction. TODO74 owns that decision.

## Hypothesis check (2026-09-22): obstruction, scoped

C124 carries the result. The scope is the clean logical code at u_a block
boundaries, with exponent-only instruments and exact phases. The interface
hypotheses are (H2) membership in H, plus (H3) phase access or (H3') a known
multiple of exp(H).

Under (H2), a quotient walk gives d. Under (H2)+(H3), a scalar sampler with the
same interface draws the exact outputs without the sector state (dominance).
(H3') with an unfactored multiple does not by itself give that sampler, but a
sector method cannot run without the phases (H3). The sector representation
therefore saves |<a> cap H| over the explicit orbit and nothing over the charged
quotient-walk baseline. An N-family with poly-time (H2) together with (H3) or
(H3'), and d <= poly(log N) on a set of a of density at least 1/poly, admits
randomized polynomial-time factoring (C124; the TX15 barrier). Among the known
N-only subgroups with cheap (H2) and (H3), {+1,-1} saves at most 2, and smooth
torsion G[k] saves gcd(r,k), which divides the smooth part of r.

J_N, QR_N, <b> for small b, G^k and factor-base subgroups fall outside the
hypotheses (membership or phases are not known to be cheap). They are not ruled
out. That no such H does better is C124's conjecture. C125 shows it cannot be
proved unconditionally and reduces it, for factored-multiple phases, to the named
assumption REA_fact: equivalent under weak phases, sufficient under point
phases. Dirty scratch, work-register operations, noise and finite-precision
phases are out of scope.
