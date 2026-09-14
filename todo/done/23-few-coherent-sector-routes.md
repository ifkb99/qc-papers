---
id: 23
state: done
title: "Can a few coherent route choices be sampled without enumerating sectors?"
outcome: "bounded coherent-history sampler implemented and audited; known gate-by-gate sampling gives a no-rejection alternative, with static regrouping and finite costs explicitly compared"
claims: [C55, C56, C58, C59]
---
# Coherent route choices, not noncommutation alone

Completed 2026-09-11 UTC at bounded first-pass scope. C59 owns the history
identity, both samplers, known-input promises and costs; §CH owns independent
evidence and retained verifier failures. The original proposal below is kept
as history, not current uncertainty about that tested construction. This
does not solve arbitrary coherent gates or numerical certification. TODO 24
owns the next question.

C58 handles a single deterministic route at each work insertion. The next
candidate is a supplied short linear combination of such routing unitaries,
with the full inserted gate promised and checked to be unitary. The connection
to coherent-component rejection in C55 is worth testing, but this is not yet a
validated sampler or novelty claim. Consult current primary literature on
linear-combination/stabilizer-decomposition simulation before positioning it.

## Candidate identity to audit before implementation

Expand k inserted gates into histories h with coefficients c_h. Each history
is a normalized routed circuit. For a chosen FINAL sector gamma, the history
must start from alpha_h=Pi_h^(-1)(gamma), not from one common initial alpha.
For initial orbit label 0, all these initial sectors carry amplitude 1/sqrt(M).
Let a_h(gamma,y) be the b-dimensional final work amplitude for output y under
that normalized component circuit, expressed in the SAME final-sector basis.
Proposed target and proposal:

    p(gamma,y) = ||sum_h c_h a_h(gamma,y)||^2 / M,
    B = sum_h |c_h|,
    q(gamma,y) = sum_h |c_h| ||a_h(gamma,y)||^2 / (M*B).

Weighted Cauchy-Schwarz predicts p<=B^2*q. If the full circuit is normalized,
the proposal draws gamma uniformly, chooses h with weight |c_h|/B, and uses
the existing finite-work sampler for that routed component. Accept with
p/(B^2*q); predicted mean proposals B^2. Audit zero-weight cases explicitly.
Uniform gamma is a PROPOSAL, not generally the true final-sector marginal.

This costs evaluation of EVERY history for an acceptance test unless a
separate compression is proved. Charge the history count H, not merely k or B.
For two-component gates H<=2^k; factorized coefficient choices can sample a
history cheaply, but evaluating its coherent sum still costs H. Verify that
complete output amplitudes use the original Fourier-filtered branch products,
with their relative phases intact; conditional probabilities alone lose them.
No second generic propagator is needed.

## Stronger baseline and honest scope

First compare exact sector regrouping. For the N=7 physical Rz example,
the two route options close on just two old sectors. Grouping them already
solves the problem; an exponential-in-k implementation would add nothing.
More generally, a short Fourier expansion is not automatically a short
unitary-route decomposition with small coefficient norm.

If a second fixture is needed to separate these costs, one explicit indexed
family uses cell reflection R|b*m+p>=|b*(-m mod M)+p> and J_q=D_q R.
J_q is Hermitian and unitary, routes alpha to -alpha-q, and
cos(theta/2)I-i*sin(theta/2)J_q is unitary. Different q can generate a large
sector orbit even when only a few coherent insertions occur. This is an
ABSTRACT known-index gate promise, not a physical few-qubit gate or general
model of hardware noise. Charge recognizing/synthesizing that promise.
Do not substitute a non-Hermitian arbitrary D_q into the cosine/sine formula
and assume the result is unitary.

## Bounded tests and decision

Send initial formula/control tests to lower-cost agents while main audits
normalization and the regrouping comparator. Fix r,b,t and background first,
vary one angle; then vary k separately. Start r<=18,t<=6,k<=3 with capped
existing full-r `sequential_path` and explicit sector projections. Enumerate
tiny joint laws to check normalization, the envelope and accepted law before
testing actual draws. Must-fail controls: delete coherent cross terms, use
the wrong initial sector for a shared final gamma, and treat uniform final
gamma as the target marginal. Retain any vacuous control and derive why.

Only after those pass consider a fixed-k, supplied-large-M draw with no
sector/output arrays and honest scalar/bit/precision accounting. If regrouping
already gives an equal or better bound in the selected family, record that
and stop the sweep. Even a positive result would be a parameterized simulation
promise, not a generic quantum-simulation breakthrough. Numerical certification
remains unresolved and must be added separately to any approximation guarantee.
