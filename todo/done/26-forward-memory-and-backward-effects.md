---
id: 26
state: done
title: "Can forward mixing remove stored history without losing rare-output correlations?"
outcome: "C66 proves the shared-state output-law contract; full-law controls and bounded product-contraction tests completed, with quantitative uniform-phase usefulness deferred to TODO 27"
claims: [C56, C60, C63, C64, C65, C66]
---
# Forward memory loss is not permission to reset the backward effect

Completed at first-pass scope in C66 / §FM. The original pre-measurement plan
below is retained as history; current next actions live only in TODO 27.
No production approximation or storage-saving benchmark is claimed.

This is the next bounded experiment after §OB. The motivation is stronger
than another block-size or known-order size record: identify when approximate
sampling can forget old forward-state information despite correlations with
later measured outputs. Keep the distinction between matrix count, bit storage,
input description, arithmetic work and final-law error explicit.

## First establish the correct approximation contract

For a fixed deterministic route history, choose one normalized PSD forward
state rho_tilde_i for each depth i, independent of the measured output prefix.
Retain the EXACT backward effects E_a and QFT branch operators K_(i,z)(a).
Audit the candidate error bound

    TV(approximate full output law, ideal full output law)
        <= min(1, sum_i ||rho_i-rho_tilde_i||_1).

Main and a lower-cost agent have checked its algebra, but no numerical
experiment or implementation for it exists yet. The specific completeness
step should be written explicitly, not waved away as generic adaptivity:

    sum_z K_(i,z)(a)^dag E_a K_(i,z)(a)
        = (B_i0^dag E_a B_i0 + B_i1^dag E_a B_i1)/2.

The feedback phase cancels in this sum. Inducting from the terminal identity
then gives a complete POVM over all (a,z) at fixed i. Applying measurement
trace-norm contraction to the SAME Delta_i at every node and using C64's
true-mass normalization inequality suggests the bound above. Approximate
earlier kernels reaching zero ideal nodes are charged by the same telescoping
argument. This is a stronger correlated-state approximation promise than
arbitrary unrelated per-node mass errors; it does not contradict C64's tree
coefficient or improve that oracle without additional structure.

Make a tiny full-law test using the existing b=3 contractions and PSD mixtures
rho_tilde_i=(1-epsilon)*rho_i+epsilon*I/b, varying only epsilon. Derive and
register predictions before measuring. Check the complete output law, all
normalizations, POVM completeness and a nontrivial trace-norm upper bound.
An outward Frobenius-to-trace bound may suffice for the first test; do not
call a floating eigensolver result an interval certificate. Keep numerical
kernel error separately charged if using finite arithmetic.

Required negative controls: resetting E instead of approximating rho;
concluding decoupling of the joint state from a mixed work marginal; and
allowing unrelated prefix-dependent rho_tilde_(i,a) while using a single-state
POVM bound. A counterexample must target the claimed implication, not merely
compare different measurement depths. The existing binary scalar-tail fixture
can hide state sensitivity, so it is a useful null control, not the main case.

## Then test whether an inexpensive approximation actually exists

The forward maps are unital channels. For a fixed repeated b=3 gate schedule,
can products contract traceless differences enough that rho_i can be replaced
by I/b after a short warm-up? Treat both rapid contraction and a conserved
mode as live hypotheses. A maximally mixed fixed point alone proves neither
mixing nor decoupling. A one-step contraction coefficient can be one while
a multi-step product contracts; test products before rejecting the direction.

Use lower-cost initial algebra/edge testers while main audits. Begin with
tiny exact-input fixtures and existing density-channel helpers, not a new
propagator or a random gate ensemble unrelated to the actual circuit. Charge
the cost of checking any contraction certificate across sectors/histories;
sampling a few phases is not a uniform large-order bound. Preserve C56's
special cancellations and C65's dyadic null control.

A summable error tail could make the number of retained FORWARD checkpoints
independent of total exponent width at fixed accuracy. It would not remove
the gate description, backward effect, precision-dependent scalar storage,
coherent-history acceptance, or input order/index costs. Recomputing windows
can exchange storage for work; count both. Changing the production sampler is
deferred until the law bound and a useful regime survive their controls.

## Primary sources and stopping rule

Main has read the centered contraction definition, Proposition 2.2's
trace-norm argument, Proposition 2.4's product bound, and the full good-block
condition/proof of Proposition 3.9 in
[Pathirana, arXiv:2605.00157v1](https://arxiv.org/html/2605.00157v1).
These are directly relevant prior results on deterministic channel products;
do not present their replacement/mixing mechanism as our novelty. Its MPS
application section remains to be read if making that comparison. No random
environment hypothesis is supplied by our fixed arithmetic circuit.

First settle the POVM/full-law test. If it fails, find the missing assumption
before a mixing sweep. If it holds but no cheap contractive regime follows,
record that specific boundary and stop the sweep. Do not expand to arbitrary
noise models, a general tensor framework or a publication claim merely to
avoid an informative negative result. TODO 24's backend-wide constant remains
deferred separately.
