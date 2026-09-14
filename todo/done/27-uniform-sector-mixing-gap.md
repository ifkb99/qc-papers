---
id: 27
state: done
title: "Can a uniform sector-phase gap make forward-history compression useful?"
outcome: "C66 certifies a uniform twelve-step gap; C67 implements exact checkpoint savings while mixing warm-up remains impractical; TODO28 pursues C68's reverse-instrument lead"
claims: [C56, C64, C65, C66, C67, C68]
---
# One phase interval, not an order-sized sector sweep

**Completed 2026-09-11 at bounded scope.** C66/§UG own the quantitative
certificate and limits; C67/§UG own the exact checkpoint baseline. C68/§RI
record the stronger reverse-instrument proof and tiny-law audit. Current
next work lives in TODO 28. The plan below is historical; no useful
approximation warm-up or general simulation breakthrough is claimed.

This is the bounded follow-up to TODO 26 / §FM. It concerns the repeated
embedded Rx(pi/7) three-state mixer with fixed sector and NO routes. Do not
replace a general known-order output problem by a different noise model.

## Derive the phase reduction before scanning

Write theta=alpha/M in [0,1], T_theta^3=exp(2*pi*i*theta)I. The averaged
FORWARD channels alternate between conjugation by T_theta and T_theta^2,
with the same W after each control. Raw branch phases remain in inverse-QFT
interference and must not be stripped there. A gauge with

    q=exp(2*pi*i*theta/3), S=diag(1,q^-1,q^-2),
    T_theta=q*S*C*S^dag

uses the ordinary three-cycle C and W_theta=S^dag W S. It reduces the gap
question to one compact phase parameter, not a growing arithmetic period.

C66 now proves a qualitative uniform gap, using norm-equality conditions,
the three-state eigendirections and compactness. Main and a lower-cost agent
audited the proof. Read that claim rather than redoing a search for exceptional
sector phases in this exact fixed-mixer family. What remains missing is a
quantitative rate; do not generalize the dimension-three argument blindly.

## Smallest decisive quantitative experiment

Use the existing exact-input builder and density helper. Register predictions
before measurement. Either derive a usable analytic gap or certify a finite
phase mesh with a proved between-mesh error bound. Floating phase samples
alone cannot certify every alpha/M. A Frobenius bound on a sufficiently long
channel product may be simpler than certifying individual singular values.
Changing the word length after a failed threshold must retain that failure
and be explicitly a new experiment, not a silently selected best row.

Keep the no-background conserved-mode control and C65's binary null. Test
the actual quantitative gap on the claimed norm and distinguish a
Hilbert–Schmidt coefficient from a trace-norm coefficient. A discontinuous
modular phase implementation must not be assigned a global derivative bound;
first strip only channel-irrelevant scalar phases using the exact identity.
Charge the interval mesh/certificate construction, not just one chosen sector.

## Cost gate before production changes

Translate the rigorous gap into C66's summed-state-error warm-up. Compare
retained forward scalar entries, setup and total arithmetic against the
unchanged verified finite-work baseline at the same requested full-law TV.
Numerical mass precision, backward effects and QFT branch phase arithmetic
are still required. A bound needing more checkpoints than the current
width<=63 API permits is mathematical progress but NOT an implemented memory
saving. Do not increase the width cap solely to manufacture a size record.

If the generic worst-case warm-up is too conservative, first test a fixed
state-specific summable tail bound or a checkpoint/recomputation baseline,
with all costs charged. Stop if these do not expose a useful regime. Do not
expand to arbitrary routes, large blocks, physical compilation or general
noise simulation before this contract is understood. TODO 24's backend-wide
bit-cost question remains deferred separately.

Consult C66 for primary-source positioning. The general channel mixing and
block-contraction mechanism is known; novelty and broad simulator impact are
not established by a uniform bound for this fixed mixer.
