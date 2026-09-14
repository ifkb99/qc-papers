---
id: 30
state: done
title: "Can two-reflection routing replace exponential history enumeration by a short coherent walk?"
outcome: "C71 proves the support and accepted-law identities; audited sparse float sampler passes full laws and improves higher-k implementation timings"
claims: [C56, C58, C59, C63, C68, C70, C71]
---
# Merge coherent histories by reached sector, not by classical dephasing

**Completed 2026-09-11 at mathematical/float scope.** C71 owns the proof and
implementation contract; §CM owns support, full-law, control and matched cost
reports, including main-audit corrections. No finite-TV certificate is assigned
to this changing sparse state. TODO 31 owns the stronger merged-prefix oracle
question; TODO 24 remains deferred. No default or manuscript change.

## Original plan (historical)

TODO 29 / C70 / §QD resolve the bounded-coordinate experiment with a useful
precision theorem but a negative runtime comparison. Do not repeat that
optimization. This candidate instead targets C59's exponential coherent
history count, under a narrower and explicitly recognizable routing promise.

## Candidate mechanism and first algebra audit

Assume every coherent reflection uses one of TWO fixed labels q0,q1:
f_q(alpha)=-alpha-q mod M. Two involutions generate reduced words which
alternate. There are at most 1+2k reduced words of length <=k, so a vector
started at one final gamma can reach at most

    S = min(M,1+2k)

coarse labels after k such gates. This does not require q0-q1 invertible or
an affine change of basis. Arithmetic and within-sector work mixers preserve
coarse labels, so they cannot enlarge this support. Histories arriving at
one label must ADD THEIR COMPLEX FINE-WORK VECTORS, not probabilities. Gates
between the reflections may change those vectors; merging is a dynamical
state update, not combining bare history coefficients in advance.

Main's proposed link is a finite-time walk on a two-generator reflection
graph. Static regrouping can still contain all M sectors when q0=0,q1=1;
that does not require a point-started reverse state to occupy all of them.
This is not a claim that a general coherent orbit or arbitrary q labels
have small support, or a novelty assertion about group theory/quantum walks.

A lower-cost agent independently checked the argument and ran an exact
route-set diagnostic: all ordered binary q words through k=8, with M=7,16
and every gamma, and M=101,1009 with gamma=0,1,17. Closure was updated in
reverse insertion order by R <- R union f_q(R). Every set obeyed S; at
M=1009 the largest supports for k=0..8 were [1,2,4,6,8,10,12,14,16],
against history counts 2^k. The initial test was read-only algebra, not a
stored sampler report. Reproduce it in the harness before promoting a claim.

Its generic-label control was M=1009,gamma=17 and
q=(0,1,0,3,0,9,0,27). Main caught an ORDERING distinction in the initial
report: applying this list directly gives support 88 and contains the 16
negative subset-sum translations; treating it as chronological insertions
and applying it BACKWARD gives support 108 and contains the 16 positive
subset-sum translations. Both fail the proposed two-label bound, but only
the latter is the stated backward test. The exhaustive binary-word maxima
are unchanged because word reversal permutes that entire sweep. More
generally independent pair displacements give exponentially many labels.
This is a control against importing the two-label bound into arbitrary
coherent routing and a warning to audit time ordering in the amplitude test.

## Derive the sampling identity before implementing

Use the exact full-space reverse instrument of C68, with actual coherent
reflection gates in their original time order, rather than selecting one
history. Start from a uniform final (gamma,j). Store only the sparse reached
coarse labels with their b-dimensional fine vectors. The total reverse vector
is normalized when computing conditional branch probabilities.

For the pure initial work state let w0=W0|0>. Its initial coarse amplitudes
are uniformly coherent, not dephased. Candidate terminal acceptance is

    a = |sum_alpha <w0,v_alpha>|^2 / S

for normalized terminal v, or with denominator S*||v||^2 for raw v.
Cauchy-Schwarz gives a<=1 on the promised support. Averaging fresh uniform
gamma,j would then give accepted joint (gamma,y) mass p(gamma,y)/(b*S),
with exact mean b*S attempts. Audit the global normalization from the real
feedback tree and all endpoint gates, not an arbitrary adaptive instrument.

IMPORTANT: redraw gamma as well as j after EVERY rejection. The true
coherent final-sector marginal need not be uniform. C70's fixed-initial-
sector retry rule does not carry over. Use one GLOBAL S bound, not an
uncharged gamma-dependent support denominator which would reweight outputs.

Derive actual operation/storage costs including reached-set construction,
dictionary merging, original gate ordering, phase evaluation and rejected
attempts. A candidate arithmetic improvement is polynomial in k for fixed
b and this two-label promise, while C59 explicitly enumerates 2^k histories.
Do not assume a Python speedup, a lower bound against other tensor methods,
or a certified finite-bit sampler from that candidate arithmetic statement.

## Smallest discriminating experiment

Use lower-cost agents for independent initial tests while main derives and
audits. Start from the existing harness template; keep k within the current
API cap and budget sparse and reference allocations before executing.

1. Reproduce the exact route-set sweep and generic-label must-fail control
   with a retained structured report. Freeze other inputs as k changes.
2. On a tiny b=3 instance with initial W0, several noncommuting work mixers,
   q0/q1 reflections and endpoint insertions, compare EVERY final output
   probability and accepted submass against existing coherent-history
   contraction and an independent full-r reference. Reuse existing branch,
   work-gate and inverse-QFT helpers; no second generic circuit propagator.
3. Negative controls: merge probabilities rather than amplitudes; retain
   gamma through rejection; omit boundary acceptance; or replace the two
   reflection labels with distant labels. Include fixed points, modular
   wrap and collisions of different histories on the same fine-work state.
4. Only after those pass, compare task-matched costs with C59/C63 history
   methods and static regrouping where cheap. Distinguish mathematical/float
   diagnostics from certified same-TV timings. If C70's local-error proof
   extends to changing rectangular sparse supports, derive that extension
   separately before assigning its certificate to a new implementation.

Consult primary sparse-trajectory, finite-time orbit and tensor-network
prior work before any novelty/impact claim. A linear route-set count alone
does not yet establish a simulator. Keep broader q families, large block
dimensions, backend-wide primitive constants and manuscript changes deferred.
