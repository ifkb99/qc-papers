---
id: 72
state: open
title: "Can the compact divisor state be constructed from the composite alone?"
outcome: "First bounded round complete (2026-09-21): C118/C119; constructive aggregation remains open"
claims: [C117, C118, C119]
---
# 72: Can the compact divisor state be constructed from the composite alone?

User-authorized continuation toward the public RSA Factoring Challenge numbers.
Start from C117's constructive structured sampling lesson; do not assume that
a low-rank final factor state supplies its own construction.

## Contract

Input: an odd composite N and declared register bounds, with no factors, order,
factor-dependent advice or precomputed solutions supplied to the candidate.
Output: a nontrivial factor, or exact counts/samples of the bounded relation
p*q = N with all preparation charged. Separate proof, exploration and measured
performance. Public solved RSA instances are eventual benchmarks, not evidence
of a new factoring method merely because their known factors fit a compact state.

## Completed decision

Note FC records the independent slates and selected derivation. C118 settles
the specified explicit low-bit residual constructor by proof; C119 states the
construction/access criteria and modular-state obstruction. TX39-TX43 retain
the competing directions and their hypotheses. No benchmark or RSA challenge
factorization was needed to settle that constructor.

## Next deciding step

Preserved for later by arb decision M8c90bbdf49a94689 (topic:rsa72). The wider
reviewed brainstorm is note FD/TODO73 and TX44-TX47; no constructor experiment
was restarted during that round. The original question below remains open.

Derive a specific N-only algorithm to aggregate many modular-inverse residues
inside dyadic interval rectangles, or select a small useful family of high-bit
boxes before low-bit expansion. Charge how the boxes/counts are constructed,
integer precision, all rejected work and extraction. Supplying leading factor
bits or enumerating all boxes is not an escape. Alternative: derive a genuinely
stronger zero-set simplification for the compact word-level product (TX41).

If a constructive recurrence closes with a proved manageable number of states,
freeze its exact output and design a bounded validation against direct arithmetic
and an appropriate classical factoring baseline, including discriminating
controls. If it reduces to unknown completion counts or expands to individual
residues, record that scoped failure and redirect. Do not run another sweep
of C118's already-settled constructor. RSA-100 is a future public reproducibility
target only after the small-instance and cost gates, not a scheduled attempt.

## Budget and completion

Local CPU only; no paid compute or challenge-scale brute force. Pure derivation
first; any later diagnostic receives its own explicit cap and reviewed design.
Record the chosen mechanism, strongest applicable obstruction, baseline and
next step under each outcome; integrate reviewed evidence and preserve all
existing dirty work.
