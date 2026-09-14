---
id: 47
state: done
title: "Can interacting carry chains retain cheap signed contraction after the single-adder baseline stops applying?"
outcome: "C86/MG validate exact planar matchgate contraction of two additions; growing nonvacuous arithmetic queries and practical advantage remain for TODO48"
claims: [C8, C48, C82, C83, C84, C85, C86]
---
# Find an obstruction beyond one carry chain

Completed at bounded mathematical, source-audit and implementation scope.
C86 owns the planar carry theorem and actual dirty-macro baselines; MG owns
the exact experiments, primary-source audit and accepted board evidence.
TODO48 owns the remaining concrete geometric/query discriminator. No broad
simulation breakthrough or measured memory advantage is established.

The original brief follows, with one premise corrected here rather than
silently rewritten: Wallén Section 3.4 concerns parallel output concatenation,
not an efficient general serial-composition representation. MG records the
source audit. The implemented finite comparison fixes two addition passes
and varies their intervening permutation at fixed width, rather than varying
pass count; the selected family and predictions were recorded before its run.

C85/CP close the isolated Cuccaro prefix as a novelty candidate. Do not
repeat its width sweep, rename Wallén's carry automaton, or claim that the
failed unrestricted escape witness succeeds. A genuine conditioned escape
and a simple final character do not themselves supply a new method.

Start with a bounded mathematical and source audit of existing interacting
arithmetic, for example `toffoli_arith.ToffoliModExp.cmult_mod` and its
repeated controlled additions to one accumulator. Read the actual logical
trace and clean/full-space contract before assigning carry variables. An
alternative is two overlapping additions separated by a specified CNOT
map or bit permutation, but explain its arithmetic relevance. A trailing
CNOT layer alone merely changes Walsh query coordinates by C82 and is
not an obstruction.

Fix one exact normalized signed coefficient (or a small predetermined set)
whose observable reaches the arithmetic. For the first finite discriminator
freeze word width and all constants and vary only the number of interacting
addition passes. Derive the factor graph and explicit elimination orders
before simulation. Charge control variables, initial scratch, circuit setup,
observable/query masks and signed extraction. Do not substitute a clean-code
identity for the full-space map of the supplied circuit.

Compare chronological and bit-position elimination with a better available
order, including known arithmetic correlation composition and weighted
automaton minimization. Wallén's Section 3.4 already composes rational
representations; a Cartesian product of carry states is only an upper bound.
A growing boundary in one chosen order does not lower-bound optimized
contraction. Where affordable, exact tiny elimination-width bounds or an
explicit graph-minor argument should separate geometry from a poor order.
Even a graph-width lower bound does not prove memory hardness for a signed
scalar: algebraic cancellation and special factor values can reduce it.

If a constant or otherwise cheap known contraction remains, record that
and redirect without another stabilizer implementation or large pilot. If
the audited comparator actually leaves an opportunity, formulate one extra
prediction: an efficiently checkable signed merge certificate or a bound
on the reduced intermediate information that survives a growing interaction
parameter. Compare against symbolic cancellation, grouped controls and
phase-aware stabilizer coalescing, retaining a sign-cancellation control.
Use independent tiny gate replay to check the same output; do not infer
new memory savings from a dense-vector comparison alone.

Use the authorized board swarm with one canonical writer. Core first for
scientific execution, bounded work under TODO34's host limits, predictions
before runs, and failed scripts/logs preserved before correction. Review
actual harness resolution and exit status: an unconditional success line
is not evidence. TODO42's sampling-cost question remains separate.
