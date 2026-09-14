---
id: 28
state: done
title: "Can reverse-instrument rejection beat checkpointing at certified accuracy?"
outcome: "C69 certifies an unnormalized reverse-vector sampler; complete finite laws and edge controls pass, working storage shrinks but matched five-seed timing loses to checkpointing; TODO29 owns local-error compression"
claims: [C56, C59, C60, C61, C63, C64, C65, C66, C67, C68, C69]
---
# Certify the reverse trajectory before scaling it

**Completed 2026-09-11 at bounded scope.** C69/§RV own the implementation,
proof, matched-accuracy comparison and preserved failures. The simpler
unnormalized-vector route avoided rare normalization; it does not establish
a stronger normalized-instrument precision bound. The plan below is historical.
Current next work lives in TODO 29. No default sampler or manuscript changed.

Read C68 for the proof, §RI for tiny-law evidence/failures, and C67/§UG for
the implemented storage comparator. Do not re-scan exceptional sector phases
or optimize the impractical generic warm-up before testing this stronger
construction. It needs neither mixing nor forward checkpoints.

## First bounded experiment

Implement an opt-in pure-trajectory or normalized-density reverse component
using existing branch constructors and the QFT helper, not a new generic
propagator. Preserve raw phases, original time ordering, initial W0 and
deterministic routes. Keep b=2/3 and width<=63; charge supplied order/index,
gate/route data and reconstruction. Order/orbit discovery is not solved.

Derive the FULL accepted-law finite-precision contract BEFORE calling it
verified. One possible route is trace-distance contraction on the joint
classical-output/quantum-work process, including feedback on recorded bits.
This might avoid the full-tree coefficient of unrelated unnormalized mass
errors; it is a hypothesis to prove, not an existing bound. A simpler
rigorous initial bound is acceptable. Charge terminal acceptance and
normalization of its accepted submeasure too.

Pure-state normalization by sqrt(w_z) can be ill-conditioned at rare prefixes.
Use honest absolute-error/accepted-mass analysis, adaptive exact-input
enclosures and zero-branch handling. Do not assume a minimum positive mass,
discard rare paths or call float clipping a certificate. Precision refinement
must preserve selected bits and boundary j, with replay/reconstruction charged.
Attempt caps or resource errors cannot be conditioned away for free.

Register predictions first; fix §RI's tiny circuit and vary target accuracy.
Enumerate its finite-bit accepted law across histories, sectors and leaves,
not a histogram. Compare target intervals and the independent full-r reference.
Include bounded zero-width, zero/rare-mass, nonzero-W0, endpoint-route and
wrong-adjoint controls. Keep the rational nonnormal regression. Omitting
acceptance must fail; mixed I/b retains its constant-acceptance null. C65
remains the stronger baseline for its binary scalar family.

Only then compare SAME-TV setup/sampling costs against full storage and C67:
attempts, precision, replay, state/branch entries, native-memory limits and
input/output bits. Reuse the supplied wide §UG fixture without orbit/output
tables. Charge retries and retain failed comparisons; do not select a
favorable seed or call scalar-entry counts RSS.

## Scope and stopping rule

Direct deterministic components come first. Coherent histories require a
separate composition with C63, not reuse of this acceptance rule without proof.
General dimensions, nonunitary channels, physical compilation and more size
records are deferred. Known reversal/trajectory methods and checkpointing
are baselines, not new principles; inspect further task-matched primary prior
work before novelty claims. If numerical or rejection costs erase savings,
record the matched-accuracy negative and retain checkpointing. TODO 24 keeps
the separate backend-wide bit-complexity issue.
