---
id: 29
state: done
title: "Can local instrument error certify bounded-bit reverse trajectories?"
outcome: "C70 proves and implements the local-error bound with bounded integer work states; full finite laws pass, but matched timings lose to existing methods"
claims: [C56, C60, C61, C63, C64, C65, C67, C68, C69, C70]
---
# Local errors in a joint output/work process

**Completed 2026-09-11.** C70 owns the proof and implementation; §QD owns
the main-audited full-law, edge and serial cost evidence, including failures.
The depth-local bound succeeds; the runtime comparison is negative. Stop
polishing this representation at this scope. No default or manuscript change.
TODO 30 owns the next structural coherent-route investigation; TODO 24 still
owns the deferred backend-wide bit-runtime question.

## Original plan (historical)

C69/§RV complete the unnormalized streaming baseline. It trades memory for
rejection work and did not win the five-seed wide timing comparison. Avoid
another size sweep or a default-method change. Test whether a different
numerical representation can remove its independent absolute-mass tree factor.
This is NOT already an implementation or a theorem asserted by C69.

Initial independent lower-cost algebra audit completed at the end of TODO 28:
it supports the stacked-instrument/hybrid argument provided the classical
prefix is retained and canonical rounding preserves a maximal ±1 coordinate.
There is no integer prototype or measured finite law yet. Treat the proof
sketch below as the starting point for main's full derivation and tests.

## Derive and independently audit first

View the reverse sampler as a classical output record plus a conditional
normalized pure work state. Its exact feedback-controlled instrument is a
linear trace-preserving map on this joint classical/quantum object, despite
individual normalized branch updates being nonlinear. This may allow an
aggregate error bound per depth instead of counting every prefix.

Candidate construction: round the shared exact reverse Kraus operators A_z
to dyadic matrices L_z; compute their child vectors and norm weights with
exact Gaussian integers. For a nonzero selected child w, let r be the maximum
absolute real/imaginary coordinate, divide by r using exact rational arithmetic,
then round to a fixed dyadic grid. A maximal coordinate stays exactly ±1.
Thus the represented normalized state has a bounded-scale representative
without a square root or high relative precision on a tiny physical norm.
Keep the branch weights from BEFORE this state compression. Never replace
their physical Born weights by the compressed norm weights after the fact.

Main's candidate bound to audit: if stacked exact A is an isometry and
||stack(L)-stack(A)||<=nu<1, unnormalized joint-output/work trace error is
at most d=(2+nu)nu. Normalizing the approximate instrument's TOTAL branch
mass adds at most another d. Dyadic child compression with normalized pure
state trace error at most kappa, and a two-bin L-bit kernel, then contributes
at most 2d+kappa+2*2^-L per depth. The nonlinear approximate update must be
shown uniformly close to the EXACT linear instrument on every normalized
input; do not assume nonlinear state normalization itself is contractive.

For entrywise eta-grid rounding, conservative candidates are nu<=2b*eta
and kappa<=min(2,4b*eta), with eta=2^-p and nu<1. Audit these including negative/tied largest coordinates,
exact-zero children, normalized density errors, and positivity of total mass.
Terminal acceptance uses exact rational overlaps of normalized stored integer
vectors; separately charge the initial physical state's approximation and
downward acceptance-bit error. C60 then amplifies accepted-measure error by b.
If justified, requested bits may scale as log(b*t/delta) rather than t plus
log(1/delta); that is not yet a sufficient backend/runtime statement.

## Smallest decisive tests

Use lower-cost agents for an independent algebra audit and an initial exact
integer prototype; main audits before any claim promotion. Reuse existing
branch/QFT constructors, not a new generic propagator. Register predictions
before measuring and first enumerate complete tiny finite laws against C69's
target references. Include a rare branch, zero child, tied/negative scale,
initial W0, nontrivial routes and the frozen-j/reset/omitted-acceptance controls.
Show that decreasing eta contracts the PROVED full-law bound, not just a
selected amplitude error. Count stored integer bit lengths and live buffers.

Only after this passes, compare at identical requested TV against C69, full
storage, C67 checkpoints and C65 where applicable. Charge verified operator
construction, exact integer products/divisions, compression and all rejection
attempts. A smaller requested coordinate tolerance or word count is not a
speedup by itself. Keep width<=63 and the existing supplied-input promise.

Consult task-matched primary quantum-trajectory/approximate-instrument prior
work before novelty claims. If the local bound or exact integer overhead
fails to improve the comparison, retain that negative and stop optimizing
this representation. Coherent outer integration, arbitrary block dimensions
and backend-wide primitive complexity remain separate from this bounded test.
