---
id: 24
state: open
title: "Can prefix-amplitude sampling carry a usable numerical error certificate?"
outcome: "implementation complete for the supplied b=2 input: C60–C65 certify prefix, finite-work and scalar proposals; backend-wide primitive constants remain open and deferred"
claims: [C55, C56, C59, C60, C61, C62, C63, C64, C65]
---
# From an exact mathematical sampler to a certified numerical one

**Triage (2026-09-22, METHOD.md barrier check 4):** the samplers here are given the order (a known-order promise), so success computes nothing known to be hard. Passes.

**Progress 2026-09-11 UTC, NOT closed.** C60 proves the full-law block-kernel
budget; C61/§VP implement verified exact-input prefixes and finite-bit sparse
sampling, with tiny complete-law and supplied-wide validation. These do NOT
certify the existing complex128 sampler. §PE retains the preceding precision
audit. Main has read the BGL robustness supplement and the selected backend's
primary documentation; do not repeat that recovery without a new question.

## Current boundary and priority

C64/§FW settle the verified linear-depth proposal question below. The
law-weighted mass bound, zero branches, precision rebuild/replay and complete
accepted-law comparison are implemented and tested. C65 supplies a stronger
normalized scalar proposal for this b=2 promise and shows that the old late
work mixers were invisible to the requested output. The comparison is now
settled at this scope; neither float arithmetic nor arbitrary rounded gates
have thereby been certified.

Keep this item OPEN solely for the backend-wide primitive-error/bit-runtime
boundary. Defer constant polishing unless it changes a concrete accuracy or
cost decision; finite maxima cannot supply a uniform proof. TODO 25 is the
next experimental priority, testing certification beyond the binary scalar
tail. Do not repeat binary timing or a larger known-order size record.

## Completed question: can certification retain a linear-depth proposal?

**C63/§RA now supply the fair certified rejection comparator.** Both proposal
and acceptance errors, finite bits, retries and refinement are charged. Its
proposal uses C61's prefix sampler restricted to one routed history, not a
verified version of C59's O(t*b^3) finite-work instrument. Do not repeat the
same fixture timing to rediscover that rejection can win despite retries.

The completed bounded experiment was to determine whether the stronger finite-work
instrument can retain that linear-depth contraction at the SAME exact-input
and TV promise. Reuse `sequential_path`'s forward-state/backward-effect
identities, not a new generic propagator. Derive an absolute, law-weighted
error bound for its unnormalized branch weights BEFORE coding; conditional
normalization by tiny mass is not itself a usable uniform certificate.
Either demonstrate a certified finite-work proposal on the same complete
tiny law and supplied-wide input, or identify a precise obstruction to this
implementation route. Compare setup, effects, precision and total bit cost,
not just floating matrix-operation counts. Keep C63 as the verified baseline.

The retained implementation plan was: use lower-cost initial proof/edge testers while main audits. Keep the fixed
gates and output contract; include exact-zero/near-zero branches and a
normalization-drift control. General exact-real sampling, new gate families,
large-k scaling and manuscript expansion are deferred until this narrower
comparison is settled or deliberately delimited.

C62's backend-wide primitive constant remains a separate open issue. A
conditional primitive-error model is proved, not an unconditional bit-runtime
theorem. Resolve it from explicit error contracts/implementation analysis or
state the precise remaining boundary; do not substitute a larger sampled
maximum. The current certificate itself does not depend on guessing a
working mantissa from that maximum.

## Completed comparison plan (retained context)

C62/§NG now separate rectangular wrapping, adaptive doubling overshoot and
actual midpoint error, prove norm-aware residual transport, and give a
precision bound under explicit uniform primitive-error assumptions. The
opt-in norm implementation passes its tests but loses the measured timing
comparison, so do not change the default or repeat a sweep to rediscover
wrapping. A backend-wide numerical value for the primitive-error constant
is still unproved; this is different from the already proved error certificate.

The algorithm-level comparison below is now realized by C63; retain the
backend-constant caveat above.

Implement or decisively delimit a certified rejection comparator under
the SAME exact-input and total-variation promise. Charge proposal, acceptance,
finite random bits, failed attempts and numerical refinement together using
C60's accepted-measure analysis. Do not benchmark the verified sampler against
uncertified float rejection as though their accuracy contracts matched.
The oracle/kernel implementation portion is complete; a uniform bit bound or
a bounded, explicit obstruction plus an honest comparison limitation can close
this first pass. Keep general exact simulation outside this narrower result.

For lower-cost initial tests, freeze the same exact r10,t4 gate fixture.
Separate the proposal distribution from the acceptance comparison. In
particular, the existing floating finite-work instrument does not become
certified by evaluating its final acceptance probability with Arb. Derive the
accepted-measure error or a verified adaptive comparison FIRST, including
zero-proposal/near-zero controls and termination semantics. Cross-multiplying
acceptance inequalities may avoid dividing by tiny proposal mass, but prove
what it buys before implementing or calling the output law exact. Consult
primary prior art if using variable-length random-bit/exact-real sampling.

## Completed implementation boundary (retained plan)

The following implementation plan is now realized in C61/§VP.

Build a verified absolute-error version of the EXISTING finite prefix
contraction, not another generic propagator. Start with b=2 background Rx/Rz
and reflection angles specified as exact rational multiples of pi. This is an
explicit exact-input specialization: reinterpreting arbitrary rounded matrices
as exact unitaries is invalid. Inspect an established ball/interval arithmetic
backend and its primary documentation before using it.

An adaptive evaluator can refine its arithmetic until every returned scaled
coordinate has radius below the planner's fixed tolerance. Its stopping rule
must work for any label, including exact cancellation, without dividing by
small amplitudes. Explain why refinement terminates and charge its bit cost;
a resource cap that can fail on some labels does not automatically certify
the distribution conditioned on successful runs. Keep output midpoint
representation and its rounding error inside the same guarantee.

Then implement a categorical kernel with exact rational/dyadic weights or
verified CDF boundaries, unbiased initial integer draws and explicit zero-block
semantics. Charge finite random bits separately. A certificate on a single
successful sampled trajectory does not bound the overall output law unless
the oracle/kernel promises are uniform. Compare a tiny complete transition
law and one supplied-wide execution, retaining the original production sampler
as an unchanged numerical comparator. Do not claim certified production until
both the oracle and sampling-kernel assumptions are actually satisfied.

The original first-pass plan below is retained for context; its oracle-level
derivation and tiny diagnostic portion are completed in C60/§PE.

C59 now supplies prefix amplitudes and two specialized sampling algorithms.
Its complex128 implementation is tested, not numerically certified. Before
extending to another artificial gate family or increasing the size record,
ask whether the sampled law can be bounded at a requested accuracy without
enumerating exponentially many sectors, strings or output probabilities.

## Read and derive before measuring

Read the robustness argument, including its supplement, in
[Bravyi–Gosset–Liu, arXiv:2112.08499](https://arxiv.org/abs/2112.08499).
At the start of TODO 24, main had read Algorithm 2's proof and Lemma 1's
statement but not the supplement; that reading is now complete. It assumes global state approximation;
a small error in a few queried amplitudes does not establish that hypothesis.
Identify the exact input promise and error metric before importing a bound.

For the known block partitions in C59, derive a distribution-level update
bound that weights rare blocks by their actual mass. Dividing a pointwise
amplitude error by an arbitrarily small block norm is not a uniform sampling
certificate. Consider whether a sum of normalized-component contraction error
bounds, charged by the coefficient norm, can establish the needed global norm
promise without enumerating basis states. Include nonunitary floating-point
drift and the consistency needed between successive prefix queries.

Separate three tasks: a theorem under an explicit amplitude/state-error oracle,
a computable conservative error budget for THIS contraction, and a measured
precision diagnostic. A solution to the first is not automatically the second.
Charge mantissa precision, integer reduction, trigonometric evaluation, RNG
resolution and any failed/retried numerical branch. No biased fallback or
unrecorded amplitude deletion.

## Smallest discriminating test

Delegate bounded initial tests to lower-cost agents while main audits the
error argument. Keep the r=10,b=2,t=4 coherent fixture and gate schedule fixed.
Derive a cancellation/rare-block control first, then vary only precision or
one strength at a time. Enumerate the full tiny classical transition law,
not just sampled histograms, against the existing independent full-r
reference. Higher precision is an independent numerical check, not itself
an interval certificate; use existing contractions, not a new propagator.

Must-fail controls should expose at least one false inference: ignoring
normalization drift, certifying from pointwise error alone, or discarding
small amplitude bins without charging discarded probability. Every compared
TV law must normalize. Also test exact-zero and near-zero blocks separately.
Preflight reference arrays and runtime before any sweep. Retain every failed
prediction and distinguish theorem violations from coding/precision failures.

Compare to the rejection sampler under the SAME accuracy promise. It may
offer a simpler certified envelope despite its extra proposals. If no usable
certificate follows without an exponential global calculation, record the
specific gap and a bounded counterexample rather than calling the sampler
certified or declaring a general impossibility. A useful conservative bound
or decisive limitation is enough to close this first pass.
