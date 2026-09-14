---
id: 20
state: done
title: "Can a periodic work mixer retain strong output effects without large sampling cost?"
outcome: "a conserved coarse Fourier sector gives a rejection-free multi-defect sampler; single-defect harmonic formulas and two distinct cancellation mechanisms are verified"
claims: [C52, C54, C55, C56]
---
# Beyond a defect that hits only a few basis labels

Completed 2026-09-11 UTC at the supplied periodic-block/known-order scope.
C56 owns the proof and restrictions; §PS records initial lower-cost tests,
main implementation/audit, physical circuit validation and failed controls.
The original candidate below is retained as the pre-experiment prediction.
The stronger coarse-sector baseline covers several separated periodic defects;
TODO 21 asks what survives one defect that breaks this shared symmetry.
TODO 14 and the user's continuing research goal remain open.

C55 completes the finite-support sampler and bounds when the simple ideal
sampler is already a good approximation. The bounded impact probe in §LF has
resonant zeros and nonmonotone output TV. Do not infer a useful large-instance
advantage merely from a large supplied r or a small matrix-payload number.

## Next candidate: a periodic block mixer over the whole orbit

Use a known r divisible by a fixed b, and apply the SAME b-by-b unitary W
to each group of orbit indices bm,...,bm+b-1. This affects the whole orbit
but each column has at most b terms. It remains an abstract known-index
construction, not automatically a constant-locality gate on physical modular
work qubits. Keep that distinction explicit.

Candidate derivation to test before implementing a new helper:

    A_p(k) = sum_q W[q,p] exp(2*pi*i*k*q/r),
    f_p(k) = exp(-2*pi*i*k*p/r) A_p(k),
    h_k(l) = exp(2*pi*i*k*l/r) f_{l mod b}(k).

The last equality uses b|r; it is not valid for arbitrary incomplete orbit
blocks. Expand the period-b multiplier in b Fourier modes,

    f_p(k) = sum_{u=0}^{b-1} c_u(k) exp(2*pi*i*u*p/b).

Then h_k is a sum of b plane waves, each with an existing scalar Fourier
sampler. Unitarity suggests sum_u |c_u(k)|^2=1. The component envelope would
have T_k=L, so its rejection mean averaged over w_k would be <=b; phase
proposals cost b using C54. Verify the signs, norm, zero phases and actual
sampler normalization rather than merely measuring this hoped-for bound.

Use b=3 before a binary block. There is a sharper 2-adic control to verify:
the repeated block mixer commutes with orbit translation U^b. If b divides
L=2^s, it commutes with ALL later controlled powers and can move to the end,
where a work-only unitary cannot change exponent outputs. Thus power-of-two
b should become exactly invisible after the corresponding split threshold;
an odd b>1 never divides L. This is a candidate bridge back to Paper B's
2-adic structure, not merely a larger version of the same size benchmark.
For the preceding endpoint impact probe, r|L likewise makes every later
arithmetic power identity, explaining why its binary resonances warrant
special caution. Test these mechanisms explicitly before promoting them.

Start
with r=6 or 9, t fixed, one seeded complex W held fixed, insertion s varying
by one. Compare every output with existing prefix and original-time-order
spectral references. Only after that compare a fixed b=3, t,s,W series over
small r divisible by 3, recording TV from the ideal sampler. Do not confound
size with a newly drawn random block. Controls: identity W, s=t, wrong phase
weights, deleted harmonic interference, and b not dividing r (should invalidate
the unmodified formula, not silently change its promise).

Send the initial tests to lower-cost agents under the user's standing request;
main should derive/audit and extend the harness only if this differs from a
known simpler explanation. A sum of b product phase states already has a
bounded-bond MPS description, so Browne's terminating-QFT framework in C55 is
relevant prior art. The potentially useful angle is the combined sampled final
phase and few-harmonic representation, not a new claim that low-bond states
can be simulated. Preserve order/index discovery and numerical-precision costs.

If the effect remains measurable while storage depends on b rather than r or
2^s, identify which broader structural promise this actually supports. If it
again becomes ideal or transfers to a tiny exponent subsystem, record the
shortcut and reconsider multiple separated defects or precision certification
instead of running a larger headline benchmark. Neither a third paper nor
generality/novelty is authorized by a finite passing row.
