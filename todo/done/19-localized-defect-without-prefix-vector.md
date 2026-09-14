---
id: 19
state: done
title: "Can a localized defect remove the early-prefix vector too?"
outcome: "finite-support geometric sums plus scalar progression marginals yield an actual sampler; averaged rejection removes the prefix array under the known-index contract"
claims: [C54, C55]
---
# Localized work defect without a length-2^s prefix

Completed 2026-09-11 UTC at the stated single-finite-support-defect scope.
C55 owns the proof and limits; §LF records tests, lower-cost agent checks,
preserved verifier failures and the impact audit. The historical predictions
below preceded the experiments. TODO 20 owns the next broader-support test;
TODO 14 and the user's continuing research goal remain open.

C54 removes an orbit table under an explicit known-index contract, but its
early prefix still costs L=2^s amplitudes. Start from the already verified
two-level mixing defect; do not search for a larger headline instance.

## Candidate algebra to verify before expanding the harness

Let V equal identity outside a fixed set S of d known orbit indices. The
unnormalized conditional early row should be

    h_k(l) = exp(2*pi*i*k*l/r)
               + sum_{p in S} delta_p(k) 1_{l mod r = p},
    delta_p(k) = sum_{q in S} (V[q,p]-1_{q=p}) exp(2*pi*i*k*q/r).

Its squared norm can be computed using residue counts n_p in [0,L):

    ||h_k||^2 = L + sum_{p in S} n_p (|exp(2*pi*i*k*p/r)+delta_p(k)|^2-1).

Each output Fourier amplitude then splits into one plane-wave geometric sum
and d geometric sums over truncated arithmetic progressions. Test these
candidate formulas against the existing explicit prefix implementation before
claiming that its exponential storage is avoidable. Include non-power-of-two
periods, L<r and L>r, both angle signs, zero angle, and complex off-diagonals.
Initial tests should again go to a lower-cost model under the user's request.

One-parameter series: fix the circuit/orbit, total exponent width and defect
angle, then vary the insertion index s by one. Independently verify the actual
small coherent gate circuit at representative insertion points. Keep the
known-order and known-index assumptions visible and charge their setup.

## Do not confuse a finite formula with a sampler

A possible proposal distribution is a mixture of the Fourier intensities of
the plane wave and the d progression terms; Cauchy-Schwarz supplies an envelope
for their coherent sum. Before implementation, derive its NORMALIZATION and
the acceptance cost averaged over sampled final eigenphases. Very small
conditional norms can make a superficially good pointwise envelope useless.

The plane-wave term has an existing scalar phase sampler. The truncated
progression terms need their own correctly normalized sampling method; a
closed formula for one probability is not sufficient. Search primary prior art
for that component. If it requires an orbit/prefix-sized table, record the
limitation and do not advertise a full table-free sampler. Must-fail controls
should detect deleted interference, incorrect progression lengths and missing
inverse-QFT feedback.

This is a candidate structured-defect extension, NOT a general Shor simulator
or an established breakthrough. Reassess impact after those discriminating
tests rather than committing to a general tensor framework.
