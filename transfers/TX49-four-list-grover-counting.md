---
id: TX49
field: "Time-space algorithms for subset sum and exact counting"
status: open
effect: unknown
one_line: "Subset-sum Grover sampling is prefix counting plus a solved coin (C122); open question: multiplicity- and prefix-preserving four-list counting in quarter-list space"
source: "Belova et al. arXiv:2402.13170v1 section 1.1.2; memory74 surveyor S1"
claims: [C117, C119, C122]
notes: [FG, GR]
todo: [74]
---
# TX49 — Stream a counting interface for a nonlocal Grover predicate

## Dictionary

Use the explicit oracle f(x)=[sum_i w_i x_i=T]. Split variables into four
quarters and stream sorted sums of pairs of quarter lists through priority
queues. Match left and right streams, retaining multiplicity counters rather
than every matching pair. Conditional marked counts and complementary
unmarked counts would support sequential uniform sampling within either class.

## Hypotheses

[Belova et al.](https://arxiv.org/html/2402.13170v1), section 1.1.2, explains
the Schroeppel-Shamir pair-sum mechanism and its time/space tradeoff.
The original report was not body-audited in this round. The source's stronger
Monte Carlo decision algorithm (section 1.2.3, Theorem 1.4) is not imported as
an exact counter or unbiased sampler. This is a known-method transfer proposal.

The proposed contract is standard Grover from uniform input followed directly
by computational-basis measurement, K outputs with joint TV error epsilon.
It omits C117's final Hadamards: counting ordinary solutions does not give
coherent Fourier branch masses. Exact counts feed the scalar class-selection
probability; only that scalar needs certified approximation.

Pending a multiplicity/prefix proof, target per-count costs are
O*(2^(n/2)) time and O*(2^(n/4)) space instead of half-list storage
O*(2^(n/2)); stars hide polynomial factors in n and input weight bits.
Charge setup, up to 1+Kn count calls, scalar precision and Kn output bits.
Exact rational random draws have expected, not bounded worst-case, time.
Repeated, negative and zero weights, ties, arbitrary prefixes and empty/full
classes must be handled without retaining all matching assignments.
Compare exact half-list counting and small-range dynamic programming on the
same sampling contract. Matching time exponents alone does not establish a
moderate measured slowdown.

## Consequence for the goal

Potential exponential-factor memory reduction in a broader structured Grover
family, while both time and space remain exponential. It would not show a
classical advantage over specialized subset-sum solving. Modular multiplication
has no supplied additive monotone decomposition; discrete logarithms cannot
be introduced for free. C119's factor-count barrier is unchanged. TX48's
precision tools could be shared, but its post-H output sampler cannot simply
be copied. TODO74 owns the deciding counted-stream derivation.

## Status after C122 (2026-09-22)

C122 derives the reduction in this row's dictionary: with the coin (given H_elem
at binary t; exact at polynomial t), Contract-A sampling needs only prefix
counts. The remaining hypothesis concerns a classical counting algorithm and has
not been checked. It stays open. TODO74 deprioritizes it, because it would not
remove an exponential from simulation.
