---
id: 21
state: done
title: "One localized symmetry-breaking kick amid periodic mixers"
outcome: "bounded probe complete: coarse dephasing fails, a low-rank complete-output formula still sums over sectors, and a co-moving support cover certifies omission; exact compressed sampling remains open"
claims: [C53, C55, C56, C57]
---
# Break the actual sufficient condition, not merely commutation with U

Completed 2026-09-11 UTC at bounded investigation scope, not as an exact
compressed-sampling solution or hardness proof. C57 owns the mathematical
results and unresolved marginal step; §LC records independent initial tests,
main audit, failed controls and the certified-approximation alternative.
The pre-experiment question below is retained. TODO 14 keeps the broader
exact frontier open; TODO 22 is the next physically motivated diagnostic.

C56 explains why multiple full-orbit periodic mixers remain cheap: they share
U^b's conserved sectors. Another noncommuting periodic mixer does not challenge
that explanation. The next candidate is ONE localized orbit-basis kick K amid
otherwise periodic mixers, with [K,U^b] nonzero. This combines C55's finite
support with C56's tractable background. Neither prior sampler directly covers
the combination. Failing their assumptions is not itself a hardness result.

## Derive and test the smallest observable distinction

Keep N=7,a=3, the known six-label orbit, b=3, and the tested W01/W12 physical
construction. Fix exponent width, all insertions, background blocks and kick
support. Vary only the kick angle, including zero. A rotation confined to a
single pair of orbit labels in ONE block breaks the repeated-block promise;
construct it with the existing Circuit Pauli/projector gates, checking every
work input and invalid labels before interpreting the outputs. Choose the
insertion after some background mixing and before later arithmetic, avoiding
an end-only work gate or an immediately cancellable early-span example.

Prediction: the correct r-dimensional original-order `sequential_path`
baseline agrees with the actual small Circuit/statevec output for every tested
angle. A model that measures the coarse sector at the start generally loses
observable interference once K couples sectors. This must be measured at the
output-probability level, not just via nonzero commutators or amplitude changes.
At zero kick, the existing periodic sampler must recover the reference.
Restoring the same kick in EVERY block must restore the sector promise.
An intentionally incorrect initial coarse dephasing should fail on at least
one fixed nonzero angle; if not, derive the next input-specific shortcut before
changing parameters. Do not reorder the arithmetic/defect stages.

Send bounded initial tests to lower-cost agents under the user's standing
request while main works on the representation/proof. Use the experiment
template and capped existing propagators; do not build another generic state
simulator. Keep r<=18, t<=7 and a preflight dense payload limit for initial
diagnostics. Charge construction and retained forward/backward matrices.

## Ask whether the violation is cheap to correct

In the coarse basis, K-I has matrix rank at most its support size d but can
couple all sector pairs. Low rank alone is NOT a sampleability or memory bound.
Derive the off-diagonal sector contribution to the actual output law before
choosing tensors or Monte Carlo. Can a finite-support insertion be handled by
a few boundary amplitudes, a normalized positive mixture, or controlled
interference rejection? Charge any sum over all M sectors, preprocessing,
cancellation condition number, proposal normalization and sampling cost.
An efficient single-probability formula is not yet an efficient sampler.

Also compare the unperturbed periodic sampler as an approximation. If F is
the ACTUAL pre-kick probability in the affected support, the same unitary
distance argument as C55 gives TV<=min(1,2*sqrt(F)). After earlier periodic
mixers, F is not automatically C55's uniform residue-count formula. Compute
it from the correct pre-kick state on the tiny instance. If this bound or a
sharper angle-dependent bound already covers the effect of interest, say so.

Outcomes: a useful positive result needs a genuine sampler/cost argument beyond
the common symmetry; a negative result should identify the failed shortcut
without claiming exponential necessity. If another small exact shortcut wins,
record its promise and redirect to precision certification or a genuinely
different structural family, not more large-known-period benchmarks. Do not
promote this perturbed-circuit branch into either manuscript before assessing
its task-matched value and primary-literature position.
