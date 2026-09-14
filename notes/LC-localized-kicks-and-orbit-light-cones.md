---
code: LC
date: 2026-09-10
title: "One symmetry-breaking kick: visible sector coherence, a finite boundary formula, and an orbit-displacement omission certificate"
outcome: mixed
claims: [C57]
todo: [14, 21, 22]
---
# LC — Exact symmetry fails, but its approximation survives in a bounded regime

TODO 21 challenged C56 with a localized kick inside its periodic background.
The first test confirms that the shared-sector shortcut genuinely fails at
the measured-output level. Main then derived two different follow-ups: a
low-rank complete-output boundary formula, and a co-moving support cover that
can certify omission of the kick without evaluating its exact hit probability.
C57 owns the mathematics and limits; exact compressed sampling of the
symmetry-broken law remains unresolved.

The user's standing request for lower-cost initial testers was followed with
three `gpt-5.6-luna` agents. They handled the physical kick, actual-state
distance bounds, and independent support-cone tests. The physical tester also
checked the derived boundary formula. Main derived/implemented the integer
cover helper, reviewed the scripts, repaired controls/resource claims, reran
everything and added regressions. No existing propagator was replaced or
extended. The research skill's strongest-baseline rule made the omission bound
a required comparator; its audit rules caught the misleading initial claim
that the boundary verifier used constant sector storage.

## Fixed physical experiment

Use N=7,a=3, orbit [1,3,2,6,4,5], b=3,t=5. Background W01@s1,W12@s3,W01@s4
is held fixed, with each block an Rx(pi/2). Insert an Rx(theta) mixing ONLY
orbit labels 0 and 1 (physical labels 1 and 3) after s=2. It is a rotation on
work bit 1 conditioned on work bit 0=1 and work bit 2=0. Four commuting Pauli
rotations implement that projector using the existing Circuit engine.
All work inputs are checked, including unchanged invalid labels 0 and 7.

The one-parameter sweep is theta=0,pi/8,pi/4,pi/2,pi. The actual original-order
compiled Circuit/statevec output agrees with full-r `sequential_path` to
below 2.3e-14 throughout. The physical norm error stays below 2e-13.
This full-r baseline is NOT a scalable small-sector sampler.

Wrong initial coarse dephasing means changing only the initial orbit density
from |0><0| to (|0><0|+|3><3|)/2, then applying the SAME full-r branches.
It is valid at zero kick, but its output TV rises to about 0.160 at theta=pi.
Repeating the kick in both blocks restores C56's symmetry and its correct
coarse mixture. Treating the one-block kick as that repeated gate instead
changes outputs by TV about 0.164. These are changes of the observable law,
not merely nonzero matrix commutators.

The separate bounds agent constructed the actual pre-kick Circuit/statevec
state. Its support probability is F=1/2 within roundoff and its X_S coherence
is numerically zero. F happens to equal the bare count in THIS fixture; this
does not validate that shortcut after arbitrary mixers. The independent
support-cone control below disproves the general shortcut.

| theta | TV from unperturbed periodic sampler | vector-distance upper bound | pure-state overlap upper bound |
|---|---:|---:|---:|
| pi/8 | 0.0372 | 0.1386 | 0.1383 |
| pi/4 | 0.0834 | 0.2759 | 0.2733 |
| pi/2 | 0.1989 | 0.5412 | 0.5210 |
| pi | 0.4578 | 1 | 0.8660 |

The exact zero-angle distance is zero; the computed output TV there is about
8.5e-14 from reference roundoff. The overlap and squared-distance identities
agree with direct states within 7.4e-14 and 2.3e-16, respectively. These bounds
are informative for small angles, but loose. The generic support cover is
trivial on this tiny fixture because its expanded neighborhood covers all
six orbit labels.

## A support-only certificate that does not need the pre-kick state

The light-cone tester held one dense three-point Fourier block fixed across
r=6,9,12,15,18, s=1,...,6 and bounded mixer-count/support combinations. Each
conditional early branch was evaluated by small finite matrix products, not
by another full joint-state propagator. In 342 cases, no support appeared
outside the predicted radius; brute-residue and quotient/remainder counts
agreed exactly; actual F never exceeded the expanded cover. Norm checks are
included in the main-reviewed run. This covers L<r and L>r, zero mixers,
wrapped supports and full coverage. There is NO L=r row in this b=3 series:
a power of two cannot equal one of these periods. An initial agent summary
overstated that coverage; the raw rows are authoritative.

The decisive control uses r=6,s=2,m=1,S={4,5,0}: actual F=1/3, whereas the
unexpanded bare count is 1/4. Thus earlier periodic mixing matters even though
the early exponent labels remain uniform. This is the reason for expanding
the support rather than importing C55's old hit-mass formula unchanged.

Main's `lightcone_cover` helper merges circular intervals and returns exact
integer fractions, including the mathematical squared-TV omission bound.
Its separate counter audit checks 3,912 finite cases against brute force,
including empty support, overlaps, full coverage, wrapping and truncated
traversals. `PeriodicOrbitCircuit.localized_kick_bound` derives the radius
from preceding background blocks, with an explicit same-insertion order.

A wide execution check uses supplied r=3,000,000,021,t=63 and b=3 blocks at
0,16,48. For a kick on S={0,1} after s=32, there are two preceding mixers.
The exact cover has 16 eligible early labels, represented by two intervals;
the squared-TV bound is 64/2^32, hence the mathematical TV bound is 2^-13.
A seeded background draw was also executed. This is an omission certificate
for an EXACT background sampler, not an exact sample of the perturbed law and
not certification of the complex128 implementation's numerical error. The
large r is supplied; no orbit search, indexing or factoring was performed.

## Low rank helps a complete probability, not yet a sampling path

The boundary experiment tests C57's two-pass formula at theta=0,pi/4,pi/2,pi
for every complete output in the fixed physical fixture. The rank-two
correction matches full-r sequential probabilities within 8.4e-17 and physical
statevec within 2.3e-14. Formula normalization error is below 4.5e-16.
Deleting interference and retaining only a same-sector boundary contribution
both fail even after normalizing their incorrect nonnegative weights: the
maximum TVs are about 0.304 and 0.274.

The corrected verifier streams sectors, recomputing each early factor in its
second pass. It still enumerates all M sectors for each requested complete
output; the small diagnostic additionally stores the tiny output law for
comparison. Neither this formula nor a low matrix rank establishes an
efficient LOW-bit marginal oracle or an orbit-independent exact sampler.

## Reports and retained failures

Authoritative main-reviewed runs:

| Experiment | Checks | Report | Log |
|---|---:|---|---|
| `experiment_symmetry_kick` | 16/16 | `out/symmetry_kick_20260911T021039609927Z.json` | `out/symmetry_kick_main.log` |
| `experiment_symmetry_kick_bounds` | 9/9 | `out/symmetry_kick_bounds_20260911T021210209038Z.json` | `out/symmetry_kick_bounds_main.log` |
| `experiment_orbit_lightcone` | 869/869 | `out/orbit_lightcone_20260911T021425398513Z.json` | `out/orbit_lightcone_main.log` |
| `experiment_orbit_cover` | 5/5 | `out/orbit_cover_20260911T021206991470Z.json` | `out/orbit_cover_main.log` |
| `experiment_symmetry_boundary` | 11/11 | `out/symmetry_boundary_20260911T021426080464Z.json` | `out/symmetry_boundary_main.log` |

The following earlier problems were preserved, not silently removed:

- The physical tester initially supplied a dimension-mismatched full-r branch,
  then a wrong-sized sector block. Production guards rejected them. Failure
  reports are `out/symmetry_kick_failure_20260911T020659595899Z.json` and
  `out/symmetry_kick_failure_20260911T020717148057Z.json`.
- The bounds tester initially omitted an angle argument to `work_mixer`;
  `out/symmetry_kick_bounds_failure_20260911T020831840906Z.json` records the
  resulting TypeError. Main subsequently replaced one-sided norm predicates
  with absolute norm-error checks and added the production cover comparison.
- The first light-cone must-fail choice r=9,s=3,m=1,S={0} decreased mass
  relative to the bare count, so the chosen underbound predicate did not fail.
  The 868/869 report `out/orbit_lightcone_20260911T020959257262Z.json` and
  `out/orbit_lightcone_run.log` remain. The wrapped control above instead
  demonstrates the false general shortcut. There is no claim that mixing
  must always increase mass on every chosen support.
- The first boundary verifier retained a list of every early sector factor
  while reporting streamed storage. Main removed that list and recomputed
  factors in the second pass. The original report
  `out/symmetry_boundary_20260911T021132196792Z.json` therefore does not
  establish the claimed storage behavior. It also labeled half-L1 differences
  of unnormalized control weights as TV. The reviewed report retains those
  raw discrepancies under explicit names and adds normalized probability-law
  comparisons. Neither repair changes the correct formula's probabilities.

Core passed before science in `out/symmetry_kick_core.log`. Full affected
lab/claims suites pass in `out/symmetry_kick_test_lab.log` and
`out/symmetry_kick_test_claims.log`. The other six science suites were not
rerun for this scalar-bound addition. Dense experiments have small dimension
caps and 16 MiB reference limits; the interval counter allocates no dense
orbit/prefix arrays. Generated records were refreshed and the documentation
gate run. No manuscript/abstract edits or commits in this follow-up.

Reproduce each script from research/ with, for example:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_symmetry_kick
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_orbit_cover
```

The other names in the table use the same module command and timestamped
reports. Primary-source bodies on causal quantum walks and unitary error
composition were read; C57 records the relevant sections and the limited
connection. No priority claim was established. The main useful lesson is that
breaking a conserved quantity and defeating an accurate approximation are
different questions. TODO 22 asks a physically natural follow-up that does
not inherit the fixed-support omission guarantee.
