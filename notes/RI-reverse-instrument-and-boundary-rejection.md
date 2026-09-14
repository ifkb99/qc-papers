---
code: RI
date: 2026-09-11
title: "A backward measurement trajectory can replace forward-state storage in exact arithmetic"
outcome: confirmed
claims: [C56, C63, C64, C67, C68]
todo: [14, 28]
---
# RI — Move the initial condition into a terminal rejection test

C68 owns the exact density/pure-trajectory proofs and costs. This is an
initial mathematical/full-law audit, NOT a certified finite-precision reverse
sampler or a measured replacement for C67. The idea arose during §UG's
stronger-baseline comparison: perhaps the positive backward effect can itself
be sampled as a work state. Main derived the proposal/boundary identity
before measurement. One lower-cost agent independently audited both proofs;
another implemented the tiny raw-effect probe. Main audited and reran it.

No forward mixing gap is required. This does not justify resetting the
backward effect: the trajectory is where its measurement memory now lives.
The primary reversal formulas actually read are cited in C68. Novelty of
this application has not been established.

## Complete-law evidence

`experiment_reverse_instrument_probe.py` fixes r9,b3,t4, the earlier W1/W3
backgrounds, both deterministic reflection declarations, and observable
W0=embedded Rx(pi/4). All four histories and three initial sectors are
enumerated at P192/P256. Width and gates stay fixed across precision.

The existing branch builder and QFT helper form UNNORMALIZED reverse effects
from exact-enclosed I/3. All 24 complete laws have proposal mass enclosing
one and accepted submass enclosing one third. At every output the accepted
interval overlaps the separately contracted joint final-sector/output target.
Bijective routes and uniform initial sectors explain why that joint target
equals the conditional target divided by three. Scalar interval widths are
below 2^-120, excluding vacuously broad overlaps. Overlap is a numerical
consistency check, not an equality proof.

The independent full-r `direct_joint` floating diagnostic agrees with target
midpoints within about 2.8e-17; it is not the certificate. Expected attempts
are derived and total accepted mass checked, not estimated from a histogram.
Omitting acceptance changes at least one conditional law by outward TV at
least about 0.235673. That lower bound charges BOTH laws' interval radii.
The normalized mixed I/3 null has constant acceptance one third, checked
without dividing at zero effects; it is not a new production input option.

Authoritative report:
`out/reverse_instrument_probe_20260911T060723268281Z.json`, 4/4 checks;
log `out/reverse_instrument_probe_audited.log`. A 16 MiB structural preflight
covers both retained precision reports, scalar strings, serialization copies
and live tiny-reference arrays. It is not RSS. This diagnostic retains full
laws and uses the full-storage builder; it does not measure streaming savings.

## Preserved verifier failures

The first reference supplied a two-dimensional background on three-state
work: `out/reverse_instrument_probe_failure_20260911T055948344721Z.json`.
Strict P1/C2 failures remain in `out/reverse_instrument_probe_*.json` reports
stamped `20260911T060001445011Z`, `20260911T060026450921Z`,
`20260911T060115723212Z` and `20260911T060158744959Z`, with paired failure
reports. The verifier confused physical and measurement-order output bits,
and Python division constructed binary64 1/3 BEFORE handing it to Acb.
Precision doubling cannot repair that rounded exact input.

Main identified the initialization error and bit-reversal signature and
required checking raw negativity before clipping to known nonnegative masses.
The agent repaired these without weakening mass-one/one-third predicates;
its first pass is `out/reverse_instrument_probe_20260911T060604813187Z.json`.
Main then charged both radii in the TV control, bounded interval widths and
serialization, retained control metadata, and renamed the misleading overlap
field `exact_identity` to `enclosure_overlap`. The stronger rerun is above.

## Regression and validation

`test_claims.py` uses a rational rotation with entries 3/5 and 4/5 and a sign
branch. Its Kraus operators are nonnormal, preventing a vacuous adjoint test.
Both completeness relations hold exactly and the pure reverse accepted mass
is target/b. Dropping rejection or the adjoint gives different rational laws.
This analytic check is not a new generic propagator or three-state sampler.

Core passed before science as recorded in §UG. The claims suite was rerun
after this regression and passes with the optional backend
(`out/reverse_instrument_claims.log`). Other scoped lab/full-law checks are
in §UG; the other six science suites were not rerun. Documentation validation
shares `out/uniform_gap_docs.log`. No manuscript/abstract edits or commits.

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_reverse_instrument_probe
```

This is the more promising lead because it removes mixing assumptions and
forward checkpoints in exact arithmetic. Whether it wins at equal certified
accuracy is unresolved. TODO 28 alone owns that next experiment; C68 keeps
numerical and coherent-history boundaries explicit.
