---
code: LF
date: 2026-09-10
title: "Finite-support defects: arithmetic-progression marginals and averaged interference rejection remove the prefix vector"
outcome: mixed
claims: [C55]
todo: [14, 19, 20]
---
# LF — A finite probability formula becomes an actual sampler

The active user goal asks for continued research and lower-cost initial tests
of promising directions. TODO 19 asked whether C54's exponential early-prefix
array was avoidable for a defect with finite known orbit-basis support. The
qsim-research workflow made the actual sampling step, normalized rejection
cost, and independent circuit comparison necessary gates. C55 owns the proof
and limits; this note records the investigation and caught failures.

## Division of work and validation

Two `gpt-5.6-luna` agents received independent finite-formula and progression
tests; a `gpt-5.6-terra` agent received implementation. All eventually hit a
usage limit. Only the progression agent left a completed initial report:
`out/progression_prefix.json` and `out/progression_prefix_run.log`, with its
finite-sum/FFT/circuit and must-fail checks. The main agent inspected the actual
files, finished implementation, ran and corrected the unfinished formula test,
and added production-helper checks. No implementation was recovered from the
third agent. Failed agent turns are not counted as completed verification.

After the reported usage-reset time, both lower-cost initial testers resumed
successfully. One independently audited the production helpers without editing
them; see `out/localized_independent_audit.md` and its raw checks log. It found
no bounded counterexample in signs, gcd lifting, component normalization or
the averaged cost. Main reconciled inconsistent near-cancellation numbers and
an ambiguous displayed amplitude in its first prose summary against the
preserved log; the final audit removes those unsupported statements. The other
agent ran the separate impact probe discussed below.

The reusable additions are `lab/fourier_sampling.py` (scalar low-bit Fourier
marginals and progression sampling) and `lab/localized.py` (one finite-support
defect). They are specialized distribution algorithms, not another generic
propagator. Small independent references use the existing FFT, statevec,
SparseOrbitPrefix and original-order spectral effects.

Final main-reviewed runs:

| Experiment | Checks | What it checks |
|---|---:|---|
| `experiment_progression_prefix` | 3208/3208 | interval low-bit marginals, all forced production paths, statevec/FFT, progression gcd lifts and five must-fail controls |
| `experiment_localized_formulas` | 65/65 | finite amplitudes/norms, spectral/prefix references, actual original-order coherent arithmetic, zero/sign/complex cases |
| `experiment_localized_sampler` | 26/26 | full tiny joint distributions, normalized rejection, enumerated progression sampling branches, wide actual samples and three must-fail controls |

The progression run contains 210 interval rows and 1533 progression rows;
these are bounded coverage, not separate scaling claims. Its production-vs-FFT
error is below 2.6e-15. The direct geometric reference's maximum prefix error
is below 4.7e-14. Preserved final report:
`out/progression_prefix_production_v2.json`.

The localized formula series fixes r=6, t=8, theta=pi/2 and varies insertion
s=1..7 by one. All seven original-order physical circuit comparisons pass,
with output error below 4.6e-14. Separate abstract cases cover periods three
and five, L<r and L>r, seeded complex blocks and an actually zero-weight phase.
Final report: `out/localized_formulas_20260911T011830224671Z.json`.

The sampler's full-joint error versus the older prefix sampler is below
5.6e-17; its marginal error versus spectral effects is below 3.7e-15. The
normalized proposal/rejection error is below 5.6e-16. This includes enumeration
of positive-probability progression sampler paths, not just evaluating their
claimed output probabilities. Finite seeded sampled outputs alone would not
establish a sampling distribution.

At total width 63, the sampler draws 16 outputs each for r=6 and s=2,16,32,62,
plus a separately labeled abstract r=1,000,000,007, s=32 row. The supplied/copied
two-by-two block is 64 bytes; that is **matrix payload, not process peak RAM**.
There are no orbit-, prefix- or output-sized arrays. These wide paths do not
validate every rare output, provide a bit-complexity theorem, or simulate a
compiled billion-period modular circuit. Known order and known orbit indices
are inputs. Final report: `out/localized_sampler_20260911T012045676548Z.json`.

The must-fail controls detect deleted interference, wrong progression counts,
missing feedback, wrong bit order/gcd reduction and uniform final eigenphase
weights. A near-cancelled final phase needs about 3e8 early proposals
conditionally, explicitly refuting a tempting POINTWISE constant bound. C55's
mean bound survives because it averages using the correct phase weights.

Core, lab and claims regression suites pass under Python 3.12 / NumPy 2.4.6:
`out/localized_resume_core.log`, `out/localized_test_lab_initial.log`, and
`out/localized_test_claims_initial.log`. Other science suites were not rerun
for these isolated helpers. No existing generic circuit engine was modified.

## Preserved failures and audit fixes

The unfinished formula verifier divided probabilities by Q*r instead of
Q^2*r in BOTH its formula and direct reference. Those two routes agreed but
their total mass was Q. Independent normalized spectral/prefix/circuit routes
caught the error. Its first run passed only 21/64 checks. Preserved artifacts:

- `out/localized_formulas_initial_run.log`;
- `out/localized_formulas_20260911T011542786207Z.json` (all checks);
- `out/localized_formulas_failure_20260911T011542793795Z.json` (terminal failure).

Main corrected both denominators, made the direct row use actual unitary
columns rather than the same delta decomposition, fixed a tuple-unpacking
error in the reference sparsity bound, removed amplitude/support cutoffs, and
changed the incorrectly constructed "zero-weight" block to one whose first
column really sums to zero. The zero phase now has an explicit assertion.
Omitting zero-weight phases from the averaged cost also requires omitting
their T_k terms from the equality audit, while retaining the stated inequality.

A separate main-added allocation guard caught an unused G_(T+1) direct sum
when its prefactor e was zero. The failed production-verifier log is preserved
as `out/progression_prefix_production_run.log`; there is no complete JSON for
that interrupted run. The corrected verifier skips that zero-weight term and
keeps the original cap. It did not enlarge the reference allocation budget.

These are verifier/hygiene defects, not counterexamples to the now separately
derived normalized formula. Preserve their record: unanimous agreement between
two related implementations would have been misleading here.

## Interpretation

### Impact probe: retain the negative result

The lower-cost impact agent compared the endpoint mixer against the ideal
sampler at fixed t=8, s=5, theta=pi/2 for abstract orders r=3..16. It derived
the hit-fraction trace-distance bound now recorded in C55. Main reviewed the
code and reran it after removing support cutoffs, making the byte guard use
non-overflowing Python integers, and moving the zero-angle/endpoint controls
from a resonantly invisible order to r=6, whose interior effect is visible.

The initial report `out/localized_impact_20260911T012439735167Z.json` left P2
unresolved despite passing its other checks. The agent's second report,
`out/localized_impact_20260911T012506599599Z.json`, added a passing comparison
of just the sweep's endpoints. That does NOT establish a decreasing trend:
the last point is a binary resonance. Main replaced the predicate with every
adjacent comparison and explicitly labeled the monotone-dilution intuition
as exploratory and refuted. The corrected experiment intentionally exits 1:
**17/18 checks pass; P2 fails**. The normalized reference/bound checks and
nontrivial endpoint/zero-angle controls still pass. Preserve this as evidence,
not a reason to weaken the test until it turns green.

The observed TV is about 0.214 at r=3, below numerical resolution at r=4,8,16,
then rises again after the first two resonances. There are also increases from
r=12 to 13 and 13 to 14. Its maximum reference discrepancy is below 3.7e-15;
the bound is never violated. A small hit fraction can justify the simple ideal
approximation, but these finite rows do not establish an asymptotic rate, and
at fixed prefix length the hit fraction need not vanish as r increases.

Audited artifacts: `out/localized_impact_audited_run.log`,
`out/localized_impact_20260911T012918979096Z.json`, and the terminal
`out/localized_impact_failure_20260911T012918984381Z.json`. The failure is the
scientific monotonicity prediction, not the localized sampler identity.

### Positioning and reproduction

This is a promising structured extension: the actual sampler now removes the
remaining exponential prefix vector under an explicit stronger promise. It is
not a factoring breakthrough, and priority has not been established. The
terminating-QFT context in C55 is established prior art. The next impact and
generalization question is tracked only in TODO 20, not duplicated here.

From the research root, use the following prefix for each module:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_localized_formulas
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_progression_prefix
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_localized_sampler
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_localized_impact
```

The last command is expected to exit 1 for the recorded P2 falsification.
The others should pass. Timestamped reports preserve prior runs; the progression
script has a fixed final report path, so preserve that file before rerunning
if investigating a changed prediction or implementation.
