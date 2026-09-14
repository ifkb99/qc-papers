---
code: PE
date: 2026-09-11
title: "Prefix-sampling precision: a proved budget, with the numerical oracle still missing"
outcome: mixed
claims: [C60]
todo: [14, 24]
---
# PE — Do not divide a global claim by the smallest conditional probability

**Follow-up:** C61/§VP now implement the exact-input oracle/kernel boundary.
The account below is the preceding investigation, when that boundary was
still open; its numerical diagnostics do not become certificates retroactively.

TODO 24 made progress but remains open. C60 owns the proved sparse-block
error bound, its scaled-amplitude accuracy requirement and rejection
comparator. Main read Bravyi–Gosset–Liu's robustness lemma AND its complete
supplemental proof, then derived a direct post-prefix-oracle version for
this harness. This is standard error-analysis reasoning applied here, not
a novelty or breakthrough claim.

The useful distinction is between worst conditional error and error in the
overall sampled law. A nearly zero block can have a very wrong normalized
direction while its mass-weighted contribution is negligible. Conversely,
equal block masses do not establish correct directions. The bound uses
the FULL Born laws, not just block marginals. One agent initially interpreted
the symbols as block marginals and proposed a refutation; clarifying that
definition made its example a valid negative control, not a theorem failure.

Three lower-cost agents worked on rational block inequalities, output-grid
quantization and native numerical precision. Main implemented the rational
budget helper, audited all three scripts, completed the quantization script
after repeated verifier errors, reran the final experiments and added tests.
The qsim-research skill's control, normalization and input-contract rules
prevented several premature claims of certification.

## Exact-rational checks and what the planner actually requests

The block experiment checks 392 seeded rational partition cases with at most
12 coordinates, including the one-step Markov-kernel inequality. Its explicit
zero-block control has ideal mass 1/100, uniform-fallback weighted TV 1/200
and full-law TV 1/100. Equal-mass opposite-direction blocks have full TV one
despite zero block-marginal TV. Another 50 exact-rational cases check the
accepted-measure normalization bound for rejection, including zero mass.
They test a conditional error theorem, not actual floating rejection draws.

`lab/sampling_error.py` contains exact-rational budget and accuracy-planning
helpers. For t=63,b=2,135 block updates and requested TV<=10^-6, the conservative
planner asks for real/imaginary scaled-coordinate errors at most 2^-62 and
29 unbiased random bits per exact-CDF block draw. Its computed total budget
is about 6.071e-7. These are ABSOLUTE amplitude-accuracy and ideal RNG requests,
not a claim that a 62-bit machine mantissa suffices. The helper explicitly
returns `certifies_existing_float_sampler=False`.

The oracle requirement ranges over all labels, not only those visited in
a diagnostic. The dimension factor costs linearly more accuracy bits in t;
obtaining those bits with certified gate inputs, phase evaluation and
contraction arithmetic remains the implementation gap.

## Output-grid rounding: complete tiny transition laws

Keep r=10,b=2,t=4, background Rx(pi/2) at insertions 1,3, and coherent
q=0@s2,q=1@s3 rotations of strength pi/4. Quantize each real/imaginary
sqrt(M)-scaled prefix amplitude to the nearest multiple of 2^-p. Sweep
p=0,1,2,4,6,8,10,12. This is rounding the existing complex128 oracle's
OUTPUT, not changing its internal arithmetic or proving its input accuracy.

The finite reference enumerates all sampling transitions and, separately,
every unique coordinate of each idealized post-prefix amplitude function.
Partial QFT coordinates contain low remaining exponent labels and measured
output prefixes once each. The rounded kernel declares uniform-within-block
fallback when all rounded weights are zero. The production sampler is
unchanged and does NOT silently gain this fallback.

All reference/rounded transition laws normalize. The unrounded transition
law agrees with the existing independently constructed full-r joint law to
TV below 2.2e-16. The derived weighted, full-Born and global-L2 inequalities
hold on every tested row, including the scaled uniform-grid error bound.
Selected results:

| absolute grid bits p | output TV | summed weighted conditional TV | twice summed full-prefix Born TV, capped at one |
|---:|---:|---:|---:|
| 0 | 0.572443 | 1.62155 | 1 |
| 4 | 0.0378565 | 0.101046 | 0.692082 |
| 8 | 0.00216204 | 0.00612763 | 0.0447583 |
| 12 | 0.000118193 | 0.000427632 | 0.00244870 |

At p=12 the largest measured local conditional TV is still 0.5. Its recorded
parent mass is about 3.08e-34: a numerical-null event in this floating reference,
not a certified nonzero physical probability. It illustrates why that maximum
alone is a poor distribution-level diagnostic. No monotone-rounding theorem
is inferred from the finite precision rows.

The deliberately discard-only chain is actually propagated; it loses all
mass at p=0 and retains about 0.903889 at p=2. Those subnormalized arrays are
NOT called probability laws or compared using TV. An independent normalized
64-coordinate control has amplitude error 1/8 but TV 1/2, refuting the false
dimension-free bound TV<=maximum coordinate error.

A normalized rounded-final-amplitude table is also compared to its global
L2 bound. This is explicitly not the output law of an implemented approximate
rejection sampler. The rejection comparison in C60/its rational experiment
states the separate accepted-measure error promise still needed.

## Native precision and the input-error floor

The separate fixed fixture uses r=10,b=2,t=4, Rx(pi/7) at insertions 1,3 and
q=0@s2,q=1@s3 rotations of strength pi/5. Rational multiples of pi are the
mathematical input intent; trigonometry is still numerical. Identical finite
full-r products are evaluated at complex64, complex128 and native clongdouble.
On this platform NumPy reports 23,52,63 fraction bits respectively (not the
storage-format name as a precision guarantee). The script checks matrix and
phase dtypes. Long double is a stronger numerical reference, not exact truth.

Normalized output TV against the reconstructed long-double reference is
about 3.33e-8 for complex64 and 4.69e-17 for complex128. The complex128
normalized marginal matches production within 4.17e-17 per output. Raw
complex64 mass is 1.0000001769512892 and its raw half-L1 gap is about 8.85e-8;
that latter number is NOT the normalized TV. Native-precision sums are kept
through normalization before serializing report numbers to ordinary floats.

Widening stored complex128 WORK-gate coefficients to clongdouble, while
evaluating QFT phases at target precision, leaves normalized TV about 1.69e-17.
Its matrix unitarity defect is about 1.08e-16 versus 4.29e-20 for reconstructed
long-double inputs. More decisively, the stored Rx column's squared norm
minus one is a nonzero EXACT rational number:

    6149843566922181 / 324518553658426726783156020576256.

Merely widening storage cannot remove that input error. This does not imply
that all effects of coefficient rounding are observable or establish a
universal error floor for every circuit/output.

## Failed verifiers and corrections retained

- The quantization scaffold originally shared a cache whose first unrounded
  query prevented later rounding, and shifted work labels regardless of the
  controlling exponent bit. Main caught both. Early reports
  `out/prefix_precision_20260911T031501417682Z.json`, `...031522965930Z.json`
  and `...031536900639Z.json` include failed controls and an undefined-variable
  exception. They are not evidence for the final algorithm.
- `...031733684132Z.json` aborted on a zero rounded full-prefix vector.
  Explicit normalized zero-block/global-zero semantics are now stated.
  The later scaffold also overwrote its fallback flag and scaled the reference
  by sqrt(M) twice, quantizing a different function than the sampler. Main
  took over, corrected both and reran. The original single-stage error control
  was unsupported; the final independent analytic control is specified above.
- A discard-control mass was initially assigned from whether fallback occurred
  rather than measured. Main replaced it with actual discard-chain propagation.
  The script now gates the intended inequalities and calls its final-amplitude
  comparator by its actual task, not an approximate rejection sampler.
- The precision-input experiment's first report
  `out/precision_inputs_20260911T031604103421Z.json` fails normalization and
  production comparison. Later early PASS reports, including
  `...031641306995Z.json`, still had extra background insertions absent from
  production and mislabeled raw half-L1 errors as TV. Main required matching
  schedules, separate raw/normalized metrics and native-precision normalization.
  Their apparent agreement did not justify overlooking the different circuits.
- The rational audit initially returned a worst-case zero-block error while
  labeling it the measured uniform-fallback error. Main corrected that to the
  actual conditional TV and added timestamped reports and normalization guards.

## Reproduction and next boundary

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_block_error_bound
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_prefix_precision
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_precision_inputs
```

Main-audited final reports:

| experiment | checks | report |
|---|---:|---|
| rational block/error budgets | 6/6 | `out/block_error_bound_20260911T032244838032Z.json` |
| full rounded transition laws | 5/5 aggregate | `out/prefix_precision_20260911T032216641658Z.json` |
| matched native precision | 4/4 | `out/precision_inputs_20260911T032137214322Z.json` |

Logs: `out/block_error_bound_main_audited.log`,
`out/prefix_precision_main_audited.log`, `out/precision_inputs_main_audited.log`.
Core/lab/claims pass in `out/prefix_precision_core.log`,
`out/prefix_precision_test_lab.log`, `out/prefix_precision_test_claims.log`.
The other six science suites were not rerun for this isolated helper addition.
The documentation gate passes all ten checks in `out/prefix_precision_docs.log`;
indexes were regenerated and `git diff --check` is clean.
No manuscript, abstract workshop or production sampler semantics changed;
no commits or publication. TODO 24 stays open: next implement and test the
verified amplitude oracle under explicit exact input semantics, then connect
its certified absolute errors to a finite-bit categorical kernel. The broader
user research goal remains active.
