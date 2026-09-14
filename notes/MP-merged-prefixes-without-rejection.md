---
code: MP
date: 2026-09-11
title: "Merged intermediate prefixes remove rejection without a clear runtime win"
outcome: confirmed
claims: [C59, C60, C71, C72]
todo: [14, 24, 31, 32]
---
# MP — The intermediate oracle is enough

C72 owns the proof, input promise and charged resource bounds. Main implemented
`lab/merged_prefix.py`, reusing the sparse adjoint helpers and the UNCHANGED
gate-by-gate sampling loop. Lower-cost agents supplied initial formula,
complete-law and comparison experiments; main audited and strengthened their
predicates before the authoritative reruns below. This is an opt-in float
diagnostic, not numerical certification or a new general sampling principle.

## Complex prefixes and complete laws

`out/merged_prefix_formulas_20260911T073945582622Z.json` passes 6/6 checks.
The 2,745 valid API requests include every tiny arithmetic/background/reflection
boundary and partial-QFT request, endpoints t=0/1, odd/binary blocks, generic
three-label routes, a pi-limit fixture and genuinely exact zero amplitudes.
Largest merged/history complex residual is 2.84e-16; the existing independent
full-r prefix reference agrees within 2.80e-16. Full joint laws also agree.
History enumeration was patched to raise during merged calls. Wrong boundary,
missing prepared-control tail factor and wrong internal QFT phase produce
large complex discrepancies. Guard checks precede dense allocations.

`out/merged_prefix_sampling_20260911T074213009424Z.json` passes 7/7 checks.
The existing transition enumerator gives complete joint sector/output laws,
not sampled histograms, for five tiny fixtures. All agree with the independent
full-r circuit law, with TV at most 1.67e-16; merged joint evaluation gives
TV at most 2.01e-16. An independent reached-SET recurrence verifies every
reported prefix/block/QFT/reflection counter. Removing intermediate work
mixers changes the main law. Wrong magnitudes, boundary inclusion and INTERNAL
QFT phase change the complete law by TV approximately .0161, .0920 and .0600.

An initially failed negative control was instructive: taking componentwise
absolute values of the COMPLETED oracle output preserves all Born weights.
It therefore cannot detect a sampling error, although it breaks the complex
amplitude API. An agent temporarily replaced it with a constant-coordinate
corruption to force failure; main removed that contrived replacement and
retained the original transformation as a positive invariance check. Internal
phase corruption remains a meaningful must-fail control. C72 states the
distinction; the skill's nonvacuous-control rule materially changed this test.

## Frozen serial comparison

Main's authoritative serial report is
`out/merged_prefix_comparison_20260911T074906529943Z.json` (4/4 checks).
It includes 120 timed rows, eight preflights and a tiny static full-r law
check. The latter is a correctness baseline, NOT a sampler timing result.
All scientific test jobs had finished before this timing run.

Inputs are the separately frozen adjacent and separated fixtures inherited
from CM: supplied r=3*(2^40-1), b=3, t=15, noncommuting work blocks, fixed
insertion positions and pi/4 coherent angles. Only k varies within each
fixture. Three seeds are retained individually; proposal caps are 2,048.
Medians below are seconds for setup, adapter construction and ONE returned
sample, with counters/stat extraction included. These are matched FLOAT tasks,
not algorithms certified to the same requested TV.

Separated labels q0=0, q1=679535556937:

| k | Merged prefix | Sparse reverse | History gate-by-gate | History rejection |
|---|---:|---:|---:|---:|
| 0 | .00423 | .00125 | .00318 | .00111 |
| 2 | .00502 | .00283 | .01093 | .00259 |
| 4 | .00795 | .00983 | .04095 | .01506 |
| 6 | .01361 | .01464 | .15250 | .09772 |
| 8 | .02321 | .02405 | .60058 | 1.08407 |

For adjacent q1=1 at k=8 the respective medians are .02188, .02268, .57401
and 1.48938 seconds. Do not pool the alphabets or interpret three seeds as a
runtime theorem. Sparse rejection varies materially with the random retry
count; the separated k=8 range is .01052–.04618 seconds versus merged
.02259–.02554. The result supports rough parity with the strongest tested
exact-arithmetic-method float implementation, not a clear new speedup.

At separated k=8 both merged methods retain at most 16 reached labels, below
the global bound 17. The declared working allowances are 1,134 complex slots
for batched prefixes versus 624 for reverse vectors, each with 153 supplied
gate slots. These include arithmetic temporaries but not Python overhead,
integer storage or native RSS. All merged calls enumerate zero histories;
prefix sampling also has zero rejection proposals. Counters instrument the
OWNED input snapshot and every actual query/retry. Charging one prefix or
ignoring retries is independently demonstrated to understate total work.

C71's neighboring-route approximation remains a stronger simple alternative
when its error is acceptable. The low-k input/static simplifications recorded
in CM also remain relevant. Large supplied r alone is not evidence of hardness;
neither construction includes order/index discovery or arbitrary physical
gate compilation. Defaults and manuscripts are unchanged.

## Preserved audit failures

- The first formula report `out/merged_prefix_formulas_20260911T073711825148Z.json`
  miscounted a diagnostic extra query. A subsequent agent pass still called a
  history joint calculation "direct" and classified tolerance-small entries
  as exact zero. Main added the existing independent direct joint reference,
  strict zero counts, a genuine zero fixture and pre-allocation guards.
- The sampling reports preceding the authoritative main run retain the false
  phase-control prediction and subsequent verifier revisions. In particular,
  `out/merged_prefix_sampling_20260911T074023086104Z.json` is the agent pass
  with the replacement control, not the final audited evidence above.
- `out/merged_prefix_comparison_20260911T074456055776Z.json` preserves an
  unexplained preflight history-rejection TypeError involving None and complex.
  Its 120 full rows passed, as did an agent repeat and main's serial rerun.
  No per-row traceback was recorded then; the cause is NOT established as an
  instrumentation error. Main added per-row traceback capture for recurrence.
- `out/merged_prefix_lab.log` retains a test failure from trying to modify a
  read-only gate array. The intended input-snapshot test now replaces the
  dictionary entry and verifies that the adapter's owned copy is unchanged.

No failed raw report was deleted or silently relabeled as a passing run.

## Validation and next question

Core ran first. Main's passing logs are `out/merged_prefix_core.log`,
`out/merged_prefix_lab_final.log`, `out/merged_prefix_lab_no_flint.log` and
`out/merged_prefix_claims.log`. The backend-absent run explicitly skips verified
arithmetic checks. Lab includes all tiny prefixes, exact operation counts,
history-forbidden sampling, snapshot ownership and invalid requests. Claims
includes an exact integer-complex batched-adjoint identity and phase controls.
The other six science suites were NOT rerun. Documentation validation is
recorded separately in `out/merged_prefix_docs.log`.

Reproduce from research/ with Python 3.12 and one BLAS thread:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_merged_prefix_formulas
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_merged_prefix_sampling
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_merged_prefix_comparison
```

Main read the primary BGL Algorithm 2 induction, adaptive extension and complete
Lemma 1 robustness proof linked in C72; not every application in that paper
was audited. Observed float residuals do not establish that robustness contract.
Prior sparse-simulation positioning remains in C71.

TODO 31 is complete at this bounded scope. TODO 32 owns the next discriminating
test: fixed-alphabet support growth and a safely separated larger-k input API.
Its initial algebra is a lead, not an implemented result at this checkpoint.
The open-ended research goal remains active.
