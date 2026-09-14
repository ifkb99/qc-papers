---
code: PS
date: 2026-09-10
title: "Periodic defects: coarse-sector symmetry handles several noncommuting mixers; odd blocks still have input-specific cancellations"
outcome: mixed
claims: [C56]
todo: [14, 20, 21]
---
# PS — A stronger baseline appears before the rejection sampler

TODO 20 proposed extending the localized sampler to repeated blocks covering
the whole orbit, testing whether output effects survive while the conditional
row remains a short Fourier sum. The lower-cost initial tests supported that
formula. Main then noticed the stronger common symmetry: all the repeated
blocks preserve eigenspaces of U^b. Sampling that coarse phase first retains
a small work register even with several separated noncommuting mixers.
C56 owns the derivation, sufficient conditions, costs and limits.

The qsim-research workflow changed the implementation plan here. The stronger
simple baseline displaced the planned single-defect rejection helper, and the
independent-reference rule prompted replacement of an agent's new full-state
loop with the existing Circuit/statevec engine. A probability-level control
also replaced a vacuous amplitude comparison. Passing formulas alone were not
treated as evidence for a new simulation barrier or breakthrough.

## Work split and authoritative runs

Three `gpt-5.6-luna` agents handled bounded initial formula, invariance and
coarse-sector tests. One then tested main's production sampler against an actual
modular circuit. Main implemented `lab/periodic.py` and `sequential_path` in
`lab/semiclassical.py`, audited the agents' scripts, strengthened allocation
guards and controls, reran the experiments and added regression tests. A
follow-up agent challenged the input-specific cancellation discovered during
that review. The scientific code and records remain uncommitted.

Final main-reviewed experiments:

| Experiment | Checks | Coverage |
|---|---:|---|
| `experiment_periodic_formulas` | 41/41 | harmonic amplitudes, phase weights, normalized accepted law and interference/feedback/divisibility controls |
| `experiment_periodic_invariance` | 48/48 | binary thresholds, fixed odd-block impact series, spectral/prefix and extended-precision checks |
| `experiment_periodic_sectors` | 38/38 | explicit sector projections, existing Circuit/statevec references, production conditional paths and invalid-dephasing controls |
| `experiment_periodic_sampler` | 12/12 | physical modular circuit with three defects, ideal/end controls, original schedule, wide actual draws and resources |
| `experiment_periodic_input_cancellation` | 17/17 | reached-span cancellation, arbitrary complement blocks, divisibility/leakage controls and a retained commuting leakage example |

Their authoritative report/log pairs are:

- `out/periodic_formulas_20260911T014441596588Z.json` and
  `out/periodic_formulas_main_run.log`.
- `out/periodic_invariance_20260911T014433811255Z.json` and
  `out/periodic_invariance_main_run.log`.
- `out/periodic_sectors_main_corrected.json` and
  `out/periodic_sectors_main_corrected.log`.
- `out/periodic_sampler_20260911T014818251729Z.json` and
  `out/periodic_sampler_main_run.log`.
- `out/periodic_input_cancellation_20260911T020207788747Z.json` and
  `out/periodic_input_cancellation_main.log`.

Earlier agent reports are preserved but do not supersede these reviewed runs.
In particular, `out/periodic_sectors.json` predates the reference/control audit,
and `out/periodic_sampler_20260911T014520164314Z.json` predates the schedule
control correction and fuller payload accounting.

## What the bounded tests actually found

The harmonic experiment fixes seeded complex W, b=3, r=6, t=8 and steps
s=0..8. Separate controls include r=9 and identity W. Joint errors against
SparseOrbitPrefix stay below 6.5e-15, and spectral marginal errors below
5.5e-15. The normalized accepted early-output law is additionally checked by
enumerating component proposals and acceptance, with maximum error below
3.8e-14. It is not just a restatement of the averaged rejection bound.
The pointwise rejection cost reaches about 56 despite an averaged cost of 3;
this supports the need for the averaging qualification in C56.

The invariance experiment resolves the binary threshold one insertion at a
time. The b=2,r=6 and b=4,r=12 rows are visibly nonideal before their predicted
thresholds and agree with ideal afterward to about 1e-14. A fixed b=3 block at
t=8,s=3 gives output TV between about 0.465 and 0.592 for r=3,6,...,18.
These are finite visible examples, not an asymptotic non-dilution theorem or
a universal statement about odd blocks. The commutator, harmonic formula,
spectral/prefix and float64/longdouble comparisons test the mechanism, not
merely the output size.

The revised sector experiment compares explicit orbit-space projections of
U^n and W against the sector matrices. Specialized b=2 and b=3 gates in the
existing Circuit engine supply an independent statevec reference, with at most
t=6,r=12. Single-defect mixtures agree with original-order spectral effects
within 4.5e-15; the two-defect full-orbit comparison is below 1.7e-16.
Production conditional probabilities match the independently constructed
sector circuits within 8.4e-15. The wrong-twist test compares matrices at a
nontrivial phase: simply reversing alpha's sign can relabel a uniform mixture
and is not necessarily an observable-level falsifier.

The physical test uses N=7,a=3 and its clean orbit [1,3,2,6,4,5]. W01 mixes
local block positions 0/1, W12 mixes 1/2, both by Rx(pi/2). They have exact
three-work-qubit Pauli-rotation constructions; full matrices and statevec
actions are checked, including unchanged invalid labels 0 and 7. With t=5
and W01 after control count 1, W12 after 3, W01 after 4, summed coarse-sector
joints match the actual compiled modular circuit within 2.24e-14. The
physical reference's norm error is 1.69e-13. This is a small physical
realization of these particular blocks, not a locality theorem for arbitrary W.

Removing inverse-QFT feedback changes that physical distribution by TV 0.347;
retaining only alpha=0 changes it by 0.256. Reversing just the arithmetic/defect
schedule, with preparation and final QFT fixed, changes it by 0.0868. The
defect-free arithmetic-reversal control remains valid within 2.9e-14.
These comparisons specifically guard the original time-order contract.

At t=63, eight actual draws each were made for r=6 and abstract supplied
r=3,000,000,021, using b=3 and defects after 16,32,48 controls. Per path,
forward-factor payload is 9,216 bytes and the reported branch-matrix payload
is 18,144 bytes; supplied blocks plus identity occupy another 576 bytes.
Those are payload categories, NOT peak process memory. NumPy temporaries,
interpreter overhead and input/setup must not be erased by quoting only the
smallest category. The implementation stores O(t*b^2) scalars, not constant
total memory. The large period is given; no factoring/order search was done.
Wide seeded draws exercise execution and guards, not every rare probability
or an arbitrary-precision guarantee.

The independent input-cancellation challenge uses seven seeded genuinely
complex block-diagonal cases, including L=1 and L=r and arbitrary unitary
complement blocks. They agree with ideal to below 3.1e-15. A first-span block
with r=6,L=4 (non-dividing) instead gives TV about 0.553. At r=8,L=4, a
rotation between labels 0 and 5 gives TV about 0.0754. These controls exhibit
failures when the sufficient hypotheses are removed; they are not necessary
conditions for invisibility. Main reviewed the spectral reference use, added
pre-allocation guards, and retained the initially invisible leakage case as
the explicit commutation test described below.

## Failed predictions and audit corrections

1. The first invariance run indexed an early row beyond its length when L<b.
   Its preserved failure report is
   `out/periodic_invariance_failure_20260911T013734915333Z.json`.
   The periodic multiplier is now evaluated directly at each residue rather
   than inferred from a prefix too short to contain it. A second run failed
   JSON serialization of complex metadata, recorded in
   `out/periodic_invariance_failure_20260911T013750986410Z.json`; real/imaginary
   payloads repaired reporting, not the scientific formula.
2. The initial sector script used a second Q-by-b propagation loop. It was
   replaced with the existing Circuit/statevec reference, as required by the
   repository's independent-reference rule. Its b=2 inverse-shift control had
   a large amplitude difference but unchanged measured probabilities. The
   revised b=3 probability-level counterexample differs by about 0.0193.
   Agent audit reports `out/periodic_sectors_audit_20260910_v2.json` and its
   matching log retain this correction; the main rerun agrees. Dropping the
   twisted-cycle wrap phase is a separate control from dropping inverse-QFT
   feedback, despite both being called feedback in some raw labels.
3. The initial physical sampler control reversed ALL gates, including state
   preparation and QFT. Its large discrepancy was not clean evidence about
   arithmetic/defect reordering. Main replaced it with the schedule-only
   comparison above and added the defect-free positive control. Main also
   removed a duplicate instrument implementation from the no-QFT-feedback
   control, using the actual Circuit with only its Fourier phase gates removed.
4. Main's first claims regression assumed a specific odd b=3 W01 at s=1
   must be visible. It failed, preserved in `out/periodic_test_claims_initial.log`.
   For r=6,t=4, that insertion is ideal within roundoff, whereas s=2 is
   visibly different. The spectral reference agrees with the sampler. The
   reached two-label span is preserved at s=1 and L=2 divides r: this led to
   C56's separate input-specific cancellation proof. The corrected regression
   tests BOTH cancellation at s=1 and visibility at s=2, rather than deleting
   the failed scientific prediction from the record.
5. Main initially invoked the revised sector experiment as a script instead
   of the documented module command. It stopped before science with
   `ModuleNotFoundError: circuits`; `out/periodic_sectors_main.log` is retained.
   The corrected `python -m experiments.experiment_periodic_sectors` invocation
   passed without a code change.
6. The input-cancellation agent's first leakage control rotated labels 0 and
   4 at r=8,L=4 and was invisible. Its failed 15/16-check report is
   `out/periodic_input_cancellation_20260911T015834046319Z.json`; the associated
   exception report is
   `out/periodic_input_cancellation_failure_20260911T015834050294Z.json`.
   Changing the selected unreached label to 5 gives the visible control above.
   Main additionally explained and retained the ORIGINAL 0/4 example: its
   two-level rotation commutes with U^4, so moving it to the traced end is
   valid despite its violation of span preservation. The final run verifies
   both the zero commutator and ideal output. This is a sufficient-versus-
   necessary condition distinction, not numerical noise or a sampler bug.

## Validation and interpretation

Core passed before the investigation (`out/periodic_test_core.log`). Full
affected lab and claims suites pass in `out/periodic_test_lab_initial.log`
and `out/periodic_test_claims_corrected.log`. The other six science suites
were not rerun for these helper additions. Dense diagnostic arrays are bounded
before allocation, with small r/t caps and 16 MiB reference ceilings. Production
uses its separate dimension/width and conservative payload guards. No GPU is
needed for these tiny checks.

Reproduce from research/ using the module invocation, for example:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_periodic_sampler
PERIODIC_SECTORS_REPORT=out/periodic_sectors_reproduction.json OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_periodic_sectors
```

The other experiment names in the table use the same module invocation.
Timestamped reports preserve earlier evidence; choose a fresh sector-report
path because that script's default path is fixed. Main regenerated the ledgers
and ran the documentation gate after writing this record. Neither manuscript
nor either abstract workshop was edited in this follow-up.

The interesting link is to symmetry sectors plus sequential ancilla/MPS
sampling and forward-state/backward-effect inference, whose primary sources
and exact scope are recorded in C56. This is a useful enlargement of the
known-order tractable family and a stronger baseline for future tests, not
evidence of a new general-purpose simulation method. TODO 21 owns the next
discriminating question, beyond the common symmetry; larger period records
alone would not test it.
