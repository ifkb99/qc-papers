---
code: BC
date: 2026-09-09
title: "Balanced nonlinear controls project exactly; biased controls have a certified geometric residual"
outcome: established
claims: [C47]
todo: [15a]
---
# BC — Balanced controls and the reduced channel

First investigation in the user-requested four-part sequence. The small toy
was initially checked interactively; the saved experiment then tested the
derived full-support and reduced-error formulas on an actual arithmetic block.
C47 is the authoritative theorem and scope statement.

`experiments/experiment_balanced_controls.py` passed **54/54 checks**. It uses
the existing Walsh replay to extract both branches of N=7's identity-tail
block, verifies that the inactive branch is identity and the active branch an
involution over all work states, and checks the one-control reduction against
the existing gate-level PPS. Both O+ and O- are nonzero, so agreement is not
an invariant-observable or zero-operator artifact.

For each fixed work map, 3-bit parity, majority, OR and constant controls are
compared at consecutive window counts: 1..4 for the CNOT toy, 1..3 for the
arithmetic block. Every reduced entry matches the formula with error zero,
and every complete Walsh support count matches the product-factorization
formula. The arithmetic calculations include arbitrary scratch inputs, not
only the legitimate initialization.

At three windows the toy full supports are 4 (parity), 130 (majority), and
1024 (OR); parity and majority contract to the same nonzero two-term operator.
OR's exact projection error is (3/4)^K, with alternating residual sign;
the constant-control negative case never decays. This supplies an explicit
tail-approximation error guarantee in the stated block model, not a general
incremental-PPS truncation guarantee.

For the actual arithmetic block at three windows, parity and majority have
full support 2606 and 84695 respectively, but the same 1303-term reduced
operator. OR has full support 667136 and reduced support 1303, with a different
reduced vector and exact projection error 0.421875. Equal reduced counts alone
would therefore have missed the OR negative control.

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_balanced_controls
```

The complete rows, predictions, check results, interpreter and NumPy versions
are saved in ignored `out/balanced_controls.json`. The harness now supports
`finish(report_path=..., rows=..., metadata=...)` and writes failed runs before
raising; the behavior is covered by `test_lab.py`.

What remains beyond this completed first pass: synthesize matched low-level
window-control circuits and measure their temporary peaks and construction
cost. The present claim is about exact block maps and final/reduced spectra.
