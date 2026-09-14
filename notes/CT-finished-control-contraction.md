---
code: CT
date: 2026-09-09
title: "Exact finished-control contraction: bounded reduced support in both order branches and an idempotent β=1 tail"
outcome: established
claims: [C45, C46, C42, C18]
todo: [12c]
---
# CT — Exact finished-control contraction

The independent review reopened TODO 12c. Naive pruning whenever a term
acquires exponent support is wrong, but a qubit's **last reverse use** is a
safe certificate. C45 contains the proof and scope; C46 gives the separate
idempotent-tail simplification. Both were derived before the sweep.

## Predictions and controls

The new experiment declared that the final reduced coefficient vector would
equal terminal projection of the full vector, that peak retained support
would be bounded by 2^(m+1) in both order branches, that β=1 tail compression
would preserve the entire reduced vector, and that expectations would match
classical enumeration of the actual input. Improved sizes/times were a
measured hypothesis, not part of the proof.

Two controls distinguish this from an accidental or zero-expectation match:

- `CNOT(0,1); CNOT(0,1)` is the identity. For observable Z₁ and |+⟩ on
  control 0, premature pruning after just one reverse gate gives zero,
  whereas correctly delayed contraction returns Z₁.
- Keeping only α+1 blocks in the β>1 arithmetic control must change the
  reduced operator at some width. It does, even when the scalar expectation
  is unchanged.

## Baseline sweep

`experiment_control_trace` passed **67/67 checks**, with no capped runs.
For each of (N,a)=(7,6), (7,3), (5,2), widths t=2,3,4,5,6 were run with
and without contraction. Every reduced coefficient agreed exactly (maximum
error 0) in all 15 rows. The nine rows with t≤4 additionally agreed exactly
with independent Walsh transforms. Width is the only varied parameter within
each series. Across bases, loaded constants and gate counts can also change.

All these instances have m=13 work qubits, so the predicted peak bound is
16,384 terms. The table reports retained terms before each contraction,
including temporary within-block support. It is not a measurement of bytes.

| N | a | β | t | Full peak | Reduced peak | Full seconds | Reduced seconds |
|---|---|---|---|---:|---:|---:|---:|
| 7 | 6 | 1 | 2 | 24,369 | 12,186 | 1.456 | 0.751 |
| 7 | 6 | 1 | 6 | 24,369 | 12,186 | 5.228 | 2.568 |
| 7 | 3 | 3 | 2 | 24,412 | 12,203 | 1.510 | 0.808 |
| 7 | 3 | 3 | 3 | 48,855 | 12,217 | 3.937 | 1.254 |
| 7 | 3 | 3 | 4 | 98,018 | 12,264 | 12.009 | 2.220 |
| 7 | 3 | 3 | 5 | 196,060 | 12,275 | 29.788 | 3.203 |
| 7 | 3 | 3 | 6 | 392,458 | 12,284 | 70.036 | 4.224 |
| 5 | 2 | 1 | 2 | 24,386 | 12,193 | 1.427 | 0.740 |
| 5 | 2 | 1 | 3 | 48,972 | 12,255 | 5.180 | 1.656 |
| 5 | 2 | 1 | 6 | 48,972 | 12,255 | 8.056 | 3.110 |

The largest odd-order baseline has about **31.95× fewer peak retained terms**
and a **16.58× shorter measured propagation time** with contraction. These
are single CPU wall-time measurements on this machine (CPython 3.12.10,
NumPy 2.4.6), not repeated performance benchmarks. Timing excludes circuit
construction but includes the contraction-schedule scan.

The β=1 compressed prefixes matched every reduced coefficient at all tested
widths. At t=6 they took about 0.77 s for N=7,a=6 and 1.72 s for N=5,a=2.
Applying this compression to N=7,a=3 at t=6 gave maximum coefficient error
0.0856475830078125, as the negative control required.

## Extended check, separate from the baseline

At t=16 (29 total qubits), the saved extended experiment passed **6/6 checks**:

| N | a | Reduced peak | Reduced final | Seconds | Two-block compression error |
|---|---|---:|---:|---:|---:|
| 7 | 6 | 12,186 | 3,883 | 6.954 | 0 |
| 7 | 3 | 12,286 | 4,064 | 13.981 | 0.14261846244335175 |

Expectations matched enumeration of all 65,536 legitimate exponent inputs.
The β=1 row also matched its compressed prefix coefficient-for-coefficient;
the β>1 compression still failed. No full uncontracted operator or Walsh
vector was constructed at t=16, so this is not an additional full-vector
baseline comparison.

## Regression checks and reproducibility

The implementation is the opt-in `trace_plus` path in the existing
`perm_pps.py`, not a second propagator. `test_perm_pps.py` includes reused and
unused controls, invalid-qubit validation, 24 seeded random classical circuits
with independent Walsh-vector checks (maximum error 0), and independent
state-vector expectations (maximum error 5.33×10^−15). Thirteen random cases
had nonzero, non-full reduced support, guarding against vacuous agreement.

```bash
uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_control_trace --max-width 6
uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_control_trace --extended-width 16
uv run --no-project --python 3.12 --with 'numpy<2.5' python test_perm_pps.py
```

Generated raw rows live in ignored `out/control_trace.json` and
`out/control_trace_extended.json`. The extended command does not rerun the
expensive full baseline. The Python 3.12 invocation avoids the documented
intermittent long-job failure in the project's pinned Python 3.14 environment.

## Interpretation

This resolves 12c in a way its original framing missed: early contraction
works in the odd-order branch too. The relevant constraint is the number of
simultaneously live variables in this contraction order, not whether β=1.
The β=1 involution structure gives the additional idempotent-tail shortcut.

Variable elimination is standard, and the broader relationship between
contraction order and simulation complexity is established in
[Markov and Shi, *Simulating quantum computation by contracting tensor networks*](https://arxiv.org/abs/quant-ph/0511069).
That is context, not a claim that our specialization is a new general method.

The measured task is **pre-inverse-QFT work-register Z**, with independent
|+⟩ exponent inputs. Full Shor output sampling introduces different observables
and dependencies. Nor is this competitive evidence against direct classical
sampling for an expectation already expressed as an average of a classical
function. Its impact is on the interpretation and implementation of PPS's
chosen operator representation; it motivates carrying this distinction into
the inverse-QFT investigation.
