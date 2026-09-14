---
code: SB
date: 2026-09-12
title: "Mask-free blocks in the exact range-Walsh contraction"
outcome: sparse-mask block contraction proved and implemented with independent finite checks; runtime and process-memory effects remain unmeasured
claims: [C96]
todo: [50]
---
# SB — Contracting unsigned bit gaps without losing joint carries

C96 owns the exact block formula, state invariant and charged cost bounds.
TODO50 owns the next performance discriminator. This work supplies an
optional contraction kernel; it does not revise the completed IS benchmarks.

Root noticed that a mask-free bit block has two integer thresholds, for
addition carry and threshold-subtraction borrow. Their intersection counts
can replace its individual signed-bit steps. An Astra author derived the
joint formula and hand fault fixtures before implementation; a fresh Astra
referee accepted it. Root then implemented the optional module and ran the
bounded candidate test, followed by a second fresh Astra code/evidence review.

The scoped primary-source check found established carry automata and
bit-parallel correlation algorithms in
[Wallén, Section 3.3](https://research.ics.aalto.fi/publications/bibdb/HUT-TCS-A84.pdf).
That section uses an n-bit register cost model for its stated averaged
carry object. C96 is a block summation of C89's fixed-shift, thresholded
DP; this search does not establish priority or the fastest competing
thresholded evaluator. No novelty or optimality claim is made.

## Frozen proof and candidate evidence

The accepted pure derivation is task T7c9d42e5395646f8, attempt
A8cf5cf5d38b24f23, submission S023deae7d4044fda, digest
3ac42234bff0be6b8ac9c54b4f9be83f463e6d614d8ce62d272e57de62c32d83.
Fresh review V983b541304cc4183 accepted it; closure Cbb4069c25d6948da
preceded candidate claim. Its `derivation_v1.md`, `contract_v1.json`,
hash record and lint log are in that attempt directory. The derivation
SHA256 is b0db9e3fb72d3458cb235b988cbf688e11325305fbf6d95b98e9460702decc1b.
No scientific execution occurred in the proof task.

Candidate task T6c486b9da006455c, attempt Aa965875ab97c4eed, binds that
accepted proof and accepted blind native reference Sd3de4ca7b85e451b.
The new `experiments/experiment_sparse_walsh.py` SHA256 is
658d77c670ec6419bf23dcf0f393842cde615829387b54cfaa8f50e63e7ce1e9.
The three original interval modules retain their accepted hashes.

The exact run is

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout 180s uv run python -u -X faulthandler out/agent-board/workers/Aa965875ab97c4eed/correctness_v1.py

M0b30e9ed495f400c declared the contract, and M55d219b44e884263 froze
the hashes/predictions before run Rbd06bddc95694cd7. Core had already
passed first in the same continuing scientific session; its log is
`out/packed_interval_core_initial.log`. Actual candidate exit was zero.
All evidence is frozen in submission S2d68fa78a8de4209, digest
7512f09f95b4b7f39b3071e84f72e9681e5d6be91ef2e96289ee123d894c53cb.
Fresh review Vb2af416848684bbc accepted it; closure C27200bb33c1345ba.

The attempt owns `correctness_v1.py`, `correctness_report_v1.json`,
`correctness_run_v1.log`, source/run lint logs and `provenance_v1.md`.
Report SHA256 is 2d34c10a7e2cbc7cf9387b393997401a32e471aa2e2d00ff8d602dc1e262ed3c.

## What the finite check establishes

All 12 harness checks passed:

* Direct signed sums matched all 74,586 prefixes for bits 0 through 4,
  every mask pair, threshold and translation residue represented by
  a centered signed shift.
* Direct block-word arithmetic matched 5,456 joint transition rows
  for block widths 1 through 5, every local shift/limit and all incoming
  carry/borrow states.
* The unchanged C89 bit DP matched 84 deterministic wider cases through
  4,097 bits, including sparse/dense masks, empty/full limits and signed
  shifts spanning multiple wraps. Checked transfer bounds held.
* Each of the packed and streamed interval paths matched all 66,304
  frozen native-gate numerators using the new contraction.
* All 28 invalid-input fixtures rejected before shortcuts.

The five fixed faults were rejected with their predicted wrong results:

| Fault | Arguments (n,L,s,A,C) | Exact | Faulty |
|---|---|---:|---:|
| Independent carry/borrow marginals | (3,6,2,0,4) | -2 | 0 |
| Inclusive borrow threshold | (2,2,0,0,0) | 2 | 3 |
| Strict carry threshold | (3,3,2,0,4) | 1 | 3 |
| Borrow-only full-limit terminal | (2,4,-3,1,1) | -4 | 0 |
| Drop final addition carry-one paths | (2,4,-3,1,1) | -4 | -3 |

Direct sums and block-word arithmetic are elementary independent references
at this layer, so no new blind-reference task or gate interpreter was
needed. The wider recurrence shares C89's four-state model. The earlier
native reference retains its disclosed interpreter/builder limits; its
frozen report was reused, not rerun. Both interval representations still
share accepted local macro arithmetic. The mutant driver deliberately
shares the unaffected signed-bit step to isolate the named fault.

Source/report lint has zero errors and one informational dynamic-control-ID
finding; all five controls are explicitly resolved in the report. Run-log
lint is clean. Submission lint's additional informational entries are the
mechanical lint logs. An initial run.start API request used prose instead
of the prediction message ID and was rejected before any run or execution;
the corrected request and provenance record preserve that administrative
error. There was no failed scientific run or discarded fixture.

The prefix bound is implemented separately from the surrounding conditional
width guards. No allocation, RSS or runtime measurement was performed here.
Both packed and streamed competitors can use the same improvement; no
full outside-average or broader quantum-simulation advantage follows.
Documentation checks and the final integration review are attached to
task T4f60b8fc8ebd4c18, attempt A7a996e99e6f044dc.
