---
code: WE
date: 2026-09-11
title: "A setup-only word envelope tightens merged sampling"
outcome: confirmed
claims: [C59, C60, C71, C72, C73]
todo: [14, 24, 32, 33]
---
# WE — Count the supplied word once

C73 owns the proof, implementation contract and resource limits. Main added
opt-in word support to the existing reverse-vector and merged-prefix adapters;
lower-cost agents supplied independent support, law and comparison harnesses.
Main strengthened their predicates before rerunning. Alphabet mode remains
the default. This is an improved global rejection envelope, not new group
theory or a new numerical certificate.

## Independent support and complete tiny laws

Main's `out/word_envelope_support_20260911T083245217209Z.json` passes 4/4.
The verifier separately evolves actual reached sectors and affine E/O maps,
comparing exact sets at every reverse depth for all sectors at M=1,...,31
and selected sectors at M=101,1009,1000003. It also checks setup counters,
scalar-only returned metadata, input snapshots and pre-growth rejection.
Its 164,184 charged reference updates/maps are below the frozen 500,000 cap;
these are not all CPU instructions. Simultaneously retained gamma/prefix
sets and their temporary updates have a separate preallocation guard.
Production releases its offset sets and never retains these verifier tables.

The 64-insertion, nine-label structural input in the raw report is accepted
by word mode with envelope 226 while alphabet mode rejects its conservative
payload cover. This is an input/setup check, NOT a complete 64-insertion law
or a proof that actual memory necessarily exceeds a cap in alphabet mode.
The genuinely growing powers-of-three alphabet must raise MemoryError from
the word setup's live-offset cap; a ValueError does not satisfy that control.

Main's `out/word_envelope_laws_20260911T083252349572Z.json` passes 4/4.
Ten fixtures in both modes give 20 reverse law rows and 6,170 forced
sector/fine-boundary/output cases, including a genuine nine-insertion circuit.
Every accepted submass and normalization is checked against the existing
independent full-r contraction and sequential instrument. Maximum accepted-law
TV discrepancy is 3.59e-16; success-mass error is at most 5.56e-17. These are
floating diagnostics, not certified error bounds on unseen executions.

There are 2,765 canonical prefix requests, evaluated in both modes and against
the full-r reference: mode discrepancy is exactly zero in this run, and the
largest reference discrepancy is 3.04e-16. The fixtures include odd/even coarse
moduli, fixed points, noninvertible differences, initial/terminal insertions,
zero angles and genuinely nonzero even-modulus rotations. Twenty-four forced
rows have exactly zero terminal norm and correctly return zero proposal,
acceptance, numerator and submass; these are not caught numerical exceptions.

On the unsaturated r=14,b=2 chronological word (0,1,0), the global envelope
drops from seven to six, and the mathematical success probability rises from
1/14 to 1/12 with unchanged accepted law. Gamma zero alone reaches three
sectors, while other gammas reach six: the local-only global-envelope control
fails as required. This is a set-bound counterexample, not a claim that this
specific wrong envelope necessarily produces acceptance greater than one.

Separate UNMEMOIZED returned samples in both modes have actual attempt work
summed independently, and prefix work checked against existing direct-set
bookkeeping. History calls are forbidden. The shared reference guard runs
before dense allocation; its largest nine-insertion operation allowance is
212,362,884 entry touches and simultaneous dense payload allowance 336,384
bytes. These are reference bounds, not production costs or native RSS.

## Comparison checkpoint: long runs aborted, not a completed benchmark

The frozen fixture uses M=2^40-1,b=3,t=16 and FA's separated three-label
cyclic word, pi/4 rotations, initial Rx(pi/4) and alternating subsequent
Rx(pi/7)/Rz(pi/5) blocks. k varies by one from zero through sixteen, with
seeds 7430–7432. Four primary methods use reverse/prefix contraction and
alphabet/word envelopes. Each reverse row has an 8,192-proposal cap and every
row a soft 60-second limit. Main froze an additional 600-second sweep budget
before running; optional legacy prefix anchors at k=0,4,8 require a cheap
preflight. Setup, independent verifier work, every attempt/query and working
complex slots are separate. Word sets are constructed once and forbidden
from being rebuilt during sampling. Instrumentation uses constant-size
running counters, not an O(attempts) list.

Lower-cost initial high-k preflight
`out/word_envelope_comparison_20260911T083708470044Z.json` passes its four
row checks. At seed 7430, word mode reduces reverse attempts from 1,366 to
100; both prefix modes use 81 queries and equal support peaks. This is ONE
seed, not a typical speedup or a completed four-method comparison. The reduced
word bound changes rejection, whereas the prefix contractions remain the same.

Main's serial `out/word_envelope_comparison_main.log` exits 139 with a
segmentation fault after the k=12 groups. Faulthandler reaches sparse_norm's
NumPy vdot call. The harness had not yet written a final JSON, so numerical
rows from that run were lost; the log is preserved. An isolated k=13,
reverse_alphabet, seed-7430 reproduction passes in
`out/word_envelope_comparison_20260911T084200760406Z.json`. The failed run's
exact in-group seed was not logged, so this is a nearby reproduction, not
proof that the identical failing call passed.

Main added per-row checkpointing and explicit start labels, then repeated
the same environment. `out/word_envelope_comparison_checkpointed.log` exits
132 with Illegal instruction during unit_phase's NumPy trigonometry at
k=10, reverse_alphabet, seed 7431. Its
`out/word_envelope_comparison_20260911T084241528152Z.partial.json` preserves
134 completed rows (four preflight, nine legacy anchors, 121 sweep rows),
explicitly marked INCOMPLETE, not a scientific verdict. These different native
fault sites do not establish a NumPy, interpreter or hardware root cause.
Both runs use Python 3.12.10 and NumPy 2.4.6; neither is the documented 3.14
failure assumed to recur. No norms, phases or production loops were changed
to suppress the crashes.

System Python 3.12.3 with the SAME NumPy version passes core and all new tiny
laws, and the four-row high-k diagnostic in
`out/word_envelope_comparison_20260911T084458560301Z.json`. This short preflight
does not prove long-run stability. A brief exact-integer follow-up overlapped
its execution, so its timing is not an authoritative serial comparison either.
Do not pool these environments or cherry-pick completed rows from the aborted
sweeps. Host/runtime evidence is saved in `out/word_envelope_runtime_environment.log`;
TODO 34 owns the diagnostic plan and hardware lead. No system settings changed.

TODO 32 remains open. Do not infer runtime improvement from the envelope ratio
alone. The tiny full-r/static baseline remains inexpensive, single-label input
simplifications remain available, and C71's neighboring-label error bound does
not certify the separated replacement. Failure of that sufficient bound or
the huge static group is not hardness. Long scaling is postponed pending
reliability; the more consequential coordinate-promise follow-up is C74/PC.

## Audit trail and scoped validation

- Early law reports preserve an import of a nonexistent reference helper and
  a verifier division by zero on valid exact-zero branches. Production did
  not throw those verifier errors. Later passing versions still reported an
  underestimated reference payload and checked only summary sample counters;
  main replaced those reports with shared-guard metadata and actual internal
  counter instrumentation, and added a nonzero even-modulus fixture.
- The first support cap control used a word that did not grow enough and
  correctly failed to fail. Subsequent stronger reference versions exhausted
  their operation budgets by rerunning direct recurrence at every depth;
  caching each gamma's depth sequence once removed that redundant work without
  dropping any of the 34 rows. All reports remain under `out/word_envelope_support_*`.
  Earlier passing predicates compared only the E/O representation, not an
  independent reached-set recurrence. Main also charged the comparison maps
  and guarded combined all-gamma retention before the final rerun.
- The comparison code's initial draft charged an extra independent E/O audit
  only to word-mode setup, checked call counts but not summed internal work,
  and had no history-item accounting for the legacy anchor. These were audit findings
  before timing, not measured performance failures.

Core ran first. Logs `out/word_envelope_core.log`,
`out/word_envelope_lab.log`, `out/word_envelope_claims.log` and
`out/word_envelope_lab_no_flint.log` pass. The last explicitly skips optional
verified-arithmetic tests. The existing FA full-prefix/law regression passes
in `out/word_envelope_fa_regression.log` with report
`out/fixed_alphabet_prefixes_20260911T083001429761Z.json`; the MP complete-law
regression passes in `out/word_envelope_mp_regression.log` with report
`out/merged_prefix_sampling_20260911T082938440507Z.json`.
The other six science suites were NOT rerun. No manuscript, abstract workshop
or sampler default was changed; no commit was made.
Documentation gates are recorded in `out/word_envelope_and_coordinates_docs.log`.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_word_envelope_support
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -u -m experiments.experiment_word_envelope_laws
```
