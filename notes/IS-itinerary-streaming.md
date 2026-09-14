---
code: IS
date: 2026-09-12
title: "Itinerary replay removes the conditional interval table with a measured allocation saving"
outcome: exact all-mask conditional evaluator validated against blind native gates; frozen task-allocation thresholds passed against both current and packed C89, with no process-RSS saving established
claims: [C94]
todo: [50]
---
# IS — Streaming exact conditional Walsh queries

C94 owns the proof, normalization, charged bounds and scope. TODO50 owns
the remaining mathematical and engineering questions. This is progress
toward the user's broader simulation objective, not its completion.

## Derivation, takeover and independent review

The user explicitly transferred coordination after Claude exhausted its
tokens. Root preserved the dirty workspace and used the authorized board
swarm, with Sol for bounded reference/provenance work and fresh Astra
contexts for substantive reviews. Root alone edited canonical sources.

The all-mask replay and local cap formula had already been derived and
independently audited. Sol archived their immutable messages and corrected
the benchmark contract into `Saba4c3b08a6841d5`, task `Tf47390e292054b28`,
attempt `A050923a6a5524fc3`. The authoritative frozen contract is
`out/agent-board/workers/A050923a6a5524fc3/reconciled_contract.md`.
Before candidate execution, root added online coalescing of adjacent equal
final shifts, recorded as `M28b9da97c6e94d02`. It changes the contraction
work without introducing a global table or assuming fixed output flags.

The implementation is `experiments/experiment_interval_stream.py`, SHA256
`b5e33f4d4406a24cdc4da24c1f534d8783403ce631bcb3dfeb6548853ff9d52a`.
It reuses C89's actual `macro_image` and four-state `prefix_walsh`; no
second propagator was written. Draft implementation task `T966e4c2e5b6143d6`,
submission `S68b1a640104e4253`, received independent code review
`V8f303ea2448a4369` before candidate execution.

The completed scientific submission is `S126923357c974d5a`, digest
`ae53c249790590350a1fa59affb7b206847bbfb720f7e86176894c7a2f59d72e`,
under task `T67102e9799a6441b`, attempt `A401ea237ec804a96`. A fresh Astra
scientific referee accepted it as `Ve80ccfd3cb5a43c8`; root's shadow audit
concurred. Its full derivation is the frozen attempt's `derivation_v1.md`.
The measured source is unchanged from the reviewed implementation.

## Native correctness and meaningful faults

Core ran first and exited zero; the full log is
`out/stream_replay_core_initial.log`. A fresh Sol worker received only the
actual native circuit contract and permitted gate/reference sources,
without the candidate derivation, expected witnesses or implementation.
Its accepted submission `Sd3de4ca7b85e451b`, task `Tf6fb639605854835`, was
frozen before candidate execution. Reference source/report/log live in
`out/agent-board/workers/A29335e4e20864fc8/` as `reference_v2.py`,
`report_v2.json` and `run_v2.log`.

This reference uses the actual native builder and existing logical gate
interpreter. It exercises X, CNOT and Toffoli, retains all outside wires
and verifies unchanged spectators. It independently constructs conditional
maps and Python-integer signs; it does not use C89's word identity or DP.
A shared builder/logical-interpreter bug remains a limitation. Three
fixed-family mutations detect omitted incoming carry, reversed complete
macro chronology and omission of an actual flag-setting CNOT.

Candidate evidence is in `out/agent-board/workers/A401ea237ec804a96/`:
`correctness_v1.py`, `correctness_report_v1.json`, `correctness_run_v1.log`.
Prediction `M500ba3cc1cdb4a16` preceded run `R73e88559734448c8`, whose actual
exit was zero. All six checks passed without harness warnings. They cover
864 conditional maps, 65,536 exhaustive n=2 all-mask numerators and 768
fixed n=3 mask numerators. Every emitted piece matches C89's maximal table,
and every raw itinerary count obeys the proved bound. Zero-depth identity,
ascending and completely reversed native orders are included.

Both pre-derived candidate boundary mutants were detected locally:

| Fault | Fixed native witness | Reference versus mutant |
|---|---|---|
| Replace the retained minimum by the latest cap | n=2, t=6, h=u=1, x=3, q=2, lo=5 | Correct [5,6), wrong [5,7), shift -2; gate images 3,12 give local flag sum 0 versus 2 |
| Omit local output-wrap cuts | Same outside values, second constant c=2, local z=1 | Correct cap 1, wrong cap 2; gate images 15,8 versus extrapolated image 16 outside the domain |

The second native macro is recovered from two gate-prefix maps as
F2 composed with inverse(F1). A full-domain zero coefficient was not
substituted for these local sensitivity checks.

## Matched task-allocation and runtime campaign

The frozen contract fixes n=128,m=129, odd Mersenne modulus, a=1, packed
scratch/control words and four arbitrary input/output mask pairs. Only
the nested prefix length q varies. The baseline is the unchanged current
C89 compiler, compiled once for four queries. Both methods receive the
same immutable lazy constant/control views; neither expands a q-word
list. Both return four unnormalized integers, with coefficient denominator
2^130. At q=127 the returned tuple is (-4,-268,-324,36).

`benchmark_v1.py`, `benchmark_report_v1.json`, `benchmark_children_v1.jsonl`
and `benchmark_run_v1.log` are in the candidate attempt directory. Prediction
`Md0f36f855c2c4816` preceded run `Rb9a531d7b58d42ce`, which exited zero.
All 49 child processes completed, all seven checks passed, and all four
numerators agree between methods and across metric processes.

Allocation, RSS and timing use separate fresh children. Allocation tracing
starts after common imports, before constructing packed task inputs,
masks, views, method workspace and retained outputs. Reference data and
JSON serialization are outside measured children. The table includes its
current construction, merging and validation costs. Timing is untraced,
warmed, and includes full input setup and four-output evaluation in each
of three batches, with alternating method order across q. All batches used
one invocation and exceeded clock resolution by ample margins.

| q | Table peak traced bytes | Stream peak traced bytes | Table median seconds | Stream median seconds | Final pieces | Raw stream pieces |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7,636 | 3,632 | 0.011171 | 0.011049 | 11 | 16 |
| 3 | 15,016 | 3,712 | 0.025591 | 0.026020 | 23 | 34 |
| 7 | 30,540 | 3,712 | 0.061443 | 0.063643 | 53 | 82 |
| 15 | 59,664 | 3,596 | 0.133995 | 0.138942 | 109 | 170 |
| 31 | 116,780 | 3,876 | 0.272334 | 0.284707 | 221 | 346 |
| 63 | 233,460 | 3,928 | 0.548083 | 0.606645 | 445 | 698 |
| 127 | 472,828 | 3,884 | 1.127865 | 1.334633 | 893 | 1,402 |

Derived directly from the fixed packed input, enabled/disabled counts are
(1,0), (1,2), (4,3), (8,7), (16,15), (32,31), (64,63) in this q order.
Both kinds are present from q=3. The script checks that fact; the explicit
counts here are input bookkeeping, not an additional experimental sweep.

All four preregistered engineering thresholds passed: stream/table traced
peak at q=127 is 0.0082144 (about 122-fold smaller); stream growth from
q=31 to 127 is 1.002064 while table growth is 4.048878; the median runtime
ratio is 1.183327. Thus the observed task-allocation reduction costs about
18% runtime at the largest fixture. Tracer bookkeeping is separate:
228,112 bytes for the table and 8,016 for streaming at that fixture.
Timing ranges are 1.123398–1.138057 and 1.323584–1.373223 seconds,
respectively; this is one bounded campaign, not repeated host calibration.

Absolute maximum RSS was 242,432 KiB for both methods and the noop at
q=127. A few other table children reported 242,688 KiB. No repeatable
process-memory reduction is established, and no high-water subtraction is
reported as auxiliary storage. These process values include the common
runtime/import footprint; traced task allocation is a separate measurement.

The campaign used 75.425 seconds wall time and 78.207 seconds aggregate
reported CPU time; the longest child wall time was 11.913 seconds. The
hard caps were 600 seconds aggregate CPU/wall and 60 seconds per child,
with conservative reserves. The host was Linux x86_64, CPython 3.14.0;
BLAS thread counts were fixed at one. No GPU, dense state or long native
permutation job was used.

The Q0 congruence B=4q mod8 from the archived boundary audit guarantees
nonzero results for these odd prefixes; all 42 non-noop results reject a
zero-output shortcut. That congruence is not the mutation-sensitivity
evidence. The native local controls above supply that evidence.

## Interpretation and provenance limits

The result is a measured allocation improvement for exact conditional
scalar queries over both current shared-lazy C89 and the packed comparator
recorded below. The proof also removes the stored global interval table.
Other packed/in-place algorithms remain possible, so neither numerical
factor is an optimality or literature-wide comparison. C93 remains the
stronger special fully averaged time comparator. No full outside-average
improvement, m-scaling experiment, process-RSS saving or general CNOT/PPS
memory breakthrough follows. TODO50 alone owns the next discriminator.

The reference worker retained its original v1 source/report/log; v2 only
replaced dynamic harness control IDs by literal resolutions, with unchanged
science. Root's first linter invocation used an invalid combined-file
`--exit-code` command and exited 2 before evidence inspection. Its usage
log `correctness_lint_v1.log` is preserved; new separate source/report and
terminal-log commands passed. No science was rerun for that CLI correction.
Initial frozen submission lint has zero errors and four expected information
messages for the linter logs' lack of scientific harness verdicts.

The old unsubmitted Claude residual-rank attempt was cancelled after the
user's takeover instruction and a process check; its artifacts remain.
The earlier Claude survey/reference submissions remain separate historical
work, not accepted support for this checkpoint. The new accepted blind
reference avoids their unresolved evidence issues.

Canonical integration task `T7e33f41f1e3246f7` owns the final documentation
checks and review. Core plus the bounded candidate/reference experiments
were the proportional science validation; the other suites were not
rerun. Production helpers, manuscripts, host settings, commits and
publication were outside this change. Frozen evidence files must not be
overwritten when reproducing; use new attempt/report paths.

## Packed comparator follow-up

B1 is now implemented as `experiments/experiment_packed_intervals.py`, SHA256
`b9f697b224189eb64132e3aa845d15e073ecb38fef57d181791cd15bc7a05f50`.
Each flat record stores its unsigned high endpoint and signed shift; the low
endpoint is the preceding high. Fields occupy ceil((m+2)/8) bytes so the
flag displacement and terminal endpoint fit exactly. Composition retains
current and next byte buffers, using the same local C89 macro pieces and
equal-shift merging. A separate packed image buffer, heapsorted in place,
validates image coverage. All buffers, conversion and validation are charged;
there is no global Python tuple/int table in this comparator.

Core passed before new scientific code (`out/packed_interval_core_initial.log`).
The frozen proof/contract and complete evidence are in
`out/agent-board/workers/A67a3a4f35d384855/`: `contract_v1.md`,
`correctness_v1.py`, `correctness_report_v1.json`, `correctness_run_v1.log`,
`benchmark_v1.py`, `benchmark_report_v1.json`, `benchmark_children_v1.jsonl`
and `benchmark_run_v1.log`.

Prediction `M0a111a377da848e0` preceded run `Rd59b270f69844f77`. The accepted
blind reference was reused without execution. All 864 maps, 66,304 supplied
numerators and C89 maximal pieces match. Both new faults are detected on
the fixed n=2,q=1,ascending,t=h=u=x=0 all-mask table: truncating shifts modulo
M gives -6 instead of 0 at masks (0,8); merging unequal shifts gives 16
instead of 0 at (1,1). All six checks passed. Unchanged streaming-native
and boundary checks were not repeated.

Prediction `M24973b77036344a4` preceded benchmark run `Rc9aafb6261ce469d`.
All 49 fresh metric children and seven checks passed, using the same fixed
inputs, q values, four outputs and metric separation as before. The same
four engineering thresholds were frozen with packed B1 as the comparator.
The streaming source is unchanged; a common measurement wrapper serves both
methods. These are paired fresh results, not cross-campaign differences.

| q | Packed peak traced bytes | Stream peak traced bytes | Packed median seconds | Stream median seconds |
|---:|---:|---:|---:|---:|
| 1 | 7,748 | 3,664 | 0.011493 | 0.011311 |
| 3 | 10,158 | 3,744 | 0.026660 | 0.026745 |
| 7 | 12,566 | 3,744 | 0.065243 | 0.065316 |
| 15 | 16,036 | 3,820 | 0.138557 | 0.138594 |
| 31 | 23,876 | 3,908 | 0.286189 | 0.296642 |
| 63 | 39,906 | 3,960 | 0.581285 | 0.620559 |
| 127 | 72,096 | 3,916 | 1.187158 | 1.363179 |

At q=127, stream/packed peak is 0.0543165, about 18.4-fold smaller. Stream
growth from q=31 is 1.002047, packed growth is 3.019601, and the median
runtime ratio is 1.148271. Thus the saving survives this stronger comparator
at about 15% additional runtime. Timing ranges are 1.183667–1.201958 seconds
for packed and 1.359189–1.369010 for stream. The four integer results match
the initial campaign exactly.

The final packed payload is 30,362 bytes at that fixture; the largest combined
composition payload is 60,452 bytes. The validator adds a 30,362-byte payload
after composition. These encoded-data diagnostics are separate from measured
peak allocation, which also includes bytearray spare capacity and Python
objects. Tracer bookkeeping is 15,088 bytes for packed and 8,080 for stream.
Absolute RSS is 242,688 KiB for both methods and every noop in this campaign;
no process-memory saving is established.

The campaign used 77.013 seconds wall and 79.801 aggregate reported CPU;
the longest child took 12.175 seconds wall (rounded upward). All CPU/wall
caps held; there were no timeouts or stderr diagnostics. Write experiment
task `T8db9c5b8e42349a9` froze submission `Sae731de7c5364537`, digest
`3cc2a397e6622b4cebe2e0e09637f2c6b98c8caf29d254f937d65ed96fe509c9`.
Fresh Astra code/science/evidence review `Vb55fcb06b9364e25` accepted it;
root's shadow audit concurred. Lint has zero errors, four information entries
for actual lint logs and one novelty advisory on explicitly negative scope
wording. The comparison validates one two-buffer packed implementation,
not every in-place representation. Integration task `Tb9fac3a1ea8e4eac`
owns this follow-up's canonical review and documentation checks. TODO50
alone owns further work; the general simulation objective remains open.
