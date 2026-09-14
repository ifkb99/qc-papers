---
code: IX
date: 2026-09-12
title: "Actual controlled-add prefixes collapse to conditional interval exchanges"
outcome: exact conditional representation and Walsh contraction proved and independently checked; the shared-variable average remains open
claims: [C89]
todo: [50]
---
# IX — Conditional interval contraction of actual controlled arithmetic

C89 owns the mathematical statements, input contract, source citations and
algorithmic limits. This is progress within TODO50, which remains open.
The full-space scalar includes dirty scratch, incoming carry and flag;
the clean logical arithmetic is a separate baseline. No general simulation
breakthrough or measured memory advantage is established. The user's goal
remains active, and TODO42 remains separate.

## Discovery and proof review

The initial question was whether shared-state controlled arithmetic leaves
the known planar carry contraction. Reading the actual `cc_add_mod` schedule
first exposed a stronger global baseline. Root and Astra independently
derived the fixed-scratch accumulator/flag interval representation, its
additive composition bound and exact range-Walsh contraction. The flag's
chronological dependence survives, but does not require retaining a tree of
its possible histories after the outside inputs are fixed.

Sol audited primary interval-exchange, weighted-automaton and decision-diagram
sources. Root read Novak's definition and preimage argument, Kiefer's
minimization theorem, and Bryant's composition/counting algorithm bodies.
These supply known frameworks, not a cheap construction for our remaining
outer average. The local shared-control matchgate tensor was not classified:
the conditional global reduction was established first. C89 records exactly
where it stops, including the controlled swaps that change register roles.

Root and Astra independently derived the nonzero single-macro probes before
the main run. A separate worker's clean active two-macro fixture then became
identity on every accumulator and flag input. Root proposed the complementary-
shift inverse rule and Astra checked its algebra independently. That proof
was added after observing the identity; it is an explanation, not a prior
prediction. It also identifies a second, fully dirty scratch sector.

## Main exact experiment

Source: `experiments/experiment_controlled_intervals.py`. Reproduce with:

```bash
timeout 60s uv run python -u -X faulthandler -m experiments.experiment_controlled_intervals
```

Prediction `M4f273739d7b74f42` preceded run `R65537c4b5abc4af9` and execution.
The fixed fixture is N=5, a=2, one exponent control, four accumulator bits
and fourteen total qubits. Only the native prefix length changes. Each
row reuses the same accumulator, scratch, carry, exponent control and flag.
Report: `out/controlled_intervals_report.json`; raw log:
`out/controlled_intervals_main.log`. The process terminated with exit zero;
all eight harness checks pass without warnings.

| prefix length | native constants | logical gates | maximum conditional pieces | two full-space coefficients |
|---:|---|---:|---:|---|
| 1 | 2 | 138 | 11 | 1/4, -1/8 |
| 2 | 2, 4 | 276 | 18 | 0, 0 |
| 3 | 2, 4, 3 | 420 | 24 | 1/512, -1/512 |

Each prefix matches actual `walsh.classical_permutation` gate replay on all
16,384 full-space labels: 49,152 map comparisons with no errors. Every one
of its 512 outside-input assignments is explicitly visited, and the full
output label comparison includes all restored wires. The finite counts do
not prove C89's general interval bound; the separate construction proof does.

The four-state digit routine also matches 69,632 independent direct signed
prefix sums, exhausting four-bit translations, input/output masks and range
thresholds. Both selected full-space coefficients agree between conditional
interval averaging and actual gate replay in every row. Only the first row's
specific nonzero values were derived in advance; the other rows test exact
agreement, not a prediction of those particular values. A selected coefficient
vanishing at one prefix can revive later, as this fixture illustrates.

Omitting the incoming carry annihilates both nonzero first-row probes.
Forcing the enable conjunction to zero annihilates the negative probe.
Reordering the shared-state pair changes the packed output from 11011 to
15107 at input 11008, exactly as derived before execution. All three wrong
references therefore fail equality as required. C89 owns the corresponding
symbolic mask definitions and chronological witness.

The reported 0.476-second harness time is descriptive, not a repeated
performance estimate. Allocated bytes and native RSS were not measured.
The truth tables are independent correctness references, not a claimed
memory-efficient implementation of the remaining outside average. The new
helper is experimental and accepts supplied macro parameters; it does not
recognize arbitrary circuits or implement a quantum output sampler.

## Independent actual-gate verification

The accepted verifier remains in its worker directory. Reproduce with:

```bash
PYTHONPATH=. timeout 60s uv run python -u -X faulthandler out/agent-board/workers/Aaab2d7d440a84f76/experiment_cmult_prefix_map.py
```

Its frozen N=3, a=2 fixture has one exponent control, three accumulator bits
and eleven total qubits. The native prefixes have constants [2] and [2,1],
with 108 and 216 logical operations. Actual gate replay agrees with the
worker's separately implemented C86 dirty word recurrence on all 4,096
labels, including restored outside wires. The following wrong-reference
witnesses were calculated before execution:

| defect | input | actual output | wrong output |
|---|---:|---:|---:|
| omit incoming carry | 1344 | 1347 | 1346 |
| fix input flag to zero | 1600 | 1095 | 1090 |
| make disabled dirty macro identity | 8 | 9 | 8 |
| reset shared flag before macro two | 1728 | 1728 | 1221 |

For six frozen outside-input assignments, direct sums over the actual gate
images agree with the root interval helper on all 256 input/output Walsh
mask pairs: 1,536 exact rational coefficient comparisons, no errors.

| prefix | fixed (t,h,u,x) assignments | interval counts |
|---:|---|---|
| 1 | (0,0,1,1), (0,1,1,1), (1,0,0,0) | 10, 9, 8 |
| 2 | (0,0,1,3), (0,1,1,3), (1,0,0,0) | 1, 15, 9 |

The direct side uses actual logical X/CNOT/Toffoli gate images; it does not
independently test the stored Clifford+T decompositions as unitaries. The
conditional side imports the root helper, so independence is supplied by
the direct gate sums and separate word recurrence, not two interval codes.

Prediction `M6dbbbf160ed04d58` preceded the original successful run
`Rcbdcae9f6fc44205`. After the resource-reporting correction below,
prediction `M6ac68b963d234b84` preceded revised run `R5e5d3166db3f473e`.
Both terminal exits were zero; the identical frozen science again passes
all eight checks without warnings. No larger fixture was added. The worker
directory contains `cmult_prefix_map_report.json`, `cmult_prefix_map.log`,
`findings.md` and `frozen_pilot_proposal.md` beside the source.

## Preserved review corrections and provenance

The source audit's first submission `S0dad4fe26afe4bd0` had an incorrect
chronological label for an N=3 witness. Review `V883d1e091807411d` requested
a correction. The revised artifact explicitly distinguishes constants
1 then 2, which end at (2,0), from 2 then 1, which end at (2,1), starting
at t=0,h=1,u=1,x=3,b=f=0. The original artifact remains archived. This is
a correction to the audit prose, not a scientific run failure.

The verifier's initial submission `S923dc2e00ba247db` claimed a four-array,
64 KiB cap without accounting for live arrays during replacement or NumPy
comparison temporaries. Review `Vcd5a6356215c47aa` rejected that resource
statement while accepting the exact calculations as sound on inspection.
The revision releases per-row arrays and reports only the two explicitly
retained actual maps. Total allocated bytes remain unmeasured. Original
source, report, log, proposal and findings are retained with `revision1`
suffixes. No peak-memory claim follows from the corrected evidence.

The first mathematical submission `S68a5a9c9d7e94e4c` was sound. Review
`Vcdd894f21077469e` returned it to active solely to add the post-observation
inverse supplement to frozen evidence. The final submission preserves both
proofs. Root inspected the full actual artifacts, logs and reports and
verified archived hashes/current identities. Astra independently reviewed
the root implementation's full archived source, report and validation logs.

| board task | accepted submission | accepting review |
|---|---|---|
| T65b4e1e9ecc04d56, Astra proof | Saee62f29bbce4278 | V589c6192b2c2491a |
| T9f85a1eeb33247dd, Sol sources | S46df079055c0400c | Vf7e6868811cd47ca |
| T5ac4d70bcca54dd2, Sol verifier | S522e9400faa84f64 | V5b71bb55aa1c42ef |
| T37668231982f4bc3, root implementation | S861d0324f3ae47f4 | Vb592d38789414f45 |

Astra's artifacts are `shared_macro_reduction.md` and
`complementary_shift_supplement.md` under
`out/agent-board/workers/A6d2daa401a4c4014/`. The corrected source audit is
`out/agent-board/workers/A223202275ce943ab/controlled_macro_source_audit.md`.
Workers made no canonical edits; root remains the sole source writer.

## Validation scope

Core ran first and after the new experiment. Both observed terminal exits
were zero, and both logs end in ALL TESTS PASSED:
`out/controlled_carry_core_initial.log` and
`out/controlled_carry_core_final.log`. The main and independent bounded
experiments supply the affected science checks. Production helpers were
unchanged; the other eight science suites were not rerun. These runs had
no timeout or native crash. TODO34's prior native failure is unresolved.

Documentation indexing exits zero (`out/controlled_intervals_reindex.log`),
and all ten documentation checks pass
(`out/controlled_intervals_docs_check.log`). All four bounded board tasks
are accepted and closed. The root task's first closure request listed the
coordinator's claim/note integration outside its declared source write path;
the board rejected that file list. The corrected closure records only the
authorized experimental source. This bookkeeping correction did not change
the accepted source or execute science again.

No manuscript, dependency, default or host setting was changed, and nothing
was committed or published. Preserve the dirty worktree.
