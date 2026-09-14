---
code: CF
date: 2026-09-12
title: "Lazy CNOT coordinates save allocation, but a stronger physical-output baseline wins"
outcome: bounded implementation improvement; no breakthrough established
claims: [C82]
todo: [43]
---
# CF — Direct CNOT memory investigation

The user requested continued research specifically on CNOT memory, using a
model-selected swarm and the local message board. The coordinator implemented
the opt-in frame in the existing permutation engine. Astra audited mathematics
and primary literature; Sol supplied an allocation/baseline draft; Luna
supplied independent small-law checks. Only the coordinator changed canonical
source. C82 owns the mathematical/API/resource statement and prior-art limits.

This goal turn is **progress**: it adds tested implementation, a precise adverse
baseline comparison and new evidence. The previous scientific checkpoint was
also progress (C81/NS); no running job was assumed from handoff text. The initial
board resume had no active assignments. No breakthrough or completion of the
broader goal is claimed.

## Authoritative coordinator reproductions

`experiments/experiment_cnot_frame.py` and
`out/cnot_frame_laws_main.json` / `.log` record the full-coefficient comparison
on 72 fixed-seed mixed circuits with n from 2 through 7, plus explicit boundary
witnesses. This is a correctness corpus, not a scaling sweep. Both coordinate
modes agree with independent classical-replay/Walsh coefficients. Tests cover
arbitrary Z masks, positive thresholds, nonempty weight filtering, physical
trace constraints, cap ordering and Mapping access. The deliberately broken
internal-key view disagrees at physical keys 6 and 7. All three harness checks
pass, with no warnings.

The truncation reference filters the final Walsh vector only in a fixture
where the first reverse CNOT leaves the observable unchanged, one Toffoli
performs all branching, and the final reverse CNOT changes none of those keys.
It is not a general reference for sequential truncation. The permanent core
regression separately checks two CNOTs that cancel but lose a term at an
intermediate physical weight cutoff; dense conjugation, inverse offsets,
X signs, first-gate threshold, cap-before-trace and a 65-bit mask are included.

`experiments/experiment_cnot_frame_memory.py` and
`out/cnot_frame_memory_main.json` / `.log` own allocation evidence. Six fixed
disjoint Toffolis build S=4096 terms before a CNOT prefix is encountered in
reverse propagation. Only the prefix length changes: D=0,1,8,32,128 uses
prefixes of one fixed deterministic schedule. Physical results agree between
modes at each D; different D need not represent the same operator.

At D=128 the complete-call **traced Python allocation peaks** are:

| requested result / implementation | peak bytes |
|---|---:|
| unchanged engine, physical dictionary | 982232 |
| framed engine, lazy physical-key Mapping | 662720 |
| framed engine followed by dict(view), retaining the view | 773224 |
| existing-engine nonlinear suffix plus one batched CNOT remap, physical dictionary | 690008 |

The lazy frame reduces this measured peak by about 32.5% against the unchanged
engine. The stronger batched comparator BEATS the frame when a physical
dictionary is requested. That comparator uses the existing nonlinear engine;
its one boundary remap is specialized to this fixture's CNOT prefix with no
intermediate trace/weight filtering. It is not a second nonlinear simulator.
At D=0 the frame is slightly more expensive, an expected adverse case.

Streaming items has 1576 bytes of incremental traced peak after the view
already exists; materialization adds 279064 retained traced bytes and peaks at
385800 incremental bytes. These separate incremental measurements are not
mistaken for the full-call rows above. Values are shared with the backing
dictionary, while the physical table and decoded keys are allocated anew.

Specialized CNOT-only boundary peaks (prebuilt coefficient dictionary outside
the trace) are 587632 bytes for rebuild-per-gate, 802416 for an in-place
key-snapshot/pair-swap candidate, 309760 for batched dictionary materialization,
132120 for full-array packed relabeling, and 70968 with packed temporary arrays
limited to 256 entries. Packed keys/values alone occupy 65536 bytes; they lack
the dictionary lookup and nonlinear merge interfaces. In-place dictionary
swaps were adverse because the table grew during singleton-key moves.

The operation comparison, 2*D*S dictionary insertions versus 2*D frame XORs,
is **derived from the current source**, not an instrumented work counter or a
runtime measurement. The memory experiment includes a nonempty nonlinear
weight/trace sentinel, with the trace reducing 9 terms to 5. Its seven checks
pass, including the materialization control that defeats the one-dictionary
storage claim. There are no harness warnings.

All traces exclude prebuilt Circuit objects, preexisting reference results,
interpreter/allocator arenas and RSS. Shallow layouts itemize tables, keys,
values, frame vectors and instance dictionaries; they are not complete process
memory. The run used Python 3.14.0 as recorded by the report. No timing result
or stable-host inference follows from these short passes.

## Review and failure provenance

The mathematical audit is submission `Sea0a6a0afbd14e91`; it read the actual
body of the relevant papers, including Clifft. The implementation matches its
reviewed source SHA in the sealed artifact. Initial math and technical tasks
were claimed before implementation changed; the board rejected stale-input
submission, so stopped workers were requeued and rebound to the current
source. Original artifacts were retained with explicit provenance.

Luna's first draft had a vacuous all-zero truncation comparison and its board
prediction reversed cap/trace order in prose. The coordinator required the
nonempty, first-CNOT, weight and uncapped-trace witnesses before accepting
`S7cb631e24a144e07`. Its initial script/report/log remain in attempt
`Aa733df9202e84ee4`, including `report_initial.json` and `run_initial.log`.
Main additionally stated why the fixture's final-filter reference is valid.

Sol's first run failed before measurement with an import-path error; the
retained attempt `A77fd9729df54409b` contains the raw log and structured
`cnot_frame_memory_audit_import_failure.json`. The corrected draft passed and
was submitted as `S82bbd87d633548f2` with limitations. Main added complete
propagation-plus-materialization measurement, a full-call batched comparator,
chunked packed scratch and instance-dictionary accounting before reproduction.
The stronger comparator is the material adverse finding, not a failed law.

Coordinator reproduction is `Se1bfa61adc22415e`, run `R42398b5a92454c18`,
accepted after independent archived-evidence review `V541ffa9d559d4ade` by the
Astra auditor. All four bounded board tasks are closed.

## Validation and incomplete gate

The core gate passed before science. After changes,
`uv run python -u -X faulthandler test_core.py` passed including the new frame
regressions; `out/cnot-frame-core.log` is the evidence. Core now also exits
nonzero when its accumulated FAILED list is nonempty, making failures enforce
the gate rather than only print a warning.

The full legacy `test_perm_pps.py` run terminated with exit 139 and an empty
buffered log. Session 17411 was confirmed terminal; the failure is preserved
in `out/cnot-frame-perm-initial-status.json` and its named log. No cause was
inferred, no full sweep was restarted, and this suite is **incomplete**. Its
large rotation-level default comparisons did not contain the new opt-in tests.
TODO34 remains the host-reliability record. The other seven science suites
were not run for this scoped change.

Documentation validation is recorded in `out/cnot_frame_docs.log` after index
regeneration. No manuscripts/abstracts, host settings, commits or publication
were changed. The original dirty worktree is preserved.
