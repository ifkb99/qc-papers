---
code: ME
date: 2026-09-12
title: "Essential matching support survives, but sparse Pfaffian methods explain the memory saving"
outcome: growing reduced grid proved; exact sparse implementation beats current dense allocations on bounded fixtures; known stronger baselines prevent a breakthrough claim
claims: [C88]
todo: [49]
---
# ME — Matching reductions and a charged sparse comparison

C88 owns the proof, input contract and algorithmic limits. This completes
TODO49's bounded reduction discriminator. The closed-form question is not
settled, but no identified new algorithmic advantage survives the audit.
TODO50 owns the redirection; TODO42 remains separate. The broad user goal
remains active.

## Discovery and independent review

Astra derived a two-interface matching criterion and explicit witnesses
showing that the C87 grid survives complete forced/forbidden-edge reduction.
A Sol verifier independently classified matching edges on the fixed tiny
family. The smallest instance collapses to a cycle, but the larger fixtures
retain a connected nonbipartite graph. This refutes a growing forced-edge
explanation of those rows; the all-k proof is separately in C88.

Root proposed using complete all-carry matching pairs as elimination pivots.
Extending every prefix matching by the fixed complement gives a short
nonvanishing proof. Astra checked it independently and derived exact band
and bit bounds. Root implemented ordinary sparse skew Schur elimination,
including a generic pair-min-fill comparison, and Astra reviewed the actual
archived source/report/log. This applies known algebra to the arithmetic
specialization; it is not a new general simulation method.

Sol's source audit identified stronger nested-dissection baselines and
showed why a signed mixture of the same planar matchgate alphabet cannot
escape FKT. Root read the cited primary theorem/algorithm bodies and made
the whole-pair block transfer explicit. C88 owns those citations and their
scope. Standard bipartite tiling formulas do not directly identify the
nonbipartite fixtures, but no absence of a more indirect formula is proved.

## Main exact comparison

Source: `experiments/experiment_carry_matching.py`. Reproduce with:

```bash
uv run --with networkx==3.5 python -u -X faulthandler -m experiments.experiment_carry_matching
```

Prediction `Mad313a64381c4e8f` preceded run `Re05af483a5244dcf` and source
execution. The parent runs nine fresh subprocess cases, each capped at
60 seconds, varying only k=1,2,3 in the fixed C87 query. Report:
`out/carry_matching_report.json`; log: `out/carry_matching_main.log`.
The observed process exit was zero; all six checks pass without warnings.

| k | word width | dense peak bytes | banded peak bytes | pair-min-fill peak bytes |
|---:|---:|---:|---:|---:|
| 1 | 8 | 151357 | 97996 | 100476 |
| 2 | 32 | 632205 | 397980 | 401156 |
| 3 | 72 | 3221237 | 905094 | 908270 |

All three methods give identical exact rational coefficients in every row,
matching the C87/DM values. At the largest fixture, the banded evaluator
uses about 28.1% of the dense evaluator's traced peak. Generic min-fill
comes very close; its small difference is not evidence of general ordering
superiority. Single measured elapsed times at that fixture are 0.6452,
0.0227 and 0.0311 seconds respectively. These are descriptive timings under
allocation tracing, not repeated native performance estimates.

The trace starts before graph construction and includes ordering,
planarity/orientation, exact arithmetic and coefficient extraction.
Imports/process startup and native/process-wide RSS are excluded. The
implementation stores the original sparse tail as well as the active
frontier. Its entry counter rescans rows at each pivot, and is charged in
the measurements. The report's `max_updated_fraction_bits` covers only
nonzero updated matrix entries; it is not a maximum over pivots, products,
initial entries or transient allocations. C88's bit bound is a proof,
not an inference from that statistic.

Natural bandwidths are 6,10,14. Independent endpoint-gap enumeration
confirms the predicted quadratic polynomial coefficients 1,66,435.
The wrong matching-pair certificate is rejected. Deleting the allowed
interface edge (1,11) changes the tiny partition from 80 to 64, so the
second must-fail control is nonvacuous. No main prediction was refuted.

Dense and sparse routes share the graph constructor and orientation helper.
This experiment checks elimination and cost on the accepted C87 scalar;
it is not a fresh independent quantum-word reference. DM owns the prior
carry-sensitive and permutation-sensitive arithmetic checks. Optimized
nested dissection was audited theoretically, not benchmarked. Neither a
best-known memory advantage nor a general CNOT-memory breakthrough follows.

## Independent essential-edge experiment

The accepted source remains in the Sol attempt directory:

```bash
PYTHONPATH=. uv run --with networkx==3.5 python out/agent-board/workers/Ab4b526cf94734ba8/experiment_matching_reduction.py
```

The source header omits the required `PYTHONPATH=.` prefix; use the command
above. Prediction `Ma924bd63523c4f6b` preceded both recorded attempts.
The first launch, `Rc79b1f728245401c`, exited one before measurement because
`experiments` was absent from its import path. Its preserved log is
`out/agent-board/workers/Ab4b526cf94734ba8/matching_reduction_attempt1_import_failure.log`.
No source repair was needed. Correcting the launch environment produced
run `R399b6012f65a4ab4`, exit zero, five passing checks and no warnings.
The same directory owns `matching_reduction_report.json` and
`matching_reduction.log`.

| k | allowed / forbidden / forced edges | residual after extracting forced pairs |
|---:|---|---|
| 1 | 10 / 4 / 2 | 8 vertices, 8 edges, one cycle, prefactor 4 |
| 2 | 84 / 2 / 0 | 60 vertices, 84 edges, connected and nonbipartite |
| 3 | 206 / 0 / 0 | 140 vertices, 206 edges, connected and nonbipartite |

An edge is allowed iff deleting its endpoints leaves a perfect matching;
it is forced iff deleting the edge destroys every perfect matching.
Feasibility uses NetworkX 3.5 maximum-cardinality matching with weights
ignored. Tiny independent recursive enumeration yields two matchings of
weights 64 and 16. The extracted residual has partition 20. Original,
allowed-only, and prefactor-times-residual FKT values agree in every row.
Only the tiny row has independent partition enumeration; the larger rows
reuse the existing FKT helper on different graph forms.

The false deletion of allowed carry edge (1,2) lowers the tiny sum to 16.
Keeping forced endpoints after extracting their prefactor incorrectly
produces 320. Both wrong calculations therefore fail equality as required.
Root inspected the full actual source, reports and both logs, and verified
archived SHA-256 identities. The finite checks are not used as an all-k
classification or a substitute for the proof.

## Preserved review corrections and board provenance

The initial source audit `Sc2eb3390c035421a` was returned for changes in
review `V56ce46a9d3d54c6c`. It understated the bandwidth coefficient,
reversed the endpoint permutation convention, and left the prefix proof
conditional. The revision also distinguishes LRT's scalar algebra
assumptions, the whole-pair transfer, and the signed-weight scope of the
modern generic lower bound. The original artifact remains archived; the
accepted revision is `Saa31447e54da42ed`.

| task | evidence submission | review |
|---|---|---|
| Tb7fbe773e0e84d65, Astra derivation | S6b95d87e1e054aea | V27e3c95b6a944269 |
| T73dcfa9d657440de, Sol primary audit | Saa31447e54da42ed | V38184b9512ca406c |
| T808e358d27c54dd0, Sol independent check | S6256c3db366a44cc | V8f2c2aa3ad4348c7 |
| T9b545ea306594ff5, root implementation | Se0060b9ac0b0464f | V698bae7d1f3c461c |

Astra's derivation is
`out/agent-board/workers/A23cd3ca781c347c5/reductions_and_recurrence.md`;
Sol's revised audit is
`out/agent-board/workers/A6698903e570c43eb/rotated_zigzag_dimer_audit.md`.
Workers made no canonical edits. Root remains the single source writer.

## Validation and checkpoint scope

Core ran first and after the new experiment; both terminal exits were zero
and both logs end in ALL TESTS PASSED:
`out/carry_matching_core_initial.log` and
`out/carry_matching_core_final.log`. The two bounded experiment reports
supply the affected science checks. Production helpers were unchanged and
the other eight science suites were not rerun. No native crash occurred in
these bounded runs; TODO34's previous interpreter failure is unresolved.
Documentation indexing exits zero (`out/carry_matching_reindex.log`), and
all ten documentation checks pass (`out/carry_matching_docs_check.log`).
Scoped diff whitespace, new-file whitespace and source syntax checks pass.
All four board tasks are accepted and closed, with no experiment process
remaining at this checkpoint.
No manuscript, dependency file, default, host setting or firmware changed,
and nothing was committed or published. Preserve the dirty worktree.
