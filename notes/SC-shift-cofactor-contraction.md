---
code: SC
date: 2026-09-12
title: "A final signed shift exposes a polynomial-storage contraction of enabled dirty prefixes"
outcome: exact growing-depth construction implemented for the selected target and checked against gates and direct cell arithmetic; no measured memory advantage
claims: [C93]
todo: [50]
---
# SC — From nested prefix shifts to one digit carry

C93 owns the mathematical results, charged bounds, normalization and scope.
TODO50 owns the remaining research question. This checkpoint is progress
toward the user's broader simulation goal; it does not complete that goal.

## Derivation and independent scrutiny

The authorized swarm retained Astra for algebra, Sol for the endpoint/source
audit and a second Sol worker for the quotient and independent cell checks.
Root remained the sole canonical writer. All work started from C92's frozen
dirty Mersenne equations and the TODO50 target. No further mask search or
conditional-piece-count benchmark was run.

Astra first conditioned on ternary signed enable digits. Each sector makes
every order-specific prefix shift constant while pinning only the required
scratch digits, allowing the existing affine arrangement proof to apply.
This yields a growing enabled algorithm even for arbitrary permutations.
Its unequal sector cardinalities are essential. A tempting rank-three
interpretation was rejected algebraically: the local matrix has rank two,
which neither establishes a lower bound nor preserves the useful fixed
prefix catalog in a rank decomposition.

Root then proposed conditioning on the final signed difference instead of
the full signed-digit word. One binary addition carry reconstructs earlier
differences during the low-to-high scan. Root derived the compatible backward
history for complete reversal. Astra independently audited both directions
and simplified the target phase so that the final count need not be guessed
in an outer loop. Sol independently checked the normalization, index
conventions, membership tests and carry terminal condition before execution.

In parallel, Sol derived the smaller endpoint carry alphabet and compared
it with the stronger interval and sector constructions. The other Sol
worker derived the exact P-odd quotient, scratch-top path projection and
one-macro cancellation. These results prevented an unsupported even-depth
selection law and supplied a charged finite-array competitor. C93 owns
their proof statements; neither worker's side construction was implemented.

Root inspected the primary arithmetic BDD construction and product bound
in Bartzis–Bultan, and Kiefer's supplied-automaton minimization theorem,
including their proof bodies. Sol inspected the archived Wallén body cited
in its memo. C93 links the applicable sources. This audit supports the
stated known-method comparisons, not a claim to published novelty or a
literature-wide optimal algorithm.

## Bounded candidate implementation

Source: `experiments/experiment_shift_cofactor.py`. Reproduce with:

```bash
timeout 60s uv run python -u -X faulthandler -m experiments.experiment_shift_cofactor
```

The core gate ran first and exited zero; its complete log is
`out/shift_cofactor_core_initial.log`. Under task `T175f5d83273c42d3`, attempt
`Abe0d09e5bb9e497a`, prediction `M9f62a9a5b8af4d08` preceded run
`R64537598a64a4b6a`. The actual terminal exit was zero. All seven checks
passed, with no harness warnings. Root read the entire raw log
`out/shift_cofactor_main.log` and every report row in
`out/shift_cofactor_report.json`; the independent worker also read both.

The family was fixed at N=7,a=1,n_exp=1, varying only prefix depth q=1,2,3
and checking both declared schedules. The target is C93's full dirty-space
coefficient. Every signed shift sector was compared with the actual logical
gate map, including its unsigned mass and independently derived cardinality.
All 50 sector pairs and six full-space scalars agree exactly. Each gate map
enumerated 16,384 physical labels; unused multiplicand spectators were
removed only after checking divisibility. The candidate uses integers and
exact `Fraction` geometry; it does not use a floating Walsh transform.

The one-macro zeros agree with C93's proof. The three-macro values and the
omitted-incoming-carry control reproduce the frozen DS common-mask row;
DS remains their numerical home. The newly observed two-macro coefficients
are 3/128 ascending and 7/256 descending. These values were not predicted
before measurement; only equality with the independent gates was predicted.
They are finite fixture evidence, not a nonzero growing-depth formula.

The run constructed 22,996 cells: 480 open-strip cells and 22,516 singleton
slice cells. Its peak line count was 138, peak event count 145, peak current
DP count 48, and maximum integer scratch span two. These are descriptive
internal counts, not total allocated bytes or a measured performance win.
Current/next dictionaries, geometric events, arithmetic objects and the
validation reference also occupy storage. The reference and report retain
dense/sector data only for this bounded correctness test; C93's algorithm
streams sectors and cells. No size or allocation benchmark was attempted.

Both compulsory wrong references were meaningful: complete reversal changed
the selected nonzero three-macro scalar, and omitting arithmetic h changed
it to the DS wrong-reference value while retaining the physical h phase.
Neither mutation was used to replace the actual circuit under study.

## Independent cell verification and preserved failure

Task `T4f8b5243704249b0`, attempt `A0d41a707f3324d11`, used direct integer
membership and C89 `macro_image` instead of the digit DP. Its source is
`out/agent-board/workers/A0d41a707f3324d11/verify_shift_cofactor_cells.py`.
Reproduce with:

```bash
timeout 60s env PYTHONPATH=. uv run python -u -X faulthandler out/agent-board/workers/A0d41a707f3324d11/verify_shift_cofactor_cells.py
```

The first prediction `M4d592171975541a9` preceded run `Rd086f4b1e8d54e5c`,
which exited one before any candidate/direct comparison. Its deterministic
selection policy incorrectly required each depth/direction/sign slot to
contain both open and point geometry. One slot was empty. The complete
source and log are preserved under the attempt's `revision1/` directory,
and the board also archives that failed log. This was a verifier selection
failure, not a failure or pass of the scientific prediction.

Revision two changed selection to global deterministic coverage quotas.
Root inspected the source diff: direct arithmetic, membership rules and
wrong-reference definitions were unchanged. Prediction `Mdd0522e92d8d43cf`
preceded run `R811f671c2a554fa1`, which exited zero with all six checks
passing and no warnings. The same attempt contains the complete final log,
`verify_shift_cofactor_cells_report.json`, and a source/report review memo.

All 40 selected cells matched both signed numerator and unsigned mass over
74,240 direct loop visits. Coverage included all p,h pairs, both schedules,
all shift signs, the three tiny depths, 35 critical singleton cells and
five open cells containing multiple integer scratch values. One separately
fixed N=63,q=1 cell spans 30 scratch values; this stresses geometry rather
than constituting a width sweep or growing-depth benchmark.

Making the upper boundary inclusive changed nonzero row 8 from (-2,2) to
(-2,4), where each pair is (signed sum,mass). Omitting the true signed-range
condition for negative S changed row 27 from (-2,2) to (0,4). These controls
test boundary multiplicity and signed-shift handling even when a final
full-space sum might hide a cancellation. The direct reference shares the
reviewed C89 arithmetic and candidate geometry generation; it independently
tests cell summation, while global partition correctness has separate proof
and the complete tiny gate-sector checks.

## Frozen evidence and review

Root read the complete proof/source memos, raw logs, reports and preserved
verifier revision, then verified archived evidence hashes. The historical
failed-run log has the original current-path name in its run manifest;
its bytes match the preserved revision1 file, not the later successful log.

| submission | scope | final disposition |
|---|---|---|
| S968a64ed29064891 | Astra ternary and final-shift proofs | accepted; task closed |
| Sc60c0d730508495d | Sol endpoint frontier and source audit | accepted; task closed |
| S0fdf91dfea274d08 | Sol quotient proof and charged recurrence | accepted; task closed |
| Sda9284db313b47f8 | independent cell source, both runs and full evidence | accepted; task closed |

The source audit has preserved evidence-version corrections. Its original
submission `Sb0ab226595cf4b5c` was replaced after a post-freeze baseline
paragraph addition, and `Sbbe3d9aa7acb42b3` was replaced to include the
subsequently requested implementation review. Neither replacement corrected
a mathematical failure. The original body is preserved under revision1;
the accepted final submission binds all three artifacts immutably.

Key SHA-256 identities, with the board manifest owning all artifact paths:

| artifact | SHA-256 |
|---|---|
| candidate source | 7dcc149a9a6ebc8d26c042ad40cf8cecd185f9f25d598a815eacfbffc4500016 |
| candidate report | d74c6e64fa468a60db72ccefe185b8e89599b5c2320f1f1b6ec24ec12108ecda |
| candidate raw log | 43b2cf01dd8a31ada70729752881e448b49bb2425c7b717940d166d6b8f6f927 |
| final direct verifier source | 2b6735c9bcd17c1b1c8f6be6f9ea0e6091e43eca42631e3485f33419f9228114 |
| direct verifier report | c79295103e71278e50b2185fbf0399b544bff19423da4ac5b6b2a547c2e39cbf |

Final syntax and whitespace checks pass. Regenerating indexes with
`uv run python tools/reindex.py` and running `uv run python tools/check.py`
exited zero; all ten documentation checks pass. Logs are
`out/shift_cofactor_reindex.log` and `out/shift_cofactor_docs_check.log`.
The successful science checks were not repeated after documentation-only
edits. Astra also audited the final C93 mathematical transfer in board
message `M13430e0f72904506`, with no correction requested.

Coordinator write task `T175f5d83273c42d3` is the authoritative record of
independent final review and integration. Its submission/review/closure
manifest owns the exact integrated version and final validation logs.
Worker tasks were closed before updating their snapshotted TODO50 input.
The scope is a bounded experiment and documentation; no production helper,
manuscript, dependency, host setting, commit or publication was changed.
Other eight science suites were not rerun; TODO34's earlier native crash
remains unresolved. Preserve the existing dirty worktree.
