---
code: DM
date: 2026-09-12
title: "Explicit growing carry graphs expose a stronger dimer baseline"
outcome: nonzero exact arithmetic family and grid minor proved; direct dimer reduction removes generic gadgets; independent verifier defects retained and corrected
claims: [C87]
todo: [48]
---
# DM — From an existence argument to an explicit arithmetic family

C87 owns the explicit permutation, surviving grid minor, all-bit dimer
formula, autocorrelation/complement identities and local rank obstruction.
This investigation completes TODO48 at bounded construction and exact-query
scope. The remaining reduction and cost discriminator belongs to TODO49.
The user's general simulation goal remains active: no practical memory
advantage or breakthrough is established.

## Discovery and stronger baseline

Astra found a direct rotated-zigzag construction, avoiding extraction of
the earlier published helix family. Its two nested paths have an explicit
planar embedding and a grid-minor certificate. Root then noticed that the
certificate avoids the four endpoint vertices removed by the all-bit
query, for the nondegenerate growing family. Astra independently confirmed
that stronger statement in a separate preserved supplement.

In parallel, root and Sol independently derived the all-bit specialization.
The three allowed transformed entries select exactly two incident edges.
Complementing them turns the problem into a perfect-matching count on the
architecture itself. Its sign is global and its residual edge weights can
be positive integers. This supplies both a nonzero lower bound and a
stronger baseline than C86's universal local gadgets. It does not test
internal sign cancellation, even though odd-width coefficients are negative.

Sol also audited affine/product alternatives. The equal-mask local factors
are W tensors, obstructing those invertible local basis reductions. Root
read the primary definitions, W product-decomposition argument and graph-
state theorem, and checked a direct matrix-slice rank proof. Local rank
does not exclude scalar simplification. No publication of this exact
arithmetic application was identified in the bounded audit, and no
priority or restricted nonplanar hardness claim follows.

## Main exact experiment

Source: `experiments/experiment_planar_query.py`. Reproduce with:

```bash
uv run --with networkx==3.5 python -m experiments.experiment_planar_query
```

Prediction M441d631a358748fe precedes run R71164a27988746b5. Raw report and
log are `out/planar_query_report.json` and `out/planar_query_main.log`.
The experiment fixes the permutation rule, all four all-bit masks and
uniform incoming carries, varying only k. The direct weighted graph has
the following exact results:

| k | Word width | Direct graph vertices | Exact coefficient |
|---|---|---|---|
| 1 | 8 | 12 | 5/1024 |
| 2 | 32 | 60 | 833389/140737488355328 |
| 3 | 72 | 140 | 2110472116271192473/5316911983139663491615228241121378304 |

The first row agrees with independent integer cyclic-autocorrelation
enumeration and the generic C86 implementation, which uses 72 gadget
vertices on that row. This is a graph-size observation, not an allocated-
byte or speed comparison. Larger scalar rows use the proved direct
reduction and existing exact FKT helper; they are not independent full
basis replays. No large physical circuit or statevector was simulated.

The odd-width identity control gives -1/16 and catches removal of the
global sign. The same first-row architecture with zero masks gives one,
while its all-bit query is nontrivial. All four harness checks pass with
exit zero and no protocol warnings. Integer/Fraction arithmetic is used
throughout. No performance pilot, timing comparison or peak-memory claim
was made; further matching-polynomial simplifications remain possible.

## Independent geometry evidence

Submission Sf39713394d7747fd is under
`out/agent-board/workers/A3469c04194384bc8/`, including
`explicit_family.md`, `geometry_check.py`, `geometry_report.json` and
`geometry_check.log`. Prediction M3ccbf05179e7463d preceded run
R5ec96218161949ae. The proposed larger fixture cap was narrowed before
execution, with its original prediction retained on the board.

All three k fixtures pass graph-structure, planarity and explicit branch-
set checks. Identity matching fails the same frozen certificate, as
required. All three harness checks pass, exit zero. These checks concern
the architecture, not weighted scalar cost. The endpoint-deletion result
is an algebraic supplement, not an extra all-fixture execution: its k=1
exception is explicit. `endpoint_deletion_supplement.md` is archived in
the main submission, preserving the original geometry submission.

## Independent scalar verification and retained defects

Luna's evidence is under `out/agent-board/workers/A0b16f9656409416d/`.
Final accepted submission S6bc0372b42094e72 was reviewed in
Vc84e1f7ff96c4df1 after source, report, logs and all archived hashes were
inspected. The authoritative corrected scalar source is
`experiment_planar_carry_verified.py`, with
`planar_carry_verified_report.json`, `planar_carry_verified.log` and
`findings_carry_placement.md`. Run R52e9e348a02d46fe follows prediction
M6f1ccdbed5a94eed and exits zero with all four checks resolved.

At fixed width three, corrected direct sums include both initial carries
and all 2048 word/carry assignments per query. All 48 equal-mask cases
(eight masks, six permutations) agree with the autocorrelation identity.
Nonuniform highest-bit mismatch fixtures vanish; the matched identity
case with all masks four is nonzero, -5/16. The least-significant-bit
identity query with both carry characters gives +1, while the old
omitted-carry reference gives zero. That last must-fail control directly
exposes the reference defect rather than merely checking a final sign.

The complete correction history remains available:

- R82bfa3fd1a704164 failed before measurement because the base environment
  lacked optional NetworkX. `planar_local.log` preserves that exit-one
  failure; the dependency-pinned run R20f316f154e34b44 then passed.
- Its original highest-bit zeros were nondiscriminating because uniform
  characters independently forced them. Review V015a5b93a3bf4814 requested
  stronger fixtures. `experiment_planar_strengthened.py` and run
  Rb60205d99b1242d7 preserve that first passing correction.
- Root's actual source inspection then found both independent direct
  helpers computed S=A+B outside the incoming-carry loop. Their stated
  full-space contract was false, despite all selected rows passing.
  Review Vcf5473ab24f945ad required a corrected reference and explicit
  carry-sensitive witness. The old files and reports remain unchanged;
  C87 explains the complement symmetry that hid the omission.
- `experiment_planar_carry_corrected.py` included the carry but introduced
  a second verifier defect: it compared every permutation with g(s)^2.
  Run Rcc555d6466a74080 exits one, preserving 26 autocorrelation mismatches
  and the separately successful carry witness. The final verified source
  uses g(s)g(P(s)) and corrects both defects.

The initial local tensor enumeration is unaffected by those word-reference
defects. It independently checks all eight K/R mask patterns and their
Hadamard entries. Its scalar rows are historical until confirmed by the
final corrected reference. Neither worker source is used by the main
experiment, whose independent autocorrelation and generic FKT comparison
were reviewed separately. Passing tests alone would have missed both the
nondiscriminating fixtures and the omitted carry.

## Review and validation scope

Root accepted the geometry submission in Vcf91b8941e7340a7 after inspecting
proof, source, logs and frozen hashes. Sol's source submission
S002aade93c324049, including `scalar_reduction_audit.md` under
`out/agent-board/workers/A697ac92a69164a93/`, was accepted in
Vdca4ab059eec4cc6. The rank and dimer deductions were independently checked;
the audit's general hafnian and Thue--Morse references are background leads,
not implemented algorithms or closed-form exclusions for this family.

Main submission S35efe49a35c84678 archives the experiment, report, logs
and endpoint supplement. Independent Astra review V312376d8993d4c92
inspected the contents and hashes and accepted the restricted result.
Core passes before execution and after implementation, with logs
`out/planar_query_core_initial.log` and `out/planar_query_core_final.log`.
Other science suites were not rerun; TODO34 retains the prior native crash.
The failed worker executions above are not relabeled passing suites.

Documentation regeneration and validation are recorded in
`out/planar_query_docs.log`. All four bounded board tasks are accepted and
closed. No manuscript was changed and nothing was committed or published.
The canonical worktree remains uncommitted.
