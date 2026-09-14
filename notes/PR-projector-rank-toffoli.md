---
code: PR
date: 2026-09-12
title: "Generic stabilizer measurement rank closes the proposed local Toffoli recognition advantage"
outcome: proved equivalent generic criterion and eight signed replacements; no breakthrough or end-to-end benchmark
claims: [C84]
todo: [45]
---
# PR — Stronger recognition baseline after C83

The previous goal turn made concrete progress at C83/QC; its cheaper formula
left a possible distinction in circuit recognition. This turn resolves that
distinction algebraically and numerically. C84 owns the projector-reflection
proof, three-bit tableau criterion, exact replacements and equivalence to
C83. The broader user goal remains active. Known stabilizer ingredients and
the same acceptance set prevent a breakthrough claim from these results.

## Authoritative experiments

`experiments/experiment_toffoli_stabilizer_rank.py` and
`out/toffoli_stabilizer_rank_main.json` / `.log` compare C83 against a stronger
generic recognizer on the SAME prepared input cell. Thirty-four five-qubit
cases include explicit witnesses and 24 seeded affine-transformed quadratic
cells. A Stim preparation circuit independently checks the cell's input ray;
Stim Pauli expectations give exact joint probability. Direct basis-index
Toffoli action then checks both the accepted C83 update and each proposed
physical Clifford replacement on the same signed vector.

All 34 acceptance decisions agree: 17 accept and 17 reject. Every one of the
eight exact replacement labels occurs. Every rank scan reads 15 tableau
entries and stores only three small basis columns; these are operation counts,
not measured runtime or allocation peaks. Constructing the tableau AND a
separate Stim simulator, seven expectation calls for compatibility, any
conversion and signed extraction are extra. There is no claim that this
experimental arrangement is an optimized complete simulator.

The decisive pair freezes n=5, full support and a linear target sign; adding
only the target-spectator quadratic edge changes rank two to rank three and
causes the C83 cubic-sign escape. Its input/output overlap becomes 3/4,
incompatible with a distinct stabilizer ray. The separate fixed-target/free-
controls case triggers curved-support escape. The p=1 case has overlap -1;
its exact -I replacement succeeds while dropping that sign fails vector
equality. The experiment completes 39/39 checks with a final JSON report.

These vector checks use floating normalization, not rational arithmetic.
Stim's float32-origin input-ray comparison uses tolerance 2e-7; signed
replacement/cell comparisons use 1e-12. Exactness is supplied by C84's proof
and exact Stim Pauli expectations. Zero observed vector errors do not make
the complete experiment an exact-arithmetic execution. Incidental harness
elapsed times are not a comparative performance benchmark.

`experiments/experiment_projector_dense.py` and
`out/projector_dense_main.json` / `.log` integrate Luna's independent audit.
Thirty-three deterministic Clifford circuits through n=5 prepare pure inputs;
direct basis-index Toffoli produces outputs. The detector enumerates all
Hermitian Pauli expectations, one matrix at a time, and counts magnitude-one
values. A pure state is stabilizer iff this count is 2^n; the detector does
not use the projector criterion. Projector probability is computed separately
both by matrices and by the control/target-minus projection.

Observed p classes 0,1/8,1/4,1/2,1 have counts 5,9,11,5,3. Every row agrees
with the criterion. Ten inputs have relative complex phases that cannot be
removed by a global scalar. The corrected detector's closest nonunit
expectation has distance about 0.5 from magnitude one, far above its 2e-9
tolerance. Maximum matrix dimension is 32. The projector calculation holds
three Pauli matrices plus identity and bounded temporaries; the one-matrix
bound applies only to the detector. All five integrated harness checks pass.
The global-input-minus check only tests linearity; the explicit p=1 output
minus and main experiment's signed replacement test carry the phase evidence.

## Failures and provenance

Main run R9aae56f94bfc4f38 printed successful science checks but exited one
when JSON serialization encountered a NumPy boolean. Its complete failed
source/log survive as `out/toffoli_stabilizer_rank_failed_v1.py` / `.log`.
The only correction wraps the recorded input comparison in native `bool`.
R8185697158cd4875 then exits zero with its complete report. The initial run
remains a failed run, not an accepted report.

The independent worker's import failure survives in its attempt directory
and sealed run Ra684a770b5344bbc. A pre-sign-check draft was saved later;
the exact source of that original failing execution was not sealed then.
Rbbd1ea0c734d46dc and R1b6385dfd5184dd2 are both terminal zero exits, with
the second adding an explicit p=1 sign check. The former finish was recorded
late; coordinator checked the records before acceptance. Their log text is
identical because the added assertion does not alter the printed summary.

Coordinator inspection found two misleading unused worker metrics: the
'nearest nonunit' calculation included the identity and was always zero,
and the one-live-Pauli metadata did not describe the projector workspace.
Integration corrected both, removed an overwritten unused roll calculation,
and added actual relative-complex coverage. Neither faulty metadata field
was used to establish the worker's p criterion. Failed evidence and bounded
numerical limits are retained rather than recast as exact validation.

## Swarm and implementation audit

Astra's proof task T2a446b94621a40a4 accepted S0431af745c4b40d2. Sol's primary
literature/source task T302b8787d6e84fc4 accepted S99c87e6207264a63. Luna's
independent test task T4b16b6fa36e8419a accepted S016cb24e503e4ec7 with the
metric limitations above. Main task T943f7c5ffdf5442f submitted
S3ce37fd134f54c49, reviewed against archived code/report hashes by Astra as
V8ee9d55d5f4c4a20. Only the coordinator edited canonical files.

The source audit is preserved in the baseline worker's `audit.md` and pinned
`src-audit/` checkouts. Stim supplies the direct tableau cell access needed
by the strong recognizer, but its tableau lacks a defined ket global phase.
Cirq's CH implementation explicitly carries `omega` through signed amplitude
and phase updates. The audited QuantumClifford CCZ path uses decomposition
and sparsification, while its tableau inner product returns magnitude only;
the audited Tsim pipeline targets probability sampling and drops scalars
that cancel from that output. These are source/API distinctions, not measured
performance rankings. Exact commits and source locations live in the audit;
C84 cites the supporting primary interfaces and mathematical predecessors.

Neither splitting every Toffoli without merging nor seven independent
expectation solves is the strongest recognition-only baseline. The direct
tableau rank scan already decides exactly the same one-cell boundary. A
phase-sensitive full-output comparison could still measure implementation
tradeoffs, but would not restore that exclusive mathematical advantage.

The science core passed before this work
(`out/toffoli_stabilizer_core_initial.log`). Both bounded experiments pass;
no production library, existing science suite, manuscript or abstract was
modified this turn, and no other science suites were run. The earlier native
legacy-suite failure remains with TODO34. Generated indexes and documentation
are checked by `tools/reindex.py` / `tools/check.py`, with final log
`out/projector_rank_docs.log`. Changes remain uncommitted.

TODO45 is complete at proof/bounded-discriminator scope. TODO46 owns the
next question beyond individual-cell closure; TODO42 remains separate. Do
not repeat this local rank comparison or call a different tableau layout a
new recognition principle.
