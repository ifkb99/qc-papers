---
code: QC
date: 2026-09-12
title: "Exact quadratic cells compress selected Walsh queries, but the recognized formula is cheaper"
outcome: tested restricted evaluator and closure proof; known theory and stronger-baseline loss retained
claims: [C83]
todo: [44]
---
# QC — Compressing the observable after the CNOT frame

This goal turn follows C82/CF's allocation result. It made concrete progress:
the coordinator implemented a restricted exact query evaluator, Astra audited
the closure proof and frozen implementation, Sol checked primary prior art
and stronger baselines, and Luna supplied independent bounded tests. The
broader breakthrough goal remains open. C83 owns the theorem, API, asymptotic
scope and prior-art boundary; this note owns executed evidence and discovery.

## Matched scalar comparison

`experiments/experiment_quadratic_cell_memory.py` and
`out/quadratic_cell_memory_main.json` / `.log` own the measurements. The small
sweep freezes m=4, n=14 and a deterministic 19-CNOT prefix, varying only p from
0 through 4. Registers are (u_i,v_i,a_i), shared control w and observed target
t. The nonlinear schedule is Toffoli(w,a_i;u_i) for i<p, then all
Toffoli(u_i,v_i;t). Before the CNOT prefix, the pullback phase is

    f=t+sum_i u_i v_i+w sum_(i<p) a_i v_i.

The requested result is ONE full-space signed physical Walsh coefficient,
including its normalization. Every method receives the actual Circuit and
physical masks. The cells additionally receive output cut {w}; cut discovery
is not measured. The stronger formula first inspects the circuit's layout,
observable and gate sequence and converts the query through the prefix.
Thus it is not given an expanded polynomial for free. It sums over w with
sign, and analytically eliminates u_i,a_i using the constraints
v_i=z_ui and z_ai=w z_ui for i<p (otherwise z_ai=0).

Four selected signed queries agree among cells, sparse PPS and the recognized
formula at every p. Construction plus one negative query is measured once
per method with `tracemalloc`, following garbage collection:

| p | sparse PPS peak bytes | streamed cells peak bytes | recognized formula peak bytes |
|---:|---:|---:|---:|
| 0 | 63264 | 4472 | 1000 |
| 1 | 139072 | 4320 | 1000 |
| 2 | 225088 | 4256 | 1000 |
| 3 | 235896 | 4256 | 1000 |
| 4 | 277296 | 4256 | 1000 |

Two branches suffice throughout, processed one at a time. The p=4 cell peak
is about 65 times smaller than sparse PPS, but the recognized formula is
smaller still. Sparse PPS constructs its full dictionary to answer the query;
this is a measured route to the same scalar output, not a lower bound on
memory for that output. Peaks include result/counter objects and conversion
to `Fraction`, but exclude the prebuilt Circuit, interpreter, native RSS and
prior allocations. They are single-run Python-allocation observations, not
stable benchmark constants or timing results.

One separate restricted endpoint uses m=p=32, n=98 with the same prefix-depth
rule. Cells and the charged formula both return exactly -1/2^33, with measured
peaks 25692 and 2100 bytes respectively. Cells perform two full replays:
38 CNOTs, 64 affine Toffoli transports, 64 quadratic Toffoli updates and 64
Gauss pair eliminations. Formula recognition records 84 gate inspections.
No sparse or dense run was attempted there; this endpoint changes width as
well as the gate count and is not an extension of the frozen-width sweep.

The exact support count for this family was independently derived in the
code review. Each fixed-w Fourier support has 2^(2m+1) masks; the intersection
has 2^(2m-p+1), and half the intersection cancels in the coherent sum. Hence

    S(m,p)=2^(2m+2)-3*2^(2m-p).

It matches the small enumerated supports 256,640,832,928,976. The large
report's 2^66-3*2^32 is DERIVED, not enumerated, allocated or a scalar-query
memory requirement. Known stabilizer structure explains this compression.

Both controls are nonvacuous: a tiny branch sum +1/4-1/4 is zero while the
sum of magnitudes is nonzero, and the uncut cubic intermediate sign raises
the named certificate escape. The memory experiment finishes 28/28 checks,
exit zero, without harness warnings. Tiny sparse references use exactly
representable dyadic floats; no large-coefficient threshold claim follows.

## Independent checks and failed drafts

`experiments/experiment_quadratic_cells.py` integrates the independent worker
draft with stronger coordinator assertions. Its reference is existing exact
classical permutation replay plus an integer Walsh butterfly, and direct
exhaustive upper-triangle Boolean evaluation for the Gauss/product helpers.
`out/quadratic_cells_main.json` / `.log` record:

- 961 seeded quadratic-form/extra-linear combinations through dimension five,
  and 200 affine products, including repeated-variable reductions;
- 42 circuit/cut/query cases through seven qubits, including eight seeded
  mixed circuits: 41 exact matches and one supported rejection, with six
  negative and 35 zero successful coefficients;
- full coordinate cuts required to succeed, plus a direct one-dimensional
  nonparallel-cell fixture comparing actual tangent sets and signed images;
- explicit cubic and curved-cell rejections, and a zero coefficient whose
  incoherent raw branch magnitude sum is 16.

The integrated audit passes all four harness checks. This is a finite
correctness corpus, not a scaling sweep or exhaustive coverage of arbitrary
affine cells. Core section 7 separately checks full coefficients against
independent dense conjugation, signs, three partition choices, two named
escapes, a 65-bit physical mask, and pre-allocation branch-cap rejection.

The worker needed two import-path fixes and an unexpected P2 coverage fix
because its first selected coefficients were all zero. A draft comparison
of unequal matrix rows also failed to establish different tangent spaces;
the corrected control cuts both b and t in Toffoli(a,b;t). Worker failed logs
remain sealed in run records R33bb39b1bf114121, R75e2681d7182407e and
R31740a0ee3124785. Their original draft scripts were overwritten, a provenance
limitation. The immutable third run's phrase 'intentionally failed' was
corrected in the replacement submission: that coverage failure was unexpected.

Coordinator strengthening initially referenced a helper name absent from the
submitted worker version, raising `NameError` before P3. Both its script and
log survive as `out/quadratic_cells_integration_failed.py` / `.log`. The final
version shares an independent upper-triangle evaluator across its reference
checks and completes successfully. These are verifier defects, not evidence
that unsupported cells can be approximated safely.

## Review, validation and limits

Board math task T2925f4e09c34412b accepted Sfe307f2cefa843f9; baseline task
T156f581d0eb24503 accepted S344a6f71256f47d4. Test task T0815738651e84391's
first submission had an incorrect core-log pointer; replacement
Sb04c445898594d1c is accepted with corrected archived evidence. Main task
T8b1c470bc31a498a submitted Sab54242e0eea44e3, independently reviewed by Astra
as Va2e81ab4b4ed4118 against frozen source/report hashes. That review audited
the solver, gamma, Gauss sign, normalization, formula and support derivation;
it did not execute another performance comparison.

The initial core gate passed before science changes
(`out/quadratic_observable_core_initial.log`). The final core including new
regressions passes (`out/quadratic_cells_core_final.log`), as do the two bounded
experiments above. No full legacy or other science suites were rerun. CF's
native-crash/incomplete legacy result remains with TODO34. Documentation is
validated with `tools/reindex.py` and `tools/check.py`; its final log is
`out/quadratic_cells_docs.log`. No commits, publication or manuscript edits
were made in this investigation.

The cheaper formula and known Gauss/stabilizer decomposition exclude a broad
breakthrough claim from these results. Generic phase-sensitive CCX splitting
and merging has not been benchmarked against the local certificate. TODO45
owns the next discriminator; repeating this recognized family's width pilot
would not resolve it. TODO42 remains a separate sampling-cost question.
