---
code: MG
date: 2026-09-12
title: "Interacting carries reach planar matchgate contraction"
outcome: exact restricted Pfaffian implementation validated; explicit growing nonvacuous arithmetic family and practical advantage remain unresolved
claims: [C86]
todo: [47]
---
# MG — A planar algebraic route beyond one carry chain

C86 owns the full-space coefficient contract, complement-parity proof,
explicit matchgate construction, resource accounting, geometry deduction
and arithmetic baselines. This investigation completes TODO47 at bounded
audit and implementation scope. TODO48 owns the next discriminator. The
general research goal remains active; no breakthrough or practical memory
advantage is established.

## How the direction changed

The initial audit followed the actual `toffoli_arith` repeated-addition
macros. Read-only constant additions group even with arbitrary scratch,
provided their actual dirty translation is retained. The modular flag
schedule obstructs replacing the full-space circuit by its clean-code
identity. The independent fixture below checks that distinction.

The coordinator then fixed a different explicit family: two independent-
operand additions with a bit permutation between them. This is arithmetic
connected by SWAP/CNOT gates, not the dirty-scratch controlled multiplier.
Its two carry paths expose ternary factors. Complement symmetry suggested
a local Hadamard basis change; the resulting pure-parity factors led to
the known matchgate/FKT algorithm. The source and mathematical audits
reviewed this mechanism before the main scientific run.

The geometry audit ruled out assuming every planar matching of the two
paths has bounded treewidth. Its source is the published Eppstein article,
whose final construction is absent from the older arXiv HTML. The root
independently read the relevant definition and construction and checked
the bounded-edge cleanup argument. No growing adjacency lists or actual
mask family were extracted. C86 states exactly what the deduction proves.

## Main bounded exact experiment

Source: `lab/planar_carry.py` and
`experiments/experiment_planar_carry.py`. Reproduce with:

```bash
uv run --with networkx==3.5 python -m experiments.experiment_planar_carry
```

The registered prediction is M291997fc0e8d494d and run is
Ree5cbd1efac84eb5. Raw evidence is `out/planar_carry_report.json` and
`out/planar_carry_main.log`. The run exits zero with all seven harness
checks resolved and no protocol warnings.

- Thirteen local carry/degenerate parity fixtures give 104 external-deletion
  signatures, each checked against independent tiny weighted-matching
  enumeration. This includes zero tensors and the forced-edge case.
- Four signed K4 edge-weight fixtures agree with FKT. An unoriented
  Pfaffian gives one where the perfect-matching sum is three. Replacing
  a negative edge by its magnitude changes exact cancellation from zero
  to four. Nonplanar K3,3 is rejected. These are the three must-fail controls.
- At fixed width three, every one of six intermediate bit permutations
  is compiled from actual Cuccaro logical gates and SWAPs. Existing
  `walsh.classical_permutation` replay matches independent integer word
  arithmetic, including both output-carry XORs, on all 8192 basis inputs
  per 13-qubit circuit.
- The same twenty predetermined mask fixtures are used for every
  permutation. All 120 signed coefficients agree exactly; 41 are nonzero
  and 18 are negative. These are selected masks, not exhaustive coverage
  of all mask tuples. The fixed arithmetic layout and masks leave the
  permutation as the only varying experimental parameter.

All coefficient and matching calculations use integers or `Fraction`.
There is no tensor SVD, floating-point cutoff, large unitary or performance
pilot. Tiny recursive matching enumeration is confined to local gadgets
and four-vertex controls. The six small architecture graphs are planar;
this experiment does not observe unbounded width or a memory separation.

## Independent actual-macro fixtures

Luna's accepted submission S5a46bedfeb2846e6 is under
`out/agent-board/workers/A019b54646717444b/`, including
`experiment_interaction_tests.py`, `interaction_tests_report.json`,
`interaction_tests.log` and `findings.md`. Prediction M962ffe53b92d4a57
preceded run Rc0d5402c8edf4a38, which exits zero. All four checks resolve
without protocol warnings.

The fixed ToffoliModExp layout uses N=3, a=2 and one exponent bit: eleven
qubits and 2048 basis labels. Twelve constant-add variants cover constants
zero through two and four unchanged control choices. Three modular-add
variants use those constants with fixed two-bit controls. Every full output
label agrees with the independent formulas, including restored wires.
The same-control pair of constant additions commutes on every label.
The deliberately incorrect clean-code shortcut for an uncontrolled add-one
disagrees on 1792 labels; t=2, c0=0 provides an explicit increment-three
witness. No dense full-space macro spectrum was built.

The math audit also gives signed coefficients for a single CNOT interposed
between repeated additions and a dirty modular noncommutation witness.
Those are algebraic deductions in the audit, not separately executed
experiments. The single-bridge coefficients reduce to one majority carry,
so they do not supply a new memory opportunity.

## Source scope and board review

Astra's arithmetic audit is submission S96436d34be3a4187 under
`out/agent-board/workers/A6e16adbdedcf43c2/audit.md`, accepted by review
Vb2e9d8ac09a4494f. Its later `planar_geometry_followup.md` is archived in
the main submission rather than retroactively inserted into that archive.

Sol's primary-source audit is Sa25e89af4a7f42a7 under
`out/agent-board/workers/A78a9281646b640fe/audit.md`, accepted by
V4088979f401b4412. A correction to TODO47's opening premise matters:
Wallén Section 3.4/Theorem 3.5 concerns parallel output concatenation and
its tensor-product rational representation, not a bounded representation
for arbitrary serial ARX/CNOT composition. Serial correlations require
summing signed intermediate masks. Trail selection and squared or
averaged linear potentials do not replace that exact signed hull.
Weighted-automaton minimization also charges construction of its explicit
input representation; small minimal dimension is not a free oracle.

The matchgate audit reads the parity/identity and planar realization
results of Cai–Gorenstein and the holographic contraction framework of
Cai–Lu. C86 cites the published primary body and gives the explicit gadgets
actually implemented. Genus-dependent Pfaffian sums appear as an
unimplemented literature lead in the source audit; no such extension was
built or benchmarked. The exact arithmetic specialization is our deduction
from known machinery, with no priority conclusion.

Root inspected the actual independent source, report and log before
accepting it in V4f621f7e07574ac1. Main submission S62f1bf0bfc8a48cb
archives the implementation, experiment, report, logs and geometry
supplement. Astra's independent review Ve6e62715c0fd4dd4 inspected the
archived contents and verified version hashes, sign calibration,
normalization, degenerate gadgets, controls and scope. All four bounded
board tasks are accepted and closed. Acceptance validates the recorded
restricted result; it does not establish novelty or a broad simulator claim.

## Validation and disposition

Core passes before scientific edits and after implementation, with logs
`out/interacting_carry_core_initial.log` and
`out/interacting_carry_core_final.log`. Both bounded experiments pass;
no failed scientific execution required correction in this investigation.
Other science suites were not rerun. TODO34 retains the prior native-crash
result; these small runs make no host reliability claim.

Documentation regeneration and validation are recorded in
`out/planar_carry_docs.log`. No manuscript was changed, and nothing was
committed or published. The next question belongs only in TODO48; this
checkpoint must not be relabeled a demonstrated CNOT memory breakthrough.
