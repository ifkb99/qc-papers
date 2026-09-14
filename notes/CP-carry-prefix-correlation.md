---
code: CP
date: 2026-09-12
title: "Actual Cuccaro prefixes reach the known carry-correlation baseline"
outcome: exact restricted implementation and independent checks; isolated-adder novelty closed; no-cut escape prediction refuted
claims: [C85]
todo: [46]
---
# CP — Carry structure beyond a local cell

C85 owns the mathematical contract, prefix invariant, tensor ranks, exact
carry-in reduction, conditioned escape and resource limits. This turn turns
the strongest known baseline into a reusable exact query implementation,
and closes the isolated-adder candidate at bounded scope. It does not
establish the user's requested simulation breakthrough. TODO47 owns the
next discriminator; TODO42 remains separate.

## Main exact implementation evidence

`experiments/experiment_carry_prefix.py` follows the existing harness and
uses `lab/carry_prefix.py`. At fixed width three it visits all 20 logical
prefix lengths, with the same 12 output masks and all 256 query masks.
The output masks are zero, each single physical bit and three fixed
multi-bit masks; this is **not** exhaustive over output masks.

The Wallén automaton and direct signed two-state integer transfer both
equal the existing gate-replay/FWHT reference on all 61,440 coefficient
triples. No mismatch occurs, the maximum symbolic wire support is three,
and 380 reference entries are negative. All six harness checks pass.
The three must-fail controls reject a clean-carry substitution, a
magnitude-only coefficient, and a modified first logical gate.

The reference WHT uses float64, but all butterfly values here are integers
of magnitude at most 256, represented exactly. Conversion to `Fraction`
therefore gives exact equality, with no fitting tolerance. The scalable
query implementations use integer/Fraction arithmetic. One tiny reference
vector is processed at a time. No large-width run, comparative timing,
allocation peak or speedup against an optimized carry automaton was measured.
Incidental harness elapsed time is not a performance benchmark.

Authoritative artifacts are `out/carry_prefix_report.json` and
`out/carry_prefix_main.log`. Board prediction M832062abae924dae preceded
run Rbcb099dfa7b04c8f. Submission S11b9cd6daff340cc archives the actual
module, experiment, report and logs; independent Astra review
Vc26fd393af1a4166 inspected the archived contents and verified their hashes.

## Independent fixtures and a retained failed prediction

Luna's final accepted evidence is submission Sb702c78972f14f94 under
`out/agent-board/workers/Af78ad0485338463f/`. Its authoritative current
outcome is `carry_tests_report.json` and `carry_tests_refutation.log`,
not the earlier log mentioned in the opening of `findings.md`.

Independent logical replay checks every basis input for widths one through
three, including arbitrary c0 and z, against direct integer addition.
All 16, 64 and 256 inputs match respectively. Forty-two exact integer
Walsh fixtures agree with the existing transform, including 32 zeros and
three negative carry coefficients. The clean-carry shortcut disagrees at
each width. These are small correctness fixtures, not an asymptotic sweep.

The original prediction asserted that `quadratic_cell_walsh(qc.inverse(),
B0, B0+A0+c0)` would escape with no output cut. It does not: all three
widths finish without escape and have exact coefficient +1. With the
explicit c0 partition, the walk instead rejects a non-affine cell image at
reverse step three. That is a separately conditioned witness; C85 explains
the normalized states and correct local rank. Endpoint simplicity alone
would not prove that an intermediate state remains a stabilizer.

The final independent run intentionally exits **1**: P1/P2/P4/C1 pass,
while original P3 is explicitly false. Fresh prediction M10bd376bb8cb4bbd
and run R7c827ffb74a246a0 preceded that execution. The harness resolves all
predictions and emits no protocol warnings. An accepted refutation is not
an all-passing science suite.

Failed provenance was retained before correction:

- `carry_tests_initial.log` records the launch import failure. The original
  run Rfc76756fbf604e68 grouped several attempts instead of recording each
  separately; its metadata is not a one-execution provenance record.
- `carry_tests_failed_v1.py` and `carry_tests_final.log` retain the incorrect
  reference z toggle and the unsupported escape assertion.
- `carry_tests_failed_v2.py` / `.log` retain a verifier that omitted
  `exp.finish`, printed success unconditionally and left P3 unresolved.
  Coordinator review requested changes; the final nonzero run corrects
  that reporting defect without relabeling the prediction as true.
- Astra's original `audit.md` in
  `out/agent-board/workers/A65c8f07394814fec/` incorrectly assigned rank
  three to the unrestricted X-eigenstate input. Its unchanged contents
  remain archived in Sca96812e0d6a4221. `witness-correction.md` and accepted
  Sc9d7d6a17df5440d explicitly retract that witness and retain the separate
  prefix proof. The deterministic negative target-X eigenvalue was
  incorrectly treated as a random measurement condition.

The corrected proof is algebraic. The worker's local escape check uses the
existing C83 implementation and is not a second independent all-Pauli
stabilizer detector. The logical permutation checks retain exact relative
signs; they do not measure the compiled circuit's common ket global phase.

## Strongest baselines and audit limits

Sol's source audit Sc9b749597a984524 is under
`out/agent-board/workers/Aa5937e7dcc9b4ee7/`, with actual primary-body text
for Wallén, Quipu/stabilizer frames and Markov–Shi. Wallén's Chapter 3 is
the direct signed linear-correlation source; a differential-addition result
would be a weaker, mismatched citation. The decisive point is that the
arbitrary carry-mask method already covers the intermediate-prefix query,
not merely the finished sum. C85 cites the stable primary sources.

Quipu's phase-aware cofactoring/coalescing already demonstrates arithmetic
recombination, though its clean-ancilla input and reported output differ.
The tensor-network algorithm permits contraction by bit position instead
of chronological gate order. A large chronological branch count therefore
cannot establish a memory lower bound for the present scalar task.

The math audit also explains why a fixed two-Toffoli window is too small
as an asymptotic discriminator: its six Pauli conditions leave a common
stabilizer subgroup of rank at least N-6, allowing a standard Clifford
decoding to at most six active qubits. Decoding and final extraction are
additional polynomial work. This is an audit deduction using known
stabilizer compression, not a new simulator or an executed comparison.

## Validation and disposition

Core passed before scientific edits and after implementation, recorded in
`out/carry_recombination_core_initial.log` and
`out/carry_recombination_core_final.log`. The main experiment passes;
the independent witness experiment exits one for the retained refutation.
Other science suites were not rerun. TODO34 still owns the previously
recorded legacy native crash; no host reliability conclusion follows here.

Documentation is regenerated and checked at integration, with output in
`out/carry_prefix_docs.log`. All four bounded board tasks are reviewed and
closed. TODO46 is complete at this scope; the broader goal remains active.
No manuscripts were changed and nothing was committed or published.
