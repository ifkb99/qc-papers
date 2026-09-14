---
id: 35
state: done
title: "Can additive physical phases expose a useful boundary through Gauss sums?"
outcome: "C76/GS prove full prime-sector coupling and validate physical kernels and output; alternative contraction remains open in TODO 36"
claims: [C52, C53, C56, C58, C59, C74, C75, C76]
---
# A cheap physical gate with dense orbit-Fourier mixing

**Completed at the frozen first-pass scope.** C76 owns the exact result and
scope, GS the physical/algebraic/full-output reports and corrected failures.
Steps 1–3 below have been executed and audited. Step 4 is a new investigation,
now tracked ONLY in TODO 36; no general contraction/sampler was obtained.
The earlier plan and pilot below are preserved as pre-experiment history,
not current pending work or a substitute for GS's released evidence.

C75 closes the clean repeated-W/reflection construction at bounded scope,
but that q=0 family retains C59's small static baseline. Do not polish a
large truth-table compiler to seek a speedup already excluded by that
baseline. Instead test a different physical intervention:

    G_k |x> = exp(2*pi*i*k*x/N) |x>.

In binary physical encoding this is a product of n=ceil(log2 N) one-qubit
phase rotations, plus an explicitly tracked global phase. It does not need
a classical x-to-discrete-log evaluator. Angle-description and synthesis
precision costs are additional. This is an additive residue character,
NOT C58's multiplicative cell character D_q.

For prime N and a primitive root, let r=N-1. On the full multiplicative
orbit (b=1), the Fourier coefficients are

    h_k(delta) = (1/r) sum_(j=0)^(r-1)
        exp(2*pi*i*k*a^j/N) exp(2*pi*i*delta*j/r).

These are normalized Gauss sums. For k!=0 mod N, the zero-frequency
squared magnitude is 1/r^2 and the others are N/r^2: dense support with
simple probabilities but nontrivial coherent phases. For k=0 only the zero
frequency survives. This is known mathematics, not a new project theorem;
[van Dam and Seroussi, finite-field definitions and facts](https://arxiv.org/pdf/quant-ph/0207131)
also distinguish easy norms from phase estimation. Main read those passages
and Lemma 2's discrete-log reduction proof, not the entire quantum algorithm
or finite-ring sections. C74 now credits that closer powering antecedent.

## Historical pilot, before the completed output experiment

A lower-cost read-only pilot checked (N,a)=(7,3),(13,2),(17,3), with k=1
and the k=0 control, within a requested 50,000 summand budget. It reported
magnitude-law error below 7.5e-16 and the expected null control. No standalone
harness report was written. This is an initial sanity check of the known
identity, not an observed output advantage, a hardness result, or a substitute
for the main independently reproducible experiment below.

## Frozen pre-experiment question, now tested at first-pass scope

Use the existing TEMPLATE/harness and Circuit/statevec; lower-cost initial
testers, main proof/code audit. Core first. No long timing, new propagator,
manuscript expansion or hidden high-order character oracle.

1. Freeze the three pilot primes/k values above as the elementary b=1
   formula audit; emit G_k with actual binary one-qubit phases and verify
   every complex column, including its global phase. Independently verify
   the Gauss coefficients and G_k G_l=G_(k+l), especially G_k G_-k=I.
   Removing coefficient phases must fail this composition test. No norm-
   only Markov transition rule is licensed by the flat magnitudes.
2. For actual output, freeze C75's N=13,a=2,b=3,r=12,t=3 initial/post-control
   fine-mixer schedule, REMOVE its q=0 reflection insertions, and insert
   G_k after the second post-control W. Vary only k=0,1,2. Compare the full
   eight-output law against the existing full-r sequential_path reference.
   Test zero kick and end insertion as nulls, and omission/initial coarse
   dephasing as candidate-failure controls. If the fixed interior kick is
   invisible, retain that result; do not silently select a different fixture.
3. Derive the b=3 subgroup kernel explicitly before interpreting its support.
   The b=1 flat-magnitude formula does NOT automatically hold in each of
   these smaller coarse sectors. Check closure of any surviving static
   grouping; compare C53's input-specific exponent-phase transfer and the
   ideal-output approximation, preserving original chronological order.
4. If phases affect output and the old cheap baselines fail, ask whether
   additive coordinates or character-sum identities supply a different
   exact/approximate contraction. Dense sector support alone proves no
   hardness and gives no new sampler. Do not implement a large sector
   matrix, assume Gauss phases are free, or import the paper's oracle
   hardness as a lower bound for this fixed initial state and observable.

Before execution, preflight all retained reference/state arrays (16 MiB),
at most ten physical qubits, per-circuit 30,000 Pauli rotations and explicit
cumulative gate-entry/sum budgets. The C75 compiler is only a tiny fixture;
use its exact conjugation/projector cancellations where appropriate. No
native-crash stress test or new full sampling benchmark is authorized by
this research direction; TODO 34 remains separate.
