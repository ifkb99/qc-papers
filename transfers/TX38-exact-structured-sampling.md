---
id: TX38
field: "Matrix product states and exact Born sampling"
status: imported
effect: removes-exponential-in-subfamily
one_line: "Known sequential MPS sampling applies to the explicit path-predicate Grover family; its homogeneous suffix environments admit an exact rolling Fibonacci recurrence"
source: "Ferris and Vidal, Perfect Sampling with Unitary Tensor Networks (2012), section IV, equations (21)-(28); Stoudenmire and Waintal, PRX 14, 041029 (2024), section IV.A and footnote 40; direct recurrence in C117"
claims: [C117]
notes: [GS]
todo: [71]
---
# TX38 — Exact sampling from the structured Grover representation

## Dictionary

C117's constant channel plus two-state predicate counter is an open-boundary
bond-three MPS in natural bit order. Squared-amplitude prefix contractions
give the conditional probabilities for sequentially sampling all output bits.
The local path predicate is supplied; marked strings and their count are not.
C117 owns the recurrence, bit-cost bounds, implementation and measured scope.

## Hypotheses

The initial state is uniform, the oracle is x0 AND no adjacent 11, and the
final layer consists of local Hadamards. The state is pure, real and exactly
represented by the constructed small-bond network. The output contract is a
finite collection of sampled outcomes, not the full probability vector or the
peak inside a native reversible implementation of the oracle.

[Ferris and Vidal](https://arxiv.org/pdf/1201.3974), section IV, equations
(21)-(28), derives sequential Born sampling; footnote 32 supplies the
open-boundary canonical-form connection. Its optimized network is already
given. [Stoudenmire and Waintal](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.041029),
section IV.A and footnote 40, distinguishes recursive sampling from initial
orthogonalization. Neither source makes construction or arithmetic precision
free, nor specifies the exact random-bit procedure used here.

The implementation uses integer noncanonical suffix environments, avoiding
an assumption that canonical tensor entries must be rational. It constructs
the count and coefficients, and charges all environments, random draws and
retained outputs. Its Fibonacci environment is a direct specialization of the
known contraction to this homogeneous predicate, not a new sampling principle.

## Consequence for the goal

At fixed iteration count, output sampling need not retain exponentially many
amplitudes in this subfamily. Rolling suffix recurrence reduces storage versus
the conventional cached contraction, which remains the same structured method.
Dense-array size is a theoretical comparison at 128; it is not a feasible
benchmark or the strongest structured baseline. Increasing iterations toward
Grover amplification changes the exact-arithmetic costs. This transfer gives
no generic oracle, Shor, native-gate peak or universal exact-sampling bound.
