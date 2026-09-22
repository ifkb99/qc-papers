---
id: TX42
field: "Constructive factor-state sampling and Fourier amplitude access"
status: barrier
effect: barrier
one_line: "Efficient N-only factor-state sampling or conditional counting is already factoring; Fourier reconstruction must charge the arithmetic amplitude routine"
source: "Direct reductions in C119; Stoudenmire and Waintal (2024), IV.B equation (13), VI.C equations (25)-(28)"
claims: [C117, C119]
notes: [FC]
todo: [72]
---
# TX42 — The construction and access bill

## Dictionary

A classical constructor and sampler for the bounded factor state returns
factors. Exact branch-count access also finds them. A local Hadamard basis
change spreads a sparse marked signal but asks for signed global counts,
not just cheap pointwise evaluation. C119 owns these reductions and the
extraction-accuracy example.

## Hypotheses

Input is N and explicit factor bounds only, with a nonempty relation promise.
Construction, tensor entries, normalization, accuracy and extraction must all
be charged. [Stoudenmire and Waintal](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.041029),
IV.B equation (13), separates oracle-state construction from subtraction and
sampling; VI.C equations (25)-(28) discusses Hadamard-basis amplitude access.
The amplitude routine is an algorithmic input that must be supplied here.
C117 supplies such constructive access for its local path predicate only.

## Consequence for the goal

This is an acceptance criterion, not an impossibility theorem. Retain Fourier
reconstruction only if a new N-only arithmetic count formula is supplied.
Generic whole-state approximation is insufficient unless its error remains
controlled after isolating and normalizing the rare solution component.
