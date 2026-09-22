---
id: TX35
field: "Weighted automata and structured Grover simulation"
status: imported
effect: removes-exponential-in-subfamily
one_line: "A two-state signed counter evaluates selected coherent Grover probabilities for an explicit path predicate; the best structured baseline is the same counter"
source: "Direct predicate recurrence in C114; Stoudenmire-Waintal PRX 14, 041029 (2024), sections IV.B and V.C, equations (13) and (20)"
claims: [C114]
notes: [MC]
todo: [71]
---
# TX35 — Structured Grover probabilities as weighted counts

## Dictionary

The predicate f(x)=x0 AND no adjacent 11 is a two-state path automaton.
A Fourier sign changes its transition weights. The uniform-start Grover
state is a linear combination of the constant function and f; its requested
post-H amplitudes are signed counts. C114 owns the derivation, benchmark
observations, implementation and output contract; note MC records discovery.

## Hypotheses

The source predicate, its local clauses and uniform initial state are given;
solutions and their count are not supplied. Grover oracle and diffusion steps
preserve the two-coefficient form. The final H layer is evaluated coherently.
This does not account for intermediate gates inside a reversible oracle.

[Stoudenmire and Waintal](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.14.041029),
section IV.B, equation (13), expresses a post-oracle state as a sum of product
states indexed by solutions. That small-solution route is not imported here:
the path predicate has exponentially many solutions. Section V.C, equation
(20), instead bounds construction for fixed local block size and interblock
depth, with the intrablock clause count depending only on fixed block size. The body was checked; our particular two-state recurrence is derived
directly from the predicate, rather than asserted from that bound.

## Consequence for the goal

The exponential state vector is unnecessary for this restricted query task.
The strongest direct weighted-count method is identical to the implementation;
no improvement over it, generic black-box oracle claim, search speedup or new
simulation principle follows. Exact integer/rational bit cost is charged.
Full output sampling remains a separate contract in TODO 71.
