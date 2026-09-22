---
id: 70
state: done
title: "Can arithmetic collision counts force a large MPS bond even when the clean post-modexp state is approximated?"
outcome: "Completed by proof and independent review (2026-09-21): C113 gives a fidelity-robust bond obstruction at the specified mixed cut; no numerical diagnostic needed"
claims: [C19, C107, C108, C109, C110, C111, C112, C113]
---
# Approximate bond from arithmetic collisions

## Completed at the proof stopping point (2026-09-21)

C113 owns the resulting theorem and proof. Two independent author reports and
fresh scientific reviews establish the bound for the normalized physical state
at this exact mixed cut, including its arithmetic hypotheses and fidelity
convention. TX34 records the external transfer and source caveat; note QB
records the round and evidence IDs.

The degree bound and finite-Fourier argument resolved the intended question
without a new fourth-order estimate. The optional diagnostic was not run:
its outcome would not change that decision. Sharper tails or other orders are
separate questions, not requirements left unfinished here.

## Original research contract

The following contract was selected from note QB and authorized by the user.
It is preserved to show the scope and stopping rules of the completed round.

## Question and exact object

Does a specific interleaving of the exponent and work bits still require a large
MPS bond to approximate the **physical clean post-modexp state**?

Let N be odd, N >= 5, gcd(a,N)=1, n=bitlen(N), t=2n and Q=2^t. The normalized
state is

    |Psi> = Q^(-1/2) sum_{0 <= e < Q} |e>|a^e mod N>|0_scratch>.

It is the state in C110, with its normalization supplied, immediately before
the inverse QFT in `ToffoliModExp(N,a,n_exp=t).build_shor()`
(`toffoli_arith.py`). Initialization includes X(x0); do not apply it twice.
Every scratch bit is clean here and may be factored out.

Fix LSB-indexed qubits and the order

    e0, e1, x0, e2, e3, x1, ..., e_(2n-2), e_(2n-1), x_(n-1).

The first target cut is after the first floor(n/2) triples. Both registers are
split at this cut. It avoids claiming that a result at the complete
exponent/work-register cut settles interleaved orders. No optimization over
orders, GF(2) changes of variables or Pauli-LIM equivalences is included.

Accuracy contract: a normalized approximant |Phi> with squared fidelity
|<Psi|Phi>|^2 >= 1-epsilon^2; use epsilon=0.01 as a reporting point and retain
epsilon symbolically. This suffices for total-variation error <= epsilon after
the same exact inverse QFT and measurement. It is a sufficient full-state
contract, not a necessary condition for an output-only sampler. Preparing the
approximation or applying the QFT efficiently is not implied.

## First derivation: count rectangles instead of diagonalizing a state

Across the fixed cut, let B[u,v] be 1 when the concatenation of row assignment u
and column assignment v is a supported bit string (e,a^e mod N), and 0 otherwise.
B has Q ones. Derive from the partial trace, rather than assume from an analogy:

    P_A := Tr(rho_A^2)
         = Q^(-2) sum_{u,u',v,v'} B[u,v] B[u',v] B[u',v'] B[u,v'].

These are ordered rectangles, INCLUDING u=u' or v=v'; counting only distinct
four-cycles would be wrong. Starting from Schmidt coefficients and
Cauchy-Schwarz, re-derive the proposed necessary bond bound

    D >= (1-epsilon^2)^2 / P_A.

The scientific work is to express the rectangle sum as explicit modular-power
and binary-coordinate constraints and bound it for a stated arithmetic family.
Start with prime N to isolate the mathematical issue; a prime-family result
does not establish the corresponding result for Shor's semiprimes. Examine
whether C107's small-kappa hypothesis actually controls this fourth-order
quantity. It is NOT permission to transplant the scalar output-bit theorem.
If it does not suffice, isolate the extra correlation estimate needed and look
for a counterexample before introducing a stronger assumption.

## Acceptance and stopping points

Within one derivation round (suggested cap: two hours of author work plus a
separately bounded fresh review), return:

1. The exact state/cut dictionary, purity identity, accuracy conversion and
   necessary-bond inequality, with hand-checkable product and diagonal-support
   examples. State every normalization and quantifier.
2. The explicit arithmetic rectangle constraints for the chosen interleaving.
   Compare the proposed assumptions with the source bodies behind C107/C109;
   do not replace missing fourth-order control with a fitted rank curve.
3. Either a nontrivial bound for a named growing family, a counterexample to a
   proposed transfer, or the precise unresolved estimate and why available
   results do not establish it. The third outcome is a scoped unresolved
   result, not a theorem of impossibility or a completed simulation algorithm.
4. A decision: pursue a robust lower bound if purity is provably small; inspect
   the full Schmidt tail if the purity bound is weak; stop this transfer if the
   arithmetic assumptions fail. Large purity alone does not prove an efficient
   approximation, and a single good cut does not prove an efficient MPS.

C108 bounds exact rank of a scalar observable, not this state's truncation
error. C110-C112 count exact Pauli-LIMDD classes, not MPS approximation bonds.
The familiar register-cut analysis (C19 and TX9) is the comparison case.
TX15 remains the order-finding barrier; a lower bound respects it, and a weak
lower bound does not evade it.

## Conditional diagnostic in the original proposal (not run)

Only if a numerical answer would change the next mathematical step: n=3..6,
t=2n, CPU only, at most ten minutes and 128 MiB of numerical arrays. Select and
freeze a growing family, bases, arithmetic hypotheses, controls and predictions
before the run; vary size without conflating a base change with a size effect.
Charge orbit enumeration and do not supply a discovered order for free.

Compare an exact integer rectangle count with an independently constructed
reduced-density calculation. For the smallest cases obtain supported strings
by classical replay of the actual gates, not solely from the modular formula
being checked. Use the existing gate representation, not a second propagator.
If a spectrum is inspected numerically, label it diagnostic rather than a
certified tail bound.

Controls must detect omission of degenerate rectangles and wrong normalization;
product and maximally correlated supports provide opposite endpoints. A
high-exact-rank state with a small-weight entangled component, following the
construction in TX34's Schuch source, checks the false inference from exact rank
to approximation cost. Instantiate controls before review. Run the core gate
first if this later scientific execution is authorized.

The completed round used read-only derivation/source assignments, separate
fresh referees and an independently reviewed canonical integration. No numerical
experiment assignment or scientific run was needed. Note QB records the IDs.
