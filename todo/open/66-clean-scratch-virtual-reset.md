---
id: 66
state: open
title: "Clean-scratch virtual reset: apply C105's clear rule at every point where a unitary circuit returns scratch to |0>, and measure what it removes"
outcome: "Virtual reset pilot is in VR; the charged clean-block sampler pilot is complete (C115, MC) and loses on total resources; work-Rx domain closure remains open"
claims: [C45, C50, C104, C105, C106, C115]
---
# Clean-scratch virtual reset

**Triage (2026-09-22, METHOD.md barrier check 4):** restricts to clean inputs; computes nothing known to be hard. Passes.

From a 2026-09-18 discussion of where PPS memory goes once the exponent register is
handled (C45, C104). What remains exponential is the work width q_w = 3n + 4, and
two claims place the excess on inputs the physical computation never visits:
C50 (clean-equivalent scratch extensions change full Walsh support) and C106
Theorem C (rank beyond the clean orbit is carried only by dirty inputs).

## Idea

C105 proves that a reset of a register pulls back on a diagonal observable as the
clear rule Z^z ↦ Z^(z with that register's bits cleared), colliding keys summed. A
unitary circuit whose scratch is |0> at a known point on every input of interest
behaves, on those inputs, as if a reset were inserted there. So insert the clear
rule at such points in reverse propagation (a "virtual reset"). For off-diagonal
terms the analogue drops every term with X-support on the cleared qubits.

## Derivation to check before running (ToffoliModExp, observable Z_x0)

Registers: b, t (m = n + 1 each), x (n), c0, anc. Read from `toffoli_arith.py`:
* t, c0, anc are |0> at every `cc_add_mod` boundary; b is also |0> at every u_a
  block boundary, **provided x < N** on entry to a block with control 1: the inverse
  cmult subtracts x_orig mod N from b = x_orig, which returns b to 0 only when
  x_orig < N (the modular adder's b < N precondition).
* R0 (clear at the input boundary only) is restriction to the subcube
  S = {b = t = c0 = anc = 0}; its final spectrum is that of f|S exactly (standard
  subcube-restriction identity). Peak unchanged.
* R1 (clear at every block boundary): final spectrum is that of g, where g resets
  the scratch after each block. g = f on S ∩ {x < N}; they can differ only when
  x ≥ N and at least two exponent bits are set (the first sets x < N and leaves b
  dirty; a later block adds into dirty b and the cswap moves it into x).
  With C45 contraction of exponent qubits, support at block boundaries is ≤ 2^n.
* R2: additionally clear t, c0, anc at every cc_add_mod boundary.

## Pilot (bounded, CPU, N ∈ {5, 7}, t ≤ 3)

Independent reference: classical forward evaluation of the gate list on each input
of S (with and without scratch zeroed after each block) followed by a dense WHT. Not
a second propagator. Predictions and must-fail controls are frozen in the
experiment header:
* clear t/c0 in the middle of a Cuccaro MAJ ladder (carries live in t) must
  disagree with f|S on some x < N input;
* R1 against f|S must disagree on some x ≥ N input with ≥ 2 exponent bits set
  (vacuous at t = 1, where no second block exists).

## What would matter

For plain diagonal permutation circuits R1 only rediscovers (e, x) ↦ a^e·x mod N
from the gate list, and C45 already makes the input expectation trivial. The value
is for objects that are not trivial: TODO 64 (one Rx(θ) between blocks, off-diagonal
branch 27–220× the permutation count) and TODO 14 (through the inverse QFT). The
technique is almost certainly standard ancilla-aware simulation; do not claim novelty.

## Checkpoint (2026-09-18, note VR)

Pilot done: `experiments/experiment_scratch_reset.py`, 46/46, both controls failed as
required. The virtual reset is exact where derived and collapses the output, but the
peak sits inside the forward cmult, where b is live, so clears cannot reach it.
**Then proposed (superseded by the checkpoint below):** the block-level split in VR's "Reading": certify each block's clean return
by replaying its 2^(n+1) clean inputs, then transfer the operator on x through that
permutation. Apply it first to TODO 64's Rx fixture (N = 3, 5, 7, t = 2), where the
off-diagonal branch is 27–220× the permutation count, against the exact dense
reference. Needs its own predictions and a must-fail control (a block whose
scratch does not return clean, e.g. with x >= N admitted).

## Checkpoint (2026-09-21, C115 and note MC)

The clean-block transfer now has a retained research implementation and bounded
inverse-QFT sampling comparison. C115 owns the negative total-resource result
and its validation; MC owns the frozen evidence and reproduction pointers.
The tested branch stopped under its declared rule and is not the preferred
backend. Do not repeat its resource sweep without a new discriminating design.

The earlier proposal to apply this directly to TODO 64 was too broad: a work
Rx can leave x<N, invalidating the clean-return promise. Before another
implementation, determine whether the actual insertion preserves an enlarged
certifiable subspace, or which dirty-scratch amplitudes must be retained. A
solution must charge that representation and compare with existing logical
methods. This residual keeps TODO 66 open; TODO 64 remains unresolved. The
separate structured Grover sampling proposal is TODO 71.
