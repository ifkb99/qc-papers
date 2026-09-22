---
id: 69
state: open
title: "Pauli-LIMDD width of the post-modexp state in orders that interleave exponent and x qubits"
outcome: ""
claims: [C110, C111, C112]
---
# Interleaved orders for Pauli-LIMDDs

**Triage (2026-09-22, METHOD.md barrier check 4):** a small LIMDD width in some order would not by itself factor, because applying the inverse QFT to a LIMDD is not known to stay polynomial. A result that also kept the QFT polynomial would be the TX15 barrier, so any claimed escape must say which step stays polynomial.

From the 2026-09-19 session (note LM). C110 covers exponent-first orders and C111 covers
x-first orders. An order that reads some x qubits before the last exponent qubit is
covered by neither: at that level most sub-states become single basis states, which
all merge.

## Probe (exploratory, not ledger evidence)

`out/agent-board/workers/A4e4043cb93d54b66/order_search_v1.py` (copy of out/limdd_orders/) computes exact LIM-class widths at every level:
XOR-shift classes of the supports, normalised by an element of each set. It runs a
swap local search over orders (4 restarts × 400 swaps, t = ⌈log2 r⌉ + 1). Log:
`out/agent-board/workers/A4e4043cb93d54b66/order_search_v1.log`. The best maximum width it found was 10, 10, 10,
19, 23, 19, 32 and 46 at r = 22, 28, 30, 52, 58, 60, 100 and 106, about 0.32–0.46·r.
Exponent-first gave 15–80 and x-first 12–54 at the same points. No order collapsed.
These sizes cannot separate r from r/poly(n), and local search is not a minimum.

## Questions

* **Last-exponent level in a general order.** C112 owns the exact count there, with F
  the x qubits read before the last exponent qubit. It is small when |F| is large, so
  which other level of such an order carries the cost? A plausible route is to pair this level with the level after the last x qubit,
  as C111 does. Another is the generic bound (classes ≥ distinct sub-states /
  2^(remaining qubits)) combined with an OBDD lower bound for the graph
  χ(e, x) = [x = a^e mod N] at a cut near the bottom.
* Is the √N loss in C110's fibre bound real in any order, or does an averaged
  (collision-count) bound give Θ(r)?
* Can the g > 1 case of C111 be proved for every β, and is T(β) ≤ 2 bitlen(β) in
  general (C111's conjecture)? For odd β < 64, C111's Lemma P plus run 2's shift counts
  already settle g > 1 (β = 5, 9, 17, 21, 33) for every t′ ≥ T(β).
