---
id: 15c
state: done
title: "3/4: conditional measurement branches for actual order-finding output"
claims: [C45, C46, C49]
outcome: "Conditional pilot validated against full quantum circuits and wide-output formulas; the noncommuting-reordering control fixes its scope"
---
# 3/4 — Conditional order-finding branches

Completed pilot 2026-09-10: C49 and §CF. Dense work storage/setup remain
exponential. TODO 14 is still open for a useful compressed inverse-QFT model.

After 15b, test a semiclassical terminating inverse QFT using conditional
measurement updates, preserving interference. Compare the entire small output
distribution with existing state-vector circuits and with an independent DFT
reference, including non-power-of-two orders and bit-order controls. Measure
the branch representation, not just the averaged observable. Do not call a
method efficient if it first enumerates a large orbit or all output histories.

Primary starting points: Griffiths–Niu, arXiv:quant-ph/9511007; Browne,
arXiv:quant-ph/0612021. This is a bounded pilot for the broader TODO 14, not
automatic closure of that frontier. Then proceed to 15d.
