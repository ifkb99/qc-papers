---
id: 15d
state: done
title: "4/4: vary invalid-scratch action while preserving the logical computation"
claims: [C8, C45, C50]
outcome: "Same-layout clean-equivalent guards change full support; a converse control preserves support count but changes odd-order output statistics"
---
# 4/4 — Operationally irrelevant scratch action

Completed first pass 2026-09-10: C50 and §SE. A general quotient construction
or extension optimizer is not claimed.

After 15c, construct same-layout reversible modifications that act identically
on the intended clean-input subspace but differ off it. Verify the logical
isometry (including coherent phases where applicable), then compare full
Walsh support, clean-subspace reduced observables and actual order-finding
output statistics. Include a deliberately invalid modification that changes
the legitimate computation. Avoid the original Fourier/Toffoli layout confound.

The intended payoff is a controlled separation between representation cost
and operational prediction, not a universal impossibility or hardness claim.
