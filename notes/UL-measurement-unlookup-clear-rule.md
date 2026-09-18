---
code: UL
title: "TODO 62. Measurement-based unlookup of the lookup register under diagonal PPS: the clear rule."
outcome: record
claims: [C105, C37, C38]
todo: [62]
---
# UL — Measurement-based unlookup (TODO 62, 2026-09-18)

**Origin.** Surveyor candidate R4 in slate round 2 (note SN), aimed at the gap that WD
and Paper B §7.3/§13 item 3 had named "the one real gap": the windowed analysis
assumed a unitary unlookup.

**Round.** Derivation task `T4132ea2a6c784b19` was written by a fresh deriver on the
best available model and submitted as `S3e2b424794764fa2`. A fresh referee accepted it
(`Ve11cc8f5fd444c0a`) with no blocking findings. The author and the referee both judged
that no run was needed: the derivation answers the TODO, and the conditional plan
(engine reset rule, P1–P5, control F1, mutant M-drop) runs only if channel-model
support sizes are to appear in the paper. In that case the referee requires a separate
`kind:"implementation"` task for the engine change, a frame-path mutant, and M-drop
instantiated at w = 2.

**Corrections to the slate sketch.**
* The surveyor's "deallocate ⇒ drop" was a mislabel. Deallocation gives the clear rule;
  drop is measure-and-leave.
* The surveyor expected the channel model to be "if anything, simpler". Exploration
  found the channel support larger at two of the three WD instances (C/U = 1.69 at N = 5,
  a = 2, w = 2; 1.83 at N = 7, a = 6, w = 1; 0.86 at N = 5, a = 4, w = 1; K = 0). These are exploration
  numbers and enter no claim.

**Referee corrections applied in C105.** D1: the tail triple for bit i is
s ^= (2^i mod N), not X(s[0]). D2: Theorem 2 includes act·V_0. D3: the paper text is
scoped to the lookup-register unlookup, and the support values are not reported.
D4: 2.8 is derived under C38's hypothesis. The referee's own diagnostics
(`out/agent-board/reviews/referee-derive-62/`) confirmed W² = id on the full 19-qubit
space at (5,4,w=1) and (7,6,w=1), and the witness input y = 8704 with the source's own
`replay()`.

**Observations (exploration, not claims).** The z_s = 0 slices of the two final
spectra differed at every instance tried (the derivation listed this under a proved
heading; it rests on exploration only, integration review `Vab1fa0e744c8487c` I1).
The s-marginal count at K = 1 is 7,714 at both (5,4) and (7,6),
although the full supports differ. This is an observation only.
