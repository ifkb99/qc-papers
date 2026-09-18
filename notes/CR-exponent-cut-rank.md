---
code: CR
title: "TODO 63. The exponent-prefix cut rank: a lower bound for linear representations, exact for commuting letters."
outcome: record
claims: [C106, C103, C48]
todo: [63]
---
# CR — The exponent-cut rank (TODO 63, 2026-09-18)

**Origin.** Surveyor candidate R1 in slate round 2 (note SN). The slate's prior: a few
dozen distinct pseudo-random ±1 rows over thousands of columns are almost surely
independent, so a measured rank equal to |S_k| would decide little, and the content
had to be a proof.

**Derivation.** Task `T095a76c1015e465e` went to a fresh deriver on the best available
model (`S1e0938302b39406a`). It proved Theorems A–C. The prior turned out to be wrong
in the way that matters: for β = 3, α = 0, where the letters commute, the rows are
consecutive Krylov vectors T^i s, so they are independent only up to the Krylov
dimension d. A fresh referee
accepted it (`V13f2bfb187274d58`). It recounted d with its own Möbius-product code
(40 random permutations in both orientations, 80/80 against exact Krylov rank) and reproduced the period. It
required a scale-matched must-fail control, C5 (without it, a rank routine capped at
1936 would pass P2b and P2c), and made I3–I5 record-only, since they test a statement
that was not reached.

**Experiment.** Task `T3e37f92ebd754215`, run by a falsifier (`Se701e050f6a04e80`,
accepted `V3cec8224aec7454d`). 16/16 harness checks passed and every must-fail control
failed as required. P2c, the one prediction against the prior (1937 rows → rank 1936),
is gated in code on C5 failing. Beyond the plan, the falsifier registered P2e, an exact
minimal polynomial over Z, and a bug check of it. The first attempt (`A9b6f006107614633`,
run `Rd2c648bf30414b27`) lost its lease during a usage-limit cutoff after the run had
finished. The coordinator reclaimed it and linked the late work (`Mbee26551788544fa`),
and the same script was re-run in a new attempt: 102 log lines, 0 differing readings.

**Process notes.**
* The plan's "under five minutes" was exceeded: the main run took 7.2 min, within the
  30-min cap. From the run log's own stamps, most of it was I2's exact eliminations
  (raw 1937-row elimination 99.5 s, C5 86.0 s, the two Gram ranks 39.4 s and 39.6 s)
  and the circuit replays (I4's alone 62.8 s). The I3–I5 eliminations took under a
  second each. The experiment review's attribution to the I3–I5 eliminations was wrong
  (integration review `V58130d845fce41c4` D2).
* The board refused `run.start` on a prediction message from the expired attempt
  (`invalid_prediction`); a verbatim reaffirmation (`Mfd70c990cb894be7`) was required.
* The derivation's v1 exploration had a Berlekamp–Massey value of 4098, wrong through
  int64 overflow. It was caught because it exceeded the 4096 ceiling.
