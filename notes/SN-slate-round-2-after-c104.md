---
code: SN
title: "SLATE ROUND 2. After C104: the unlookup gap, the exponent-cut rank, and TODO 56 closed on its contract."
outcome: record
claims: [C52, C99, C103, C104]
todo: [56, 58, 62, 63, 64]
---
# SN — Slate round 2 (2026-09-17/18): the next direction after C104

Second use of SWARM.md phase 0. The thread SL → TODO 55 (C104) → TODO 57 had left
TODOs 56, 58 and 59 behind it. Two fresh contexts, which could not see each other, ran
on the best available model: a surveyor (transfer, obstruction; `S900eeb6398644d93`)
and a deriver (structure inventory, symmetries, sequences, nine exact CPU runs;
`Saf6eaffcf3a0462c`). Before dispatch the coordinator sealed its own candidate by hash
(`M7763a806b65a4ca8`).

**The trigger was weaker than stated.** TODOs 58 and 59 are reproducibility items, so
the thread had one scientific residual (TODO 56), not three (surveyor, §0.2). The slate
still changed the ranking.

## Premise corrections (surveyor, §0.2–0.4)

* The brief said K2 is "the fastest admissible member ... measured". C104 Limits
  measures the non-endpoint slowdown only; that the endpoint is admissible is derived.
* The Ω(2^q) O1/O2 floor is an **emission** floor of 2^(q−1) items, not an arithmetic
  floor. For β ∈ {1, 3} the O2 arithmetic is polynomial in t below saturation, because
  T_M and T_{M^−1} commute (C101), which collapses K2's leaves to 2^α(|E|+1)(|O|+1)
  distinct spectra. This is derived and unmeasured. (Corrected 2026-09-18: an earlier
  version said this changes C104's wording. C104 never states an Ω(2^q) floor; the
  sentence the surveyor corrected was in the coordinator's slate brief. Detector: a
  search of claims/C104.md for "floor", "Omega(2^q)" and "Ω(2" finds only "its proved
  array floor", which is a different statement.)
* SL's O3 has no query axis, and without one it is vacuous: the circuit is itself a
  T_L·2^(q_w)-word exact representation with O(t·2^(q_w)) coefficient queries. O3 is
  refined to (S_store, T_build, T_query).

## Candidates and choice

The user chose **R4 and R1**. C1 is queued behind a premise check.

| id | candidate | origin | disposition |
|---|---|---|---|
| R4 | Gidney's measurement-based unlookup (X-measure, diagonal fixup, reset) has a Heisenberg adjoint that clears the lookup register's key bits: no branching, and diagonal stays diagonal | surveyor | **TODO 62** |
| R1 | the rank of the ±1 matrix F_k across the exponent-prefix cut lower-bounds the bond of every representation linear across that cut, and is at most C103's \|S_k\| | surveyor | **TODO 63** |
| C1 | one Rx(θ) between blocks: exact support is the same for every θ with irrational cos θ; the off-diagonal branch pays 27×/110×/220× at N = 3, 5, 7 (t = 2) | deriver | **TODO 64**, queued |
| R2 / C2 | TODO 56's class form generalised to every β: W(z_e, z_w) = Σ_{g∈S_t} κ_g(z_e) Ĝ_g(z_w), κ from the prefix automaton | **both, independently** | TODO 56 **closed**, below |
| C3 | distinct K2 multipliers T_L(t) = min(t, α + ν(β)), ν = ord_β(2), proved; verified on 17 families | deriver | into **TODO 58** (predicts 5 tables at F2 t = 7) |
| R3 | β ≥ 5: count K2's distinct exponent-Walsh leaves | surveyor | parked; bounded below by R1's rank |
| C5 | prove C103's equality through witness states (59–99% of dirty states separate all reachable words at t ≤ 8) | deriver | parked |
| C4 | a time–space lower bound for O2 (pebbling) | deriver | parked (the surveyor discarded the Hong–Kung route: K2 escapes it by changing the DAG) |
| R5 | a residue-automaton bound on C102's u_a ROBDD (growing N) | surveyor | parked |
| C0 | the post-QFT distribution for the real input as the output contract | coordinator, sealed (sha256 c3eccf6a…c30c4) | **killed**: C52 already gives scalar conditional sampling with no work vector once r is known, and C45 makes the input expectation trivial. Neither slate generated it, and both discarded the neighbourhood. |

## TODO 56, closed

Both slates derived the same identity independently, and both found its value
conditional. (i) Nothing in the ledger consumes Walsh-coefficient queries: the pre-QFT
expectation is the z_e = 0 slice (C45), and post-QFT sampling needs amplitudes.
(ii) C103's exponent-first ROBDD is already an O3 object answering queries in O(nodes),
and it is smaller at the one benchmark point: 8,517 nodes against 58,279 class terms at
(7, 2, t = 16). (iii) For β ≥ 5 in reach, |S_t| = 2^t for t ≤ α + μ + 1, and at
N = 11, t = 12 the compression is 3.4× (arithmetic from C103's width law). The identity
was instantiated at (11, 2, 4) and (7, 2, 4): the automaton κ equals brute-force class
sums, and the reconstruction equals the direct transform in total count (521,932;
62,330) and at 64 random keys (deriver s4). Whether the class form is minimal among
linear representations is TODO 63.

## Deriver's symmetry and seed results

* **Sign law.** W_{a^−1}(z) = (−1)^{|z_e|} W_a(z) holds at t = 2 and 4 and fails at
  t = 1 and 3 for N = 7, 9, 13. So (a, a^−1) at even t are one fixture, not two.
* **SL seed "N = 9, a = 4 / 7 equal counts at t = 3".** Refuted as a symmetry: the
  per-class supports differ (15,404 vs 15,414), the per-z_e count vectors are different
  multisets with equal sums, and the same pair at N = 7 and 13 has unequal totals. It
  is a coincidence of totals. The other SL seed (the 3/4 cover) is O1-only and was not
  reopened.
* **The dirty-space block group is not Z/r.** M_c^3 ≠ I moves 7,836 of 8,192 states at
  N = 7, and M_{c^2} ≠ M_c^2, while both hold on the ideal orbit.
* The deriver's own scripts had three defects, each caught by its own registered
  cross-check and preserved as v1/v2: a uint8 sign overflow, the ideal-orbit set fed
  through iota, and π/3 called generic.

## Process lessons

* **Convergence again.** As in round 1, the two generators converged on a candidate
  (TODO 56's generalisation), this time together with its obstruction.
* **Reading earlier slates costs independence.** Both workers read the round-1 slates,
  which the brief permitted, and disclosed it. R2/R3 and C2 are therefore not
  independent of the round-1 deriver.
* **The coordinator's sealed candidate was already answered by the ledger.** It was
  worth sealing, since it measured the coordinator against two fresh generators, and it
  lost.
* **Review statements carried into the ledger, three times.** In the integrations
  of TODOs 62 and 63 the coordinator transcribed review statements without
  re-deriving them: "support sizes differ" (review `Vab1fa0e744c8487c` I2), a float
  rank of 1913 (`V58130d845fce41c4` D1) and a timing attribution that the run log
  shows was wrong (D2). The referees caught all three. A mechanical detector was
  tried afterwards on the frozen v1 submissions (`S776d9c75bb0d4899`,
  `S688bd467869a4ed6`) and was not adopted. The phrase-overlap variant is archived
  with its log in attempt `Ac15c8fb2825b4168` (`detector_trial.py`,
  `detector_trial_n5.log`); the token-location run was not archived, and its result
  below rests on inspection. Locating each number by token in the
  sources' frozen evidence flagged none of them: 1913 is in the falsifier's own
  frozen results file, copied there from the binding review, and the other two are
  wording. Five-word phrase overlap with review bodies but not evidence flagged the
  1913 line and the timing paragraph, but it missed "support sizes differ", and
  it flagged 18 other lines, mostly binding corrections and shared phrasing. The rule is
  now a list, in SWARM.md: the integration summary names each review-sourced
  statement and what the author re-read.
* **A premise check missed by the deriver.** C1 sits beside the perturbed-circuit
  claims C52–C57 (single-defect output effects, including Rx(θ) kicks) and cites none of
  them. TODO 64 starts with that check.
