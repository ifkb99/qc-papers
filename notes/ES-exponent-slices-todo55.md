---
code: ES
title: "EXPONENT SLICES (TODO 55). The Heisenberg peak above a streamed output is removable at fixed N; memory alone does not separate methods."
outcome: record
claims: [C104]
todo: [55, 56, 57]
---
# ES — Exponent slices under a streamed output (TODO 55, 2026-09-16/17)

Chosen from slate round SL under output contract O2 (streamed exact output). This note
records how C104 was reached, what was withdrawn on the way, and the board provenance.
Calibration detail for the review rounds is in SWARM_REVIEW_2026-09-12.md section 5.

## Sequence

1. **Derivation and plan**, task `Ta86889efd4b24683`, three versions. v1
   (`Sc90d08b3801742c1`) proved the slice lemma and a depth-first walk W, and refuted two
   statements in the coordinator's brief (a target-role must-fail control that cannot
   fail; "multiplier-1 blocks are identities", which are the involution V of C23).
   Review `V947329feabd540ad` found that per-slice C45 runs B1 already deliver the
   stream with existing code, and that a planned control could not fail. v2
   (`Sec74297d42184c28`) derived B1 and the dense recursion K2; review
   `Vc95cba7568944751` found the same unfailable-control defect in a new place and a
   byte-ledger slack larger than the pitfalls it had to catch. v3 (`Sc4c99b9dd0c14e6c`),
   with every control instantiated and asserted before submission, was accepted in
   `Vaa6ee566b4454a47` with nine execution-time corrections.
2. **Phase A**, task `Tc70e07bf09c54b3e`. Run 1 (`S5d39eb493d264070`) stopped at gate
   G1 on 26,228 B against a 24,000 B precision band, traced to CPython tuple free-list
   residue; a focused review (`V2d55fc6b31404ef2`) approved a proportionate amendment
   after correcting two defects in the coordinator's amendment text. The rerun
   (`S20e03f947ba94317`, accepted `V7293e51a79fe4ac6`) confirmed exactness at 14 points
   and failed K2's leaf-slot band at F1 (N=7, a=2) t=4..6: robustly at t=5,6 (1,835 and
   5,587 B over) and by 3 B at t=4, inside the ~46 B process-to-process variation.
3. **Bug check** (exploration), task `T2074e4a4ace54998`, `S96857fe4418b49bd`, accepted
   `V6b0fdd469efc4fd4`: the growth was a bounded keyword-dict free-list refill from
   `np.take`'s wrapper, removable by calling the array method.
4. **Fix plan and repeat**, tasks `T49f2eb21610b4bb4` (`S37d7675ac15f4190`, accepted
   `V5df0e47ea3df4a6c`) and `T7ac44bd84996402b` (`S0077516c27884295`, accepted
   `V0e29250446b9478a`): every registered prediction that was evaluated held. The
   F3 (N=7, a=6) parts of the byte predictions at t=7..10 were not measured, so they
   did not hold; the registered CPU projection dropped them (`V0e29250446b9478a`
   IC-2, D(2)). The numbers are in C104. Their evidential weight is the digests, T_L,
   the tight per-t bands and the difference-scored controls, not the wide bands
   (`V0e29250446b9478a` IC-8); EQ-X, A2-X and SRC-X are sanity checks only.
5. **Prior-art survey**, task `T8c5935fbd8814ef2`, `Sef28b547805c419a`: applied known
   mathematics; it derived a careful dense route D*.
6. **D* plan**, task `T0c34d1af337d41b0`. Review `Vbfd8fd978d2f46a0` (on
   `S58deebc6e57a44b8`) showed D*'s second full-size array is avoidable (a windowed
   replay), moving the claimed crossover from t=3 to t=4. Review `Vf1d3e39b30b84bb1`
   (on `S4759788344bf4411`, accept) then showed that no crossover can be claimed: dense
   routes trade passes for memory. The user chose not to run the D* measurement and to
   open TODO 57 on time at matched memory.

## Withdrawn or corrected statements (all the coordinator's unless stated)

* "The leaf-slot excess doubles per exponent bit, a possible O(2^t) term": four points
  of a free-list filling to its cap; withdrawn after the bug check.
* "K2 is 114× below the engine at N=7 t=5": that used K2's phase-R peak; the whole
  traced peak gives 99× (K2/E = 0.010).
* "K2's memory is flat in t": the N=7 constancy is T_L = 2; at N=11 it grows with T_L.
* ru_maxrss values once quoted as MB are KiB, over a ~41 MiB baseline.
* Survey: "K2 is the partial Walsh–Hadamard transform over the exponent register"
  (derivation_v3.md:730) omits the one-representative mechanism; corrected in C104.
* Survey: "K2 is smaller than careful dense from t=3"; the windowed route moves that
  to t=4 and the multi-pass family removes the crossover altogether (step 6).
* C104's first integration draft stated that comparison backwards, saying the dense
  route would be LARGER than K2 at t <= 3, which flattered K2 exactly where it is
  1.5-5x worse. Caught by review `V99e34f8a27704055` (RC-1).
* Note SL labelled the dense recursion "B1"; corrected to K2.
* A linter false positive (`HARNESS-NO-CONTROL` on controls registered in a loop) failed
  `S0077516c27884295`; the linter was fixed and the referee confirmed the controls.

## TODO 57's round, and what it cost to state correctly (2026-09-17)

The matched-memory question was derived (`S1b25ca7326624b0b`, accepted
`V7fb5f7389f5941b8`), executed (`Sa13414787bd54716`, accepted `Va41b64e4d5064682`) and
integrated into C104's Limits. Three defects in that chain are worth remembering, all of
them mechanism claims rather than measurements:

* The plan would have measured the wrong opponent. Design review `Vb02c29cc1e094239`
  (C1) found that it excluded hybrids whose proved floors fit the budget, so its
  registered arm was 4–23x slower than those members by the plan's own rows; a later
  review (`V7fb5f7389f5941b8`, N4) separately observed that the fastest ADMISSIBLE
  member is K2's own endpoint, from the family identity rather than from measurement.
  Both were caught before any compute. The execution then measured those hybrids ABOVE
  the budget, so the excluded members were not in fact admissible.
* A mechanism was asserted from a numerical coincidence: the buffer gap equals the
  plan's registered slack, 49,152 B, by two unrelated derivations. The causal step was
  already in the earlier review's own hypothesis, which labelled itself untested
  (`V751b7d55aeed4816` C1); the coordinator's correction brief then relayed it as that
  review's account, and the executor wrote it down. Nobody ran the refuting check, which
  was four integers from the plan's own ledger (`Va41b64e4d5064682` section 5). The correct account is a MISSING term for the unwindowed transform.
* The executor then under-reported its own sharpest finding: the two counterfactual
  admissibility flips are exactly the two arms that are faster than the reported winner.

Every corrective round in that chain introduced a new defect. The pattern that produced
them is narrow and nameable: a constant was matched rather than decomposed.

## Hashes of the unversioned scientific inputs (2026-09-17)

`out/` is gitignored, so these identities are the only link between C104 and its
evidence: derivation_v3.md 52f9e25188ace435…, k2_proto_v3.py fbca01e4ea9bc97e…,
survey.md a9c6c811207bfb00…, a2p_report_v1.json 9f2fe4b30f49…, c3_report_v1.json
0d17cc4131c5…. Full list: `out/integrate-todo55/scientific_inputs.sha256`.

## What this does and does not answer for the original goal

The Heisenberg peak above a streamed exact output is removable at fixed N: K2's traced
peak is two orders of magnitude below the dictionary engine's at N=7 t=5. That is a
statement against the engine, a weak baseline. Against dense exact streaming no memory
statement is available, because dense routes can make memory small by taking more
passes; the open question is time at matched memory (TODO 57). The output itself stays
at about 2^(q−1) terms (C30); compressing it is TODO 56.
