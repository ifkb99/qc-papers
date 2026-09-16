---
code: GB
date: 2026-09-14
title: "Group-covariant differential counts and the clock-shift observable"
outcome: identity verified in four order-16 groups; the XOR pairing is exponent-2 specific; the order-8 extremal prediction was refuted and explained for R_a by a telescoping obstruction; on one modexp fixture the Pauli/mixed ratio is 0.64–1.88, in both directions
claims: [C100]
todo: [52]
---
# GB — Replacing XOR by other groups

C100 owns the statements. This note records the path to them.

## Why this was tried

The user asked whether groups or geometry could help beyond C99's XOR
result. C99's derivation uses only group translation plus one XOR-specific
fact: σ is an involution. The post-QFT measurement in order finding is built
from Z/2^t clock shifts. So the natural questions were what survives the
change of group, and whether a basis matched to the exponent's own group
shrinks the actual observable. C89's translated-interval exchange showed
that the project's strongest diagonal constructions already rely on Z/2^m
geometry.

## Runs

A helper smoke test (one random 8-element permutation, three groups)
preceded the header; it measured no prediction. A later feasibility probe
counted Part B classes (R = 256/767 at a = 1, 256/256 at a = 4) and timed
one Fourier transform, but no term counts. PB3 was already declared open,
so no prediction depended on that probe.
`experiments/experiment_group_bridge.py` run 1 then scored 8/9, exit 1;
log `out/differential_bridge/group_run1.log`.

## The refuted prediction (PA4)

I predicted exhaustive maxima of 56 for Z/8 and Z/2×Z/4. My argument
assumed some permutation reaches R_a = 7 with a full-support double class.
The measured maxima were 55 and 48, with max R_a 7 and 6.

Reading the construction again: the differences telescope to zero around
the row, so one repeated difference must equal −Σ_G g. That sum is nonzero
exactly when G has exactly one involution; it is 0 when G has several
involutions or none (odd order). So Z/2×Z/4 cannot reach R_a = 7. Z/8 can, and does. For a double class {y1, y2}, its support is full
unless −χ(π(y2)−π(y1)) is an o(y2−y1)-th root of unity; when y2 − y1 has
order 8 that always happens, which costs exactly one character. This does
NOT mean that an R_a = 7 row forces an order-8 difference. The review
Vefabe53d66914d24 found 7,168 of the 11,264 such rows with an even
difference. The coordinator's archived check
(`out/agent-board/workers/A869ad224272b43d0/z8_r7_double_class.log`) gives
best counts 52, 54 and 55 for difference orders 2, 4 and 8. So the measured
fact is that no R_a = 7 row has a full-support doubled class. Why the
root-of-unity condition always holds has no conceptual proof. Because the
enumeration covers every permutation, direction and character, 55 is
nonetheless established for Z/8 by complete enumeration (the status
relabelled under task T393699c97b9745ca). I first dismissed the Hall–Paige-type sum
argument because c = 0 is excluded; that was wrong, since the argument bites
on the repeated difference.

## Part B reading

On ToffoliModExp(7, 3, n_exp=3), the counts favor the mixed basis for
a ≠ 4 (ratios 1.4 to 1.9) and the Pauli basis at a = 4 (0.64), where the
shift is an XOR on the top bit. The a ↔ 8−a symmetry of the rows is
expected, since M_{−a} = M_aᵀ; it was noticed after the run and is not
evidence. The base is non-degenerate (ord_7(3) = 6), but n_exp = 3 and one base
make this a single fixture. Both bases carry millions of terms out of 4^16.

## Candidate direction (not opened)

On this fixture, term counts stay dense in both bases: the average class
needs about 7.5k–32k characters (T/R). This is not general. A group matched
to the map can collapse counts completely (PA2: T = 1), and the uncertainty
bound forces density only for small classes. The compression the project has actually achieved
(C89 intervals, C93 cells) comes from the *shape* of the difference classes
in their native geometry. A measurable follow-up: count geometric atoms
(maximal intervals in Z/2^t, subcubes in GF(2)^n) of the Part B classes and
compare their growth in n_exp with term counts. Define the prediction and a
must-fail control before measuring.

## Review corrections (Vefabe53d66914d24)

The first referee review (changes_requested) confirmed the identity, the
bounds, the telescoping proof, the preserved PA4 refutation and the Part B
object. It required four corrections:

* **R1:** make the HANDOFF board state accurate.
* **R2:** replace "constant factors" with the measured ratio range.
* **R3:** limit the density statement to this fixture.
* **R4:** restate the Z/8 question.

The non-blocking notes were also applied:
* odd-order groups are obstructed as well;
* the obstruction explains the maximum of 48;
* cite ord_7(3) = 6 for non-degeneracy.

Two points are recorded here:
* PB3 cannot fail, so the run's "8/9 pass" includes a check that only
  records data.
* The Part B feasibility probe is disclosed in this note, not in the
  experiment header, which is unchanged because it is frozen evidence.

## Second review corrections (Vfb370f16c33e447a)

The focused re-review of version 2 (changes_requested) confirmed R1–R4 and
the shadow check. It required three one-sentence fixes, all applied:

* **RC1:** the Z/2×Z/4 bound is R_a·|G| ≤ 48, not a singleton count.
* **RC2:** the zero-sum condition above now covers odd-order groups.
* **RC3:** the HANDOFF checkpoint had described C99's version-2 submission
  before that submission existed; it now states only settled facts.

Non-blocking notes from that review:
* The shadow check has no must-fail control. The referee's independent exact
  check supplies one: a mutant counter shifted to e+1 reaches 56 on exactly
  the 7,168 order-2/4 rows.
* The referee recorded an unreviewed case-analysis sketch suggesting that
  T = 56 − 8/o(d) holds for every χ and is provable. That conceptual proof
  remains a lead; the finite Z/8 statement itself is established by
  complete enumeration.
* C100 now gives the shadow-check path.

