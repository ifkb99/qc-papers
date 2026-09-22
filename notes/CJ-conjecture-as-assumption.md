---
code: CJ
date: 2026-09-22
title: "C124's conjecture is an assumption: reduction to REA_fact, and a parked survey"
outcome: record
claims: [C124, C125]
todo: [74]
---
# CJ — Working on C124's open conjecture (2026-09-22)

## Why and how

The user asked to work on C124's open conjecture. It says that no N-only subgroup
with cheap membership and phases has |<a> cap H| super-polynomially larger than
the B-smooth part of r, for typical RSA N. Phase 0 followed SWARM. The coordinator sealed a candidate
(SHA256 68ff3ce3..., posted in M3a47c1cec8894393, never on the shared filesystem)
and revealed it after both slates froze (out/memory74c/coordinator_reveal_v1.md,
hash verified). A fresh deriver and a fresh surveyor worked independently from
the frozen contract out/memory74c/contract_v1.md. Task inputs were only claims
and that contract, and the TODO/transfer text went in source_ref. This avoided
note GR's stale-input cycle.

## Outcome

* **C125 (deriver, accepted).** Under weak phases the conjecture implies
  factoring is not in BPP. Under point phases it implies that factoring and
  F_p^* discrete logs are not both in BPP. Either way it cannot be proved
  unconditionally.
  - For typical N it concerns |H|_R alone.
  - For factored-multiple phases it reduces to the new named assumption
    REA_fact: equivalent under weak phases, sufficient under point phases.
  - Its oracle and all-N strengthenings fail conditionally.
  - Lemma P is the central open gap.

  **The deriver's framing is the decisive one.** The sealed candidate and the
  surveyor both found the factored-order split-or-balanced dichotomy, but neither
  saw that no unconditional proof can exist.
* **Comparison with the sealed candidate.** The candidate's dichotomy (saving
  <= gcd(p-1,q-1) given a basis with factored orders) holds only for such a held
  basis. As a comparator for C124's setting it is wrong: C124's G[k] saves
  gcd(r,k), which can exceed gcd(p-1,q-1) (N=35, a=2, k=6), and G[k] supplies no
  factored basis. Review V40b63e1546734ad5 found this in the survey, which had
  independently made the same claim. The candidate made it too.

## Parked survey (task Ta19ae1ee98b24b1f, cancelled at the user's direction)

Four revisions (S09a0152534d74e2c through Sc78326c4b7914d83) were reviewed in
V40b63e1546734ad5, V993dde7d6a7f44d2, Vdf7bff3ef0324e65 and V0909a684370d49e0.
Each review verified the mathematics. Each wording round introduced new
defects. Some were sentences copied from the previous review's wording without
re-derivation, and one review miscounted a control's coverage (21 against the
correct 9: Vdf7bff3ef0324e65, corrected in V0909a684370d49e0 N17). The user chose to stop.

The results below are **leads, not claims**. Each was re-derived by the review
named. None is in the ledger:

| lead | re-derived in |
|---|---|
| Split-or-balanced for a held pair (g, m) with m = ord(g) factored: gcd(g^(m/l)-1, N) splits N unless ord_p g and ord_q g agree at every l (Seres-Burcsi 2020/402 Theorem 4 mechanism; GLMS arXiv 1511.04385 section IV.B) | V40b63e1546734ad5 |
| gcd(m, N-1) divides D=gcd(p-1,q-1) for every m dividing lambda, and gcd(lambda, N-1)=D | V993dde7d6a7f44d2 |
| For balanced N (p<2q), delta = p+q-2 sqrt(N) < (3/sqrt 2 - 2) sqrt(N), about 0.121320 sqrt(N). A known element order m divides phi, so p+q = N+1 (mod m). Whether the survey's rule L5a recovers p+q depends on m against that rule's own ambiguity window, as checked in V0909a684370d49e0; no simpler sufficient condition is recorded here | V40b63e1546734ad5, V0909a684370d49e0 |
| p+q is known mod c^2 from (N-1)/c for c dividing D. This fixes p+q when c^2 > delta, which holds for p<2q once c > sqrt(3/sqrt 2 - 2) N^(1/4), about 0.348311 N^(1/4); re-derived here, and "0.3483" is the unsafe rounding (counterexample in V0909a684370d49e0) | V993dde7d6a7f44d2, V0909a684370d49e0 |
| Under factoring hardness, an H with point (H3) and polynomial (H2) has super-polynomial index on typical N (V993dde7d6a7f44d2); and delta/ord(h) is super-polynomial for every polynomially produced h in H (Vdf7bff3ef0324e65) | V993dde7d6a7f44d2, Vdf7bff3ef0324e65 |
| \|G[N-1]\| = D^2; Pr[D >= y] <= 16 (ln X)^2/(y-1) + o(1) for uniform primes in [X,2X] | V40b63e1546734ad5 |

To promote any lead, give it a claim of its own and a fresh review, and derive
the text from the construction. Do not copy it from the survey or its reviews.

## Edits made and deferred

* TX50 gains C125 and one sentence. Note GR's deferred O1/O2 wording fixes are
  applied in the same edit.
* TODO74's TX50 paragraph records the reduction. GR's deferred N2 is applied: C124
  proves sufficiency of d <= poly(log N), not necessity.
* **Deferred:** C124's conjecture paragraph needs a pointer to C125. C124 is a
  sealed input of the accepted deriver task T86e85c8c9a484988, so it is not
  edited here. Apply at C124's next edit. GR's TX49 one_line fix also remains
  deferred.
* REA_fact is named in S6c07b6b839fd452a. No prior-art search for it was done.

## Process lessons

* Review wording is a recurring source of defects. Note GR records two instances
  (C117's row (0,B,0); "never, at any M"). A third came from applying a review's
  correction list incompletely (V68797b01a59844f2 RC1, on submission
  Se0c7fe1a42cf475a). In this round it happened again in the survey revisions,
  including an unsafely rounded constant (V0909a684370d49e0 B11), and in this
  integration, where three C125 sentences taken from the accepting review's
  correction text were wrong (V7f9024351f81417b R2-R4). Reviews also err on their
  own, as the miscount recorded above shows.
* Cap wording-only revision loops. After two wording rounds that each introduce
  new defects, park the text and record the verified content as leads, rather
  than polishing indefinitely.

## Evidence pointers

* **Deriver:** S6c07b6b839fd452a; review Vf59b4285e64444e4.
* **Survey:** S09a0152534d74e2c, Sc7b5189d23a54129, S8834db9e6b7f4c90,
  Sc78326c4b7914d83; reviews as listed above.
* **Board:** topic:memory74 messages M3a47c1cec8894393 (seal) and
  M44427a5189834009 (reveal).
* **Integration:** task T53a3660437144474.
