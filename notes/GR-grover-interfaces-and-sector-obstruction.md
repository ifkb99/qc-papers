---
code: GR
date: 2026-09-22
title: "Grover reduces to predicate interfaces; subgroup sectors are an obstruction"
outcome: record
claims: [C122, C123, C124]
todo: [74]
---
# GR — Deciding the memory74 slate (2026-09-21/22)

## Why

The user asked for an assessment of topic:memory74 (note FG), then to "continue
with the research and see this bit through", trying promising ideas and
interleavings along the way. The assessment (arb message M717c24d46b0c413d and
earlier in the topic) argued that FG's ranking was by tractability rather than by
distance to the goal:

* TX48 and TX49 looked like one elementary fact about Grover's two-dimensional
  subspace.
* TX50 looked obstructable.
* TX51 is where the factoring bar actually sits.

## What was decided

* **C122** (coordinator-authored): uniform-start Grover output at any t reduces to
  a coin plus a predicate interface. TX48 is imported, for the coin only. TX49
  stays open and deprioritized. Four review rounds.
* **C123** (coordinator-authored experiment): Contract-B output of the
  two-orientation semiprime predicate is Simon's law for p xor q, and known
  branch-and-prune then factors heuristically on average. P4 failed and was left
  failing.
* **C124** (memory74_tx50_deriver): TX50 is an obstruction under its interface
  hypotheses.

  My brief stated a wrong premise, that the sector method "saves at most
  |<a> cap H| over a charged baseline that walks the quotient". The referee
  corrected it: the saving is relative to the explicit orbit, and there is none
  relative to the quotient walk.

  The deriver also corrected my sketch in two places. (d, |H|) gives a multiple
  of r, not r. And the reduction is polynomial only for polynomial d.
* **TX51** splits into two bars. Review V641e3fc372634f78 showed that my first
  bar (a rho order walk at sqrt(r), about N^(1/2)) is weaker than ECM and Pollard
  rho, which already factor with polynomial memory.

## Observation, not in a claim

This is review V5f7a0852665441d9's K5, re-derived here. Take Contract B over the
*restricted* factor predicates f_N AND [x extends w].
* If N_1(w)=0, the output is y=0 at every t.
* Otherwise, choosing t for M_w in {1,2} makes P(y!=0) a constant, and any y!=0
  certifies N_1(w)>0.

The C119 descent therefore factors N with no p xor q step. So C123's p xor q
question concerns only the unrestricted predicate.

## Process lessons

* Twice, a sentence taken from a review without re-derivation was wrong:
  - v2 used C117's row (0,B,0); it must be (0,1,0).
  - v3 said "masses never fix prefix counts at any M", which is false at M in
    {0,D}.

  SWARM's "re-derive; do not transcribe" names exactly this. Each fix now
  carries an exact detector with a mutant (revision_detectors_v2/v3/v4 in
  attempt A5df10f8d3a234159).
* Detector bugs of my own occurred in three runs:
  - a loose approximate identity;
  - a missing p,q>=2 bound;
  - the environment applied to the empty prefix.

  All logs are kept. The v3 run-1 source was overwritten, so from v4 on no
  source is overwritten.
* C123's first run executed without a board run record: run.start needs message
  IDs, errored, and was not gated. The identical seeded rerun was recorded. A
  design review would probably have caught C123's weak controls (C0, C2).

* The first integration submission was refused as stale. The upstream tasks had
  sealed todo/open/74 and transfers/TX48-TX51 as inputs, and the integration
  rewrites exactly those files. SWARM already says to keep backlog text in
  source_ref. The recovery took three provenance rounds (reopen, re-claim,
  byte-identical evidence plus a revalidation note, focused review), because the
  integration's own row text needed correcting each round:
  - TX48 wrongly said polynomial t needs no M;
  - TX48 said "equivalent to factoring" for Contract B;
  - TX50 scoped dominance and factoring wrongly;
  - TX50 omitted "known" subgroups;
  - TODO74 had an escape overstatement;
  - C122 dropped a "given H_elem" while condensing v4.

  For future tasks, list only claims and frozen contracts as inputs.
* Text detectors (out/memory74b/fix2-fix5) evolved under review. fix5 asserts
  qualifiers and t-conditions and logs the hash of every file it reads. It is a
  phrase sweep, and paraphrases can evade it.
* Deferred edits to sealed files, not applied so that the accepted dependencies
  stay fresh. Apply at the next edit of each file:
  - TODO74 (V915f0d1a57704024 N2): "only with polynomial d" should read
    "d <= poly(log N) ... suffices". C124 proves sufficiency, not necessity.
  - TX49 one_line (V68797b01a59844f2): add "given H_elem at binary t", and say
    that sampling *reduces to* prefix counting plus the coin. It is not "is":
    C122 gives only sampling <= counting.
  - TX50 (V5d5f994666f34906 O1, O2): "the family" in the factoring sentence has no
    antecedent (say "an N-family with these properties"); lines 68 and 72 break
    the line wrapping.
  - TX50 O3 and O4 are detector and scope remarks, needing no text change. O3:
    fix5 checks phrase presence plus four specific absent phrases; it cannot
    detect arbitrary unscoped statements or paraphrases. O4: G is not listed among
    the outside-the-hypotheses subgroups, and the list is not claimed exhaustive.

## Evidence pointers

* **Submissions:**
  - C122: Sc377b4ddf8aa4d5c, S0519269042584d57, Sb0f23b23c1e544d7,
    S76412a5b408a486a.
  - C123: S3ad6a279f9f2462c.
  - C124: S3c167ac49a4b4e0a.
* **Reviews:**
  - C122: V641e3fc372634f78, Vcae2d0a4cfb1407b, V7b25dbcb654749ac,
    V5f7a0852665441d9.
  - C123: Va27665ae402a41b6.
  - C124: Vbe013aaa5e46408d.
* **Provenance revalidations** (after the stale refusal, finding M17a80c3977074b49):
  - Tb3cbc3840a654e12: S8e1ae27e09824eec, S61145e55c01c4d2b, S316fde06fd7643d0,
    S13eaf10b46fe45fd (accepted); reviews V83636ade0dc549a1, V9393387257274f16,
    Vae659f0ae0524d0d, V915f0d1a57704024.
  - T1a9534f62e444d80: S1b2fa9711c7c4341, Sa49251a7220148eb, S203b13f214304f03
    (accepted); reviews Va76dac65b8b341bc, V2ba70051548746c8, V5d5f994666f34906.
* **Integration:** task T0b6b947664694fc3; reviews V68797b01a59844f2 and later.

Board and out/ evidence is local and gitignored. The claims carry the
interpretation.
