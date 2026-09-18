# HANDOFF.md as it stood before the 2026-09-18 migration (TODO 61)

Archived verbatim below this header: the checkpoint history, 2026-08-08 to 2026-09-17. Not
loaded into sessions. Live instructions were kept in HANDOFF.md or found in their owning
records; claims, notes, TODOs and the board own everything else here.

---

# Handoff — pick up here

Written for a session with no prior context. **Read `CLAUDE.md` first** for the
conventions and traps; this file is only *where things stand*.

**Checkpoint, 2026-09-17 (Claude coordinator): process changes and TODO 55 in flight.**
Pointers only; the facts live in the records named.
* **Process.** METHOD.md and SWARM.md changed on 2026-09-16/17: idea slates
  (phase 0), design review before compute (phase 2), exploration vs test,
  host-neutral role duties, gate-stop amendments. Summary and validation:
  SWARM_REVIEW_2026-09-12.md section 10; calibration addenda for the rounds in
  section 5. Codex has no role files: its briefs must cite SWARM.md "Role duties".
  A second 2026-09-17 pass (section 11) separates scientific authorship from the
  coordinator, adds SWARM.md "Revisions and corrections", and replaces the
  open-ended shadow-review rule with a packet-mode calibration trial and the
  condition that discharges it.
* **Research line.** Note SL records the slate and the output-contract premise.
  TODO 55 is done: C104 and note ES (exponent-slice decomposition proved; K2 measured
  against the dictionary engine; applied known mathematics; no memory crossover with dense
  exact streaming is claimed). TODO 57 (time at matched memory) is DONE: its
  graded statements are in C104's Limits and its three mechanism defects in note ES. TODO 56 (class-form output) is
  queued. C104 was integrated under review `V98eddcc1b4aa4a71` (accept;
  closure C326812cbdeb04e48), and its second edit carrying TODO 57 is under
  review on T392666b372104cd4. Its evidence is
  unversioned under the gitignored out/ tree, so nothing in git reproduces its numbers:
  TODO 58 holds that promotion, as TODO 59 does for the DD pilot.
* **DD-native PPS closed out.** Note DN records both accepted rounds
  (T563881b120674641, T163d291c2cdd4cd0), the byte result, and the refutation of C103's
  crossover as the mechanism for the node-count turn; TODO 13 keeps its original target
  (a real DD simulator's state diagram) open. The promotion was never integrated because
  its registered P5 passed on its two registered points and the hold was the
  coordinator's decision (M889892e48cad42b9); the module's promotion gap is TODO 59.
* **Coordinator-run, not independently reviewed.**
  experiments/experiment_dd_ordering.py and experiments/experiment_dd_alpha_replicate.py
  (untracked; rows and verdict logs in out/dd-ordering/, ignored) are the source of DN's
  refutation table. Unlike the two board rounds they had no referee, and DN says so where
  it reports them.

**Decision-diagram checkpoint, 2026-09-15: C102 / HD / TODO54 (task Ta294b7d88a4c4015).**
A single Claude coordinator session chose this direction at the user's
invitation; the user then authorized spending remaining usage with arb kept
current. C102 owns the proved adder counts (Walsh 3*2^(m-1)-2 against a
3m+4-node ROBDD) and the measured u_a scaling to n=8 (Walsh ~2^(3n+3)
against median ROBDD 9.7-14.3*4^n, with an N mod 4 class split). HD owns the pilots, three board runs (one native
crash, one preserved failed prediction P5a, one clean n=8 run), controls and
weak points. New code: `lab/bdd_count.py`, `lab/affine_pieces.py`, and two
experiments, which need `--with dd==0.6.0` and Python 3.12 for long runs.
No bound for u_a is proved and nothing escapes the C19/C48 order barrier.
TODO54 is done. TODO13 stays open and now cites C103. Nothing was committed:
C99-C103, their notes, experiments and lab modules are all uncommitted.

**Review status, 2026-09-15: both results reviewed, closed and integrated.**
* C102: task Ta294b7d88a4c4015. Round-1 referee V06b5088b79994f1c returned
  changes_requested after confirming the adder proofs line by line and
  recomputing every u_a row with q <= 20. Its catch: ROBDD size splits by
  N mod 4 and every n=7 fixture is N = 1 mod 8, so the all-N "step rose at
  n=8" was class mixing. All eight corrections were applied and frozen as
  Sa4101c802d464f5d. Round-2 referee V6840fa9b93014e09 verified the
  unchanged-science hash claim (snapshot digest 8cc4909846dfde6a, 16 files
  byte-identical, no run metadata rewritten), confirmed all eight
  corrections are present in the bytes, and returned changes_requested for
  three wording/scope defects plus one mechanical blocker. Those four are
  now applied here: the N = 1 (mod 4) monotone fall is stated as a property
  of medians that the two n=8 fixtures do not establish; the brute-force
  check is scoped to the reducer rather than the u_a counts; HD calls the
  archived script a reconstruction of the unarchived session check; and this
  block no longer describes its own review as pending. Round 3
  (Vee711b303f314378) accepted the v3 re-freeze Sb66f293fef1d4a37 and the task
  is CLOSED (closure Cf31ebb9eaea34122). Two record-level errors in the v3
  submission text are recorded in that closure rather than edited into the
  frozen submission: a wrong line range for the HANDOFF novelty advisories
  (they start at line 134, not >=524) and a wrong reason for dropping a
  retained failure log from the manifest. Neither touches a canonical file or
  any measurement.
* The mu(beta) width theorem was accepted by referee V2e9e6c322bee4511 and
  task T58ddf9e611124d65 is CLOSED (closure Cba59956c7eca4d8b). It was a
  proposal task, so it reaches the ledger through integration task
  Tf102392bd0124c1f, which creates claims/C103.md and
  notes/HB-beta-dichotomy-decision-diagram.md, copies the two experiments into
  experiments/ under their frontmatter names, and applies the two corrections
  the review requires (the "unloaded scratch" wording, and the alpha+mu
  crossover caveat). Tf102392bd0124c1f replaced Tbe30548d3e774494, which could
  not be claimed: the board revalidates dependency inputs recursively and went
  stale on prose-only drift in claims/C102.md (the accepted v3 corrections).
  The override and what it gives up are recorded in Mf2ab5cb32c584940.
* Do not re-run the science: runs R2fe53ff3836846c2 (crash, retained),
  Rdb87a87850f64ccc, Rb9caf0557f6e4b48, R5e2b4f1518e74ec6,
  Re9eda3cc191a448a and Rf7a4f1d224f54e5c are finished and archived.
* Process note from V6840fa9b93014e09: freeze HANDOFF.md last, or keep
  review-status prose out of the submitted candidate set. A status line that
  describes its own review invalidates its own freeze, which is what blocked
  acceptance of Sa4101c802d464f5d.
* What the mu(beta) result says, and what it does not: see claims/C103.md,
  which owns the statement, the proof and the limits. The one thing worth
  carrying here because it changes how the headline reads: the bound is the
  trivial 2^k for k <= alpha + mu(beta), so at LARGE beta MOST measured levels
  are trivial and the separation from PPS's base-2 growth is asymptotic in t.
  At small beta (beta=3, mu=1) the crossover is early and the separation is
  visible throughout, so the qualifiers matter. C103's Limits own the exact
  statement. TODO13 stays open: this is a canonical diagram of the pullback,
  not DDSIM's state diagram.

**Differential/boomerang checkpoint, 2026-09-14/15: C99/C100/C101 / DB/GB/BO / TODO51-53.**
A single Claude coordinator session, at the user's request:
* cancelled two superseded Claude board tasks, with reasons on the board;
* derived C99, the DDT count for X/Y-type Pauli pullbacks, and C100, its
  extension to any finite abelian group. There, PA4 in
  `experiment_group_bridge` is refuted and left failing on purpose.
* derived C101, the X–X OTOC as a boomerang connectivity count, with an
  exact carry butterfly and a compiled tail-bit symmetry (β ∈ {1, 3} when
  t−2 ≥ α).

The user then asked for referee review, so the coordinator opened write tasks
T79116bfece0e44a1 (C99) and Tf8bb4261531241e1 (C100) and spawned fresh
qsim-referee agents. Two coordinator process slips on C101 are recorded in
its closure: one submission was frozen despite a lint error, and one with a
stale manifest.
* **C99** was accepted and closed (V88c6aae3d06240d2, then
  V803a36546f284017). Its first attempt lease had expired during an API
  rate-limit pause and was reclaimed.
* **C100** was accepted and closed after two rounds of corrections
  (Vefabe53d66914d24, Vfb370f16c33e447a, then Ve40769e7918e4d0e).
* **Follow-up task T393699c97b9745ca** was accepted and closed
  (Vc00905ba7aba40ec). It applied or dispositioned C99's wording notes and
  relabelled the Z/8 maximum 55 and the row formula T = 56 − 8/o(d) as
  established by complete enumeration, with no conceptual proof.
* **C101** (OTOCs as boomerang counts, BO, TODO53) was accepted and closed
  after three rounds: V7cf842d53f784f58 and Va55e183ff79849a1 requested
  changes, then V16ef0281b61e45a1 accepted version 3.
  `experiment_boomerang_otoc` exits 1 on purpose, because its must-fail
  control C5 did not fail.
* **Wrap-up task T54da8c3405ec4c1c** applied or dispositioned the carried
  non-blocking notes and recorded the C101 referee catches. It was accepted
  and closed after one round of corrections (Vf2cc905b879f49df, then
  V55ef18e5bb7043a4).
* **Final micro-task Tcc271e7afe344499** fixes the last wording items. Once it is
  reviewed and closed, the board has no open assignments and C99–C101 have
  no carried notes.

On ToffoliModExp(7,3,3), the Pauli/mixed ratio for the clock-shift
observable is 0.64–1.88, in both directions. TODO50 is unchanged. Nothing
was committed.

**Disjoint-mask checkpoint, 2026-09-13: C98 / DW / TODO50.**
A Claude-coordinated round (Codex idle; board notice M3baabb2970914082)
proved C98's width-independent coefficient bound for disjoint masks and its
exact overlap failure range. DW owns the review history and the withdrawn
fixed-q width study, which was designed but deliberately not executed.
TODO50 now points away from that study. No code, experiment, manuscript,
commit or publication changed; the earlier Claude submissions
S57edfaa06f754400/Saacd7b3d73a24826 were untouched then (their tasks were cancelled 2026-09-14). Preserve the dirty tree.

**Sparse contraction and carry-catalog checkpoint, 2026-09-12: C96/C97 / SB/HC / TODO50.**
C96 owns the exact sparse-mask block contraction and charged costs; SB
owns its validated opt-in implementation, direct/frozen references and
fault controls. C97 owns the synchronized-prefix bound and stopped-prefix
catalog obstruction; HC records its proof/source review. The latter does
not bound actual reachable histories or the signed scalar. TODO50 alone
owns the new word-width/process-memory cost question and chronology-sensitive
closure discriminator. The broader simulation goal remains active and unresolved.

Core ran first in this continuing scientific session; the bounded sparse
candidate passed. Fresh Astra reviewers accepted the proof, code/evidence
and catalog mathematics. SB/HC and their integration task record final
documentation checks and review. Do not repeat completed sparse correctness
or fixed-width interval benchmarks. Existing C89/C94/packed sources, C95,
production helpers and manuscripts remain unchanged. Nothing was committed
or published; preserve the dirty tree. TODO34/TODO42 remain separate.

**Packed comparison and counting-closure checkpoint, 2026-09-12: C94/C95 / IS/SG / TODO50.**
The C94 allocation saving survives a stronger packed C89 comparator; IS
owns its paired measurements, controls and explicit process-RSS limit.
C95 owns the reachable convex-cell obstruction and exact all-shift set
recurrence, including the fixed-dimensional individual counting atoms and
missing repeated-closure bound. SG owns the pure-proof/source-audit evidence.
TODO50 now prioritizes that precise closure question. No faster full dirty
average or general quantum-simulation breakthrough is established.

Core ran first for the packed implementation; its native/control/benchmark
checks passed. Separate fresh Astra reviews accepted the packed experiment
and mathematical proof. IS/SG and the integration task record final doc
validation and review. Do not repeat either completed conditional benchmark
or start a rank/cell-count sweep for C95. No production helper, manuscript,
commit or publication changed. Preserve the dirty tree; TODO34/TODO42 remain
separate. The user's broader goal remains active and unresolved.

**Conditional streaming checkpoint, 2026-09-12: C94 / IS / TODO50.**
After the user's explicit handover from Claude, Codex implemented the
all-mask itinerary evaluator and completed blind native checks and a matched
allocation/runtime campaign. C94 owns the proof and cost/scope limits; IS
owns the measured improvement over current shared-lazy C89, the absence of
an established process-RSS saving, and all evidence/review provenance.
TODO50 remains open and owns the stronger packed comparator and full dirty
average question. The broader simulation goal remains unresolved.

Core ran first. The frozen proof, native reference, implementation and
scientific result received their required reviews. IS and the integration
task record documentation validation. Do not repeat the completed native
mask/boundary checks or current-B0 benchmark. Old Claude residual-rank work
was cancelled and preserved; older survey/reference submissions are not
support for C94. No production helper, manuscript, commit or publication
changed. Preserve the dirty worktree; TODO42 and TODO34 remain separate.

**Final-shift checkpoint, 2026-09-12: C93 / SC / TODO50.**
The board swarm derived an enabled contraction based on final signed shifts,
with matching forward and backward digit scans. C93 owns the proof, scope,
normalization and stronger charged competitors; the selected target now has
a guarded implementation. SC owns finite gate/cell evidence, source audits
and the preserved verifier-selection failure. Broader endpoint/schedule
extensions remain proof-only. TODO50 stays open and alone owns the next
compression/cost discriminator. The user's general simulation breakthrough
goal remains active; no practical memory advantage is established.

Core ran first; the bounded candidate and corrected independent cell checks
pass. SC and the coordinator write-task manifest record final documentation
validation and review. Other science suites were not rerun; TODO34 retains
its prior native crash. Do not repeat the completed mask/sector/cell pilots
or restart the older unimplemented disabled-target proposal. No production
helper or manuscript changed; nothing was committed or published. Preserve
the dirty worktree. TODO42 remains separate.

**Dirty-prefix checkpoint, 2026-09-12: C92 / DS / TODO50.**
The board swarm proved a common reversing involution and a polynomial
construction for the growing disabled Mersenne sector. C92 owns these proofs,
their exact scope and the stronger charged symbolic baselines; the polynomial
sector construction is unimplemented. Root discovered a nonzero full-space
three-macro scalar, independently checked by interval contraction. DS owns
the finite values, source/proof reviews and preserved latent verifier fix.
TODO50 remains open around growing enabled prefixes and alone owns the next
discriminator. No general memory breakthrough is established; the user goal
stays active. Do not repeat the completed mask search or disabled-sector pilot.

Core ran first; the bounded discovery and corrected independent verification
pass. Other eight science suites were not rerun; TODO34 retains its prior
native crash. Documentation validation and board completion are recorded in
DS. No production helper, manuscript or host setting changed; nothing was
committed or published. Preserve the dirty worktree. TODO42 stays separate.

**Shared-average checkpoint, 2026-09-12: C90 / C91 / SA / TODO50.**
The swarm proved and checked a complete two-macro shared-variable scalar
average, and established stronger clean-scratch rotation reductions. C91
owns those exact scopes; its clean scalar construction is unimplemented.
C90 owns the raw cubic-lift obstruction and affine-branch baseline. SA owns
evidence, source audits, the cancelled unexecuted recognizer pilot and the
preserved report/metadata corrections. TODO50 remains open around growing
dirty order-sensitive prefixes; it alone owns the next discriminator.
No general memory breakthrough is established, so the user goal stays active.

Core before/after, the corrected bounded experiment and the documentation
gate pass. Other science suites were not rerun; TODO34 retains its prior
native crash. All four reviewed tasks are closed and the redundant pilot
task is cancelled. No production helper or manuscript changed; nothing was
committed or published. Preserve the dirty worktree. TODO42 remains separate.

**Conditional controlled-add checkpoint, 2026-09-12: C89 / IX / TODO50.**
The board swarm proved and checked an exact interval representation for
actual modular-add prefixes after fixing outside scratch and controls.
C89 owns the proof, source baselines, inverse cancellation rule and outer-
average limitation; IX owns finite evidence and preserved review corrections.
TODO50 remains open and owns the next signed shared-variable discriminator,
including pointers to unaudited side proposals about global quadratic
structure. No general memory breakthrough is established; the user goal
remains active. Do not restart with another conditional interval size pilot.

Core before/after, the main experiment and the corrected independent verifier
pass. The documentation gate passes. Other science suites were not rerun;
TODO34's earlier native crash is unresolved. All four bounded board tasks
are accepted and closed. No production helper or manuscript was changed,
and nothing was committed or published. Preserve the dirty worktree;
TODO42 remains separate.

**Matching-reduction checkpoint, 2026-09-12: C88 / ME / TODO49.**
The board swarm proved that the carry grid survives exact matching
reductions and established a nonzero-pivot sparse recurrence. C88 owns
theorems and stronger known baselines; ME owns the independent matching
audit, charged allocation comparison, preserved import failure and source-
review corrections. The implementation improves our dense evaluator on
bounded fixtures; it establishes no best-known memory advantage or general
breakthrough. The user's goal remains active.

TODO49 is done at this scope; TODO50 alone owns the redirection to actual
shared-state controlled arithmetic. Mask-only signed variants stay within
known planar FKT, and another larger dimer benchmark is not the next step.
Core before/after and both bounded experiments pass. Other science suites
were not rerun; TODO34's earlier native crash is unresolved. All four board
tasks are accepted and closed. No manuscript was changed and nothing was
committed or published. Preserve the dirty worktree; TODO42 stays separate.

**Explicit dimer checkpoint, 2026-09-12: C87 / DM / TODO48.**
The board swarm supplied an explicit growing planar carry family and a
stronger direct-dimer specialization. C87 owns the permutation/minor proof,
nonzero scalar formula, complement/autocorrelation identities and local
rank obstruction. DM owns exact evidence, source audits and retained
independent-verifier defects. TODO48 is done; TODO49 owns further matching
reductions and the conditional cost or signed-query discriminator. Known
methods explain the current result. No practical memory advantage or
breakthrough is established, so the user goal remains active.

Core, the main experiment, geometry checks and final corrected independent
reference pass. Earlier worker failures remain documented; passing symmetric
fixtures hid a missing incoming carry until source review. Other science
suites were not rerun, and TODO34 retains the prior native crash. All four
bounded board tasks are accepted and closed. No manuscript was changed and
nothing was committed or published. Preserve the dirty worktree; TODO42
remains separate. Do not start with a larger known-FKT width pilot.

**Planar carry checkpoint, 2026-09-12: C86 / MG / TODO47.**
The board swarm completed the bounded interacting-addition investigation.
C86 owns the exact planar matchgate reduction, implementation contract,
resource limits, geometry deduction and dirty-scratch arithmetic baselines.
MG owns independent evidence, primary-source scope and board review.
TODO47 is done; TODO48 owns the concrete growing circuit/query discriminator.
This is an arithmetic application of known machinery. No practical memory
advantage or general simulation breakthrough is established, so the user
goal remains active. Do not begin with a large Pfaffian width pilot.

Core and both bounded exact experiments pass. Other science suites were not
rerun; TODO34 retains the previous native crash. All four bounded board tasks
are accepted and closed. No manuscript was changed and nothing was committed
or published. Preserve the dirty worktree; TODO42 remains a separate question.

**Carry-prefix checkpoint, 2026-09-12: C85 / CP / TODO46.**
The board swarm completed the bounded actual-adder investigation. C85 owns
the logical-prefix invariant, exact carry-correlation implementation, resource
limits and corrected conditioned escape witness. CP owns primary-source
audits, independent evidence and retained failures. TODO46 is done; TODO47
owns the next interacting-carry discriminator. Known carry-correlation
methods explain the isolated-adder result, so no breakthrough is established
and the user goal remains active. Do not repeat an isolated-adder width pilot.

Core and the main exact experiment pass. The independent witness experiment
intentionally exits one for the refuted no-cut escape prediction; CP explains
its separately passing checks and corrected reporting. Other science suites
were not rerun, and TODO34 retains the earlier legacy native crash. All four
bounded board tasks are reviewed and closed, with no experiment processes
remaining. No manuscripts were changed. Preserve the uncommitted worktree;
TODO42 remains a separate sampling-cost question.

**Generic stabilizer baseline checkpoint, 2026-09-12: C84 / PR / TODO45.**
The board swarm closed the local recognition discriminator. C84 owns the
projector-reflection theorem, exact Clifford replacements and equivalence to
C83; PR owns source audits, bounded independent evidence, failed-run provenance
and cost/numerical limits. TODO45 is done; TODO46 owns the next question beyond
individual-cell closure. No breakthrough is established; the broader user
goal stays active. Do not repeat this local rank comparison or QC's pilot.

Core and both bounded experiments pass; the prior legacy native-crash result
remains with TODO34 and other science suites were not rerun. All four bounded
board tasks are reviewed and closed, with no experiment processes remaining.
No production library or manuscript was changed in this turn. Preserve the
uncommitted worktree. TODO42 remains a separate sampling-cost question.

**Quadratic-observable checkpoint, 2026-09-12: C83 / QC / TODO44.**
The board swarm completed the bounded exact-query investigation. C83 owns
local nonlinear closure, streamed-cell implementation and known-theory limits;
QC owns the stronger-baseline loss, independent audits, failed-verifier
provenance and validation. TODO44 is done at this scope; TODO45 owns the next
comparison with phase-sensitive stabilizer splitting/recombination. TODO42
remains open separately. No breakthrough is established; the user goal remains
active. Do not repeat the recognized-family width pilot as the next task.

Core and both bounded experiments pass. The previous incomplete legacy suite
remains with TODO34; other science suites were not rerun. All four bounded
board tasks are reviewed and closed, and no experiment processes remain.
Changes are uncommitted; preserve the dirty worktree. Use the board resume
workflow below, then the individual claim and current todo for details.

**Direct CNOT-memory checkpoint, 2026-09-12: C82 / CF / TODO43.**
The user requested a board-coordinated research swarm. Main integrated an
opt-in binary Walsh-coordinate frame into the existing permutation engine;
`final_terms` can expose physical keys lazily. C82 owns the proof, interface
and prior-art limits. CF owns independent tiny-law/core evidence, the complete
allocation comparison and its stronger-baseline loss. This is bounded progress,
not an established breakthrough. TODO44 owns the next structural-observable
question; TODO42 remains open on sampler selection.

All four bounded board assignments are reviewed and closed; no experiment
processes remain. Core and bounded independent experiments pass. The full
legacy permutation suite terminated with a native crash and remains incomplete;
CF preserves its status/log and TODO34 owns reliability. Other science suites
were not rerun. Changes remain uncommitted; preserve the dirty worktree.

**Coordination board installed, 2026-09-12.** The reusable implementation lives
in `/home/djneko/Workspace/agent-research-board`; the installed `arb` command reads
this project's `arb.toml` and keeps operational data in `out/agent-board/`.
Start with `arb --actor coordinator resume`, read the inbox, then acknowledge its
receipt. The updated `qsim-research` skill supplies the workflow. `arb board`
and `arb thread topic:general` let the user inspect it; `arb post` contributes.
One coordinator remains the canonical writer; worker output ownership is advisory,
with no host-enforced isolation. This tooling task did not resume experiments or
change the research frontier. Continue to the scientific checkpoint below.

**Resumed at user request, 2026-09-12: TODO41 completed at bounded float scope.**
`NestedPhaseProgressions` now implements C81 in the existing work-first loop.
C81 owns its input contract/resource bounds; NS owns authoritative tiny-law,
actual RNG and matched returned-sample evidence, including failed verifier
drafts and the retained exhaustion/precision controls. TODO42 owns the next
question: selecting for sampling cost rather than construction alone.
Do not repeat TODO40's pilots or TODO41's integration checks as the next task.
Core/lab/claims pass; documentation validation is recorded in NS. The other
six science suites were not rerun. No experiment processes or active delegated
tasks remain to await. All changes are uncommitted; preserve the dirty worktree.
Keep research scratch files and backups inside `research/out/`.

The earlier 2026-09-11 usage-limit pause ended with this explicit resumption.
This follow-up used the original qsim-research workflow, with main owning
implementation/integration and lower-cost initial testers as specified by TODO41.

**Updated 2026-09-11.** The user requested an independent review, delegated
corrections, and a new research direction. The lower-cost correction agent
made partial edits before hitting its usage limit; the main agent finished
and audited them. Changes are uncommitted and preserve the user's pre-existing
edits to both papers and both abstract workshops. All nine science suites
(including CUDA) passed at the consolidation checkpoint; later scoped
validation is recorded with each follow-up below.

---

## Where things stand

**Both papers are drafted in full and have had a consistency pass and a
structural edit.** `PAPER_A.md` is the exact cost model, `PAPER_B.md` the
2-adic result for modular exponentiation. `ABSTRACT.md` and
`ABSTRACT_SHOR_2ADIC.md` are abstract workshops and prior-art dossiers, **not**
superseded copies to edit in parallel.

**The user-requested two-file consolidation is now in those existing papers.**
Paper A integrates clean-code equivalence and tensor-rank limitations in
§6.4–6.5. Paper B integrates balanced reduced tails and the conditional
baselines in §10.3 and §11.1–11.4. Their abstract/conclusion/scope passages and
Appendix A maps were updated together. These are working consolidations, not
a claim that the manuscript titles or publication positioning are finalized.

Review findings and their counterexamples: `REVIEW_2026-09-09.md` and
`experiments/experiment_independent_review.py`. The original review is a dated
record; current statuses live in C17/C44 (peak relation), F12 (scalar periods),
C21 (onset), C42 (useful fractions), and the revised manuscript scope.
Paper A now explicitly credits the diagonal Pauli/Fourier identity as prior
art. The supplied Fourier and Toffoli circuits are not a same-permutation
comparison on the full space.

The new result is in **C45/C46 and `notes/CT-finished-control-contraction.md`**:
finished-control contraction for the reduced pre-QFT observable. The existing
`perm_pps.py` now has an opt-in `trace_plus` API, with regression tests and a
reproducible experiment. TODO 12c is closed. Default full-operator propagation
is unchanged. Keep its cost distinct from reduced-observable cost when editing
Paper B; do not present the latter as Shor output sampling.

`todo/INDEX.md` is the ranked list. **TODO 14 / TODO 32,34,42** now point from
C73/WE's implemented word envelope and C74–C75's physical-coordinate work to
runtime/host reliability, the wider comparison and additive physical phases.
TODO 33 and TODO 35 are complete at bounded scope; C75/CG own the clean-gate
evidence and C76/GS the additive-phase boundary. C53/UT complete the stronger
exponent-unitary diagnostic. TODO 36 is complete at bounded scope: C77/UP
own the support-based uniform-prefix certificate and independent checks.
TODO 37 / FB complete the Fourier-feedback-aware comparison at bounded scope.
C78/WF own the work-first sparse-row/progression proof, opt-in implementation
and bounded same-output comparison; TODO 38 is complete at float diagnostic
scope. C79/ER own complementary earlier-phase progression covers and their
complete-output comparison and opt-in sampler. C80/CW complete TODO39 with
an opt-in root-mass proposal and independent same-law/RNG/cost checks.
TODO40 is complete at mathematical / bounded amplitude scope. C81/NC own the
multiple-phase cover, exact truncated construction selector and independently
reproduced long-row comparison. TODO41/NS now complete production sampler
integration and a bounded matched returned-sample comparison; TODO42 owns the
next sampling-cost selector discriminator.
C60/§PE own the proved
oracle-level error budget; C61/§VP now supply the exact-input oracle/kernel.
TODO 24 remains OPEN for deferred working-precision complexity. C63/§RA
supply certified rejection; C64/C65/§FW supply verified finite-work/scalar
proposals and clarify the binary benchmark limitation. C59/§CH
own the coherent-history samplers and regrouping comparison. General exact
sampling and numerical certification remain open; C57 separately distinguishes
complete-output evaluation from omission bounds.

**The requested sequence 15a → 15b → 15c → 15d is complete at first-pass
scope**, with each item moved to `todo/done/`:

| Investigation | Claim | Working record |
|---|---|---|
| Balanced controls and certified tail residual | C47 | `notes/BC-balanced-control-channels.md` |
| Tensor ranks and running-residue realization | C48 | `notes/TR-walsh-tensor-memory.md` |
| Conditional actual order-finding distributions | C49 | `notes/CF-conditional-order-finding.md` |
| Same-layout scratch-extension equivalence | C50 | `notes/SE-scratch-extension-equivalence.md` |

The harness now writes structured reports, including failures. New reusable
helpers are `lab/tensors.py` and `lab/semiclassical.py`; they are covered by the
existing science suites. Experiments are not silently promoted into manuscript
theorems: the claims and notes state their scope, and TODO 14 remains the
broader compressed inverse-QFT frontier. In particular, the conditional sampler
relies on commutation on the reachable clean arithmetic subspace, not on the
arbitrary full-scratch branch maps.

**Further research completed after consolidation: TODO 16 / C51–C52 / §RO.**
`lab/reachable.py` uses the existing `walsh.py` classical replay kernel on
selected basis inputs, avoiding full-scratch tables. `eigenphase_path` in
`lab/semiclassical.py` supplies the explicitly order-informed spectral baseline.
The experiment compares compiled sparse, high-level modular, orbit-vector and
latent-eigenphase sampling, charging discovery/setup. The result is a useful
negative: the compiled route removes scratch overhead but does not beat these
classical baselines. C52 / Cleve et al. §6 explain why known order allows
sampling without a work vector despite a larger coherent Schmidt rank.

Reproduce with Python 3.12 as in §RO: an expanded 3.14 experiment exited 139
without traceback, whereas the complete 3.12 run passed. The partial failed
log is preserved; its crash cause was not established. Fresh science-suite
logs are `out/consolidation_test_<suite>.log`; all nine pass, with GPU checks
enabled where supported. Python 3.12 was used for the long PPS/claims suites.
No commits have been made.

**Research alignment skill installed and smoke-tested.** `AGENTS.md` now
routes Codex to the existing repo instructions and the user-level
`qsim-research` skill. See §SK for source/discovery locations, the frozen-skill
hash, three independent behavioral trials and their limitations. The skill
does not carry duplicate scientific facts or prescribe the next frontier.
`CLAUDE.md` now makes clear that conceptual discussion and document/skill-only
edits do not automatically launch science suites. These instruction changes
did not alter the scientific code or either manuscript.

**Next experiment completed: TODO 17 / C53 / §SD.** The single-defect
experiment preserves the original ascending-power gate order and identifies
which initial eigenphase coherences affect measured outputs. The old physical
work-mixture model fails, but an input-specific exponent-phase replacement
restores scalar sampling. This is a useful negative for a proposed memory
barrier, not a simulation breakthrough. `lab/spectral.py` adds capped dense
diagnostics; `eigenphase_path` accepts optional product exponent phases while
preserving its default ideal behavior. The existing circuit engines are
unchanged. See §SD for controls, precision/resource limits and validation.
Neither paper nor either abstract workshop was edited during this experiment.

**Active user goal: continue research with lower-cost initial testers.** Two
`gpt-5.6-luna` agents tested the physical work mixer, final-eigenphase prefix
factorization, sparse rejection, scaling and edge cases. Main implemented
`lab/prefix.py`, reviewed and reran their work, and added regressions. C54 / §PF
own the result and limitations; TODO 18 is closed. The mechanism is verified
within its explicit input/precision contract, not established as a breakthrough.
The thread goal remains active; do not mark open-ended research complete merely
because this experiment passed. Its next discriminating test was TODO 19,
now completed below.
The papers and abstract workshops remain untouched in this follow-up. Read §PF
for a caught reference-allocation error, its missing original log, and the
corrected pre-allocation guards; do not erase that audit history.

**Localized-defect follow-up completed: TODO 19 / C55 / §LF.** The sampler
now avoids both orbit and prefix arrays for one finite-basis-support defect.
The proof, known-order/index contract, averaged (not pointwise) rejection cost,
and simple ideal-approximation bound live in C55. Main recovered partial
lower-cost work after usage limits; an initial progression run and a resumed
independent audit provide agent verification. The unfinished formula verifier
had a shared normalization error caught by independent references; §LF
preserves that failure and the corrected evidence. Core/lab/claims pass;
the other six science suites were not rerun for these helper additions.
The impact experiment deliberately retains a failed monotone-dilution intuition;
read its header before treating its nonzero exit as a regression. TODO 20 was
the next test, completed below. No manuscript/abstract edits or commits in this follow-up.

**Periodic follow-up completed: TODO 20 / C56 / §PS.** Three lower-cost agents
tested the harmonic route, cancellation thresholds and conserved coarse sectors;
main implemented/audited a stronger multi-defect sampler and reran their work.
`lab/periodic.py` and the new `sequential_path` in `lab/semiclassical.py` retain
the original gate order. C56 owns the proof, input/precision contract and known
symmetry/MPS/measurement-calculus positioning. §PS records physical Circuit
validation, resources and failed controls, including an initially false odd-block
visibility assumption. Core/lab/claims pass; the other six science suites were
not rerun for these additions. TODO 21 was the next question, investigated below.
The user goal remains active; this experiment is not a general simulation
breakthrough or completion of the open-ended goal. No manuscript/abstract edits
or commits were made in this follow-up.

**Localized symmetry-breaking probe completed: TODO 21 / C57 / §LC.** Three
lower-cost agents tested the physical kick, actual-state error bounds and
co-moving support cone; one also verified the low-rank boundary identity.
Main implemented exact integer cover counting in `lab/periodic.py`, audited
and reran all experiments, and added lab/claims regressions. C57 owns the
proof, the distinction between approximation and exact sampling, and the
remaining sector sum in complete-output evaluation. §LC preserves a failed
bare-count control and a corrected false streaming-memory report. Core,
lab/claims and the documentation gate pass; the other six science suites
were not rerun for this helper addition. TODO 22 owns the next diagnostic.
The broader goal remains active. No manuscript/abstract edits or commits.

**Physical-phase probe completed: TODO 22 / C58 / §SR.** Three lower-cost
agents tested physical circuits, exact orbit phases and routed-sector formulas;
main implemented `RoutedOrbitCircuit`, audited and reran all three experiments,
and added regressions. C58 owns the proof and the distinction from broader
normalizer-circuit claims. §SR retains an input-specific vacuous control,
an incorrect physical work-qubit embedding and resource/precision corrections.
Core/lab/claims pass; the other six science suites were not rerun for this
helper addition. TODO 23 completed its next discriminating test below. The open-ended
user goal remains active. No manuscript/abstract edits or commits.

**Coherent-route probe completed: TODO 23 / C59 / §CH.** Lower-cost agents
tested formulas, exact regrouping, compiled indexed gates and the full tiny
sampling transition law. Main implemented `lab/coherent_routes.py`, audited
and reran all four experiments, and added lab/claims regressions. A primary
literature connection supplies a no-rejection alternative through prefix
amplitudes; this applies an established method, not a new sampling principle.
§CH records its tradeoff with rejection and the failed verifier controls.
Core/lab/claims and the documentation gate pass; the other six science suites
were not rerun for these helper additions. TODO 24 owns the next
numerical-certification question.
The broader goal remains active. No manuscript/abstract edits or commits.

**Numerical follow-up progressed: C60 / §PE.** Current next actions are ONLY in
TODO 24. Three lower-cost agents supplied initial tests; main audited/completed
them and added `lab/sampling_error.py`. This exact-rational budget helper is
not a certificate for the existing float sampler. Core/lab/claims pass; the
other six science suites were not rerun. The production sampler, papers and
abstract workshops are unchanged. The broader research goal remains active.

**Verified exact-input implementation progressed: C61 / §VP.** Main implemented
`lab/verified_prefix.py` using the shared finite contraction; lower-cost agents
tested the backend, independent reference and integer kernel. Main corrected
overstated verifier predicates, reran tests, and recorded their failed reports.
The existing float sampler preserves its API/mathematical algorithm but the
shared arithmetic loop can change last-bit rounding; it is still uncertified.
Core/lab/claims pass (verified tests run with the optional backend); other six
science suites were not rerun. Current next actions are ONLY in TODO 24.
No manuscript/abstract edits or commits. The open-ended user goal is active.

**Enclosure-geometry follow-up progressed: C62 / §NG.** Lower-cost initial
probes and main's verified implementation separate interval wrapping from
adaptive doubling overhead. Norm-aware error transport is opt-in: the matched
timing comparison does not justify replacing rectangular mode. Main corrected
an overstatement that the previous doubling endpoint was a required minimum.
Core/lab/claims and affected full-law regressions pass; the other six science
suites were not rerun. TODO 24 owns the certified rejection/primitive-constant
boundary. No paper/abstract edits or commits; the broader goal stays active.

**Certified rejection follow-up progressed: C63 / §RA.** Main implemented
verified proposal/history selection and one-sided absolute-width acceptance,
reusing the existing finite contraction. Lower-cost agents supplied initial
exact-edge/full-law tests and a proof/code audit; main strengthened vacuous
verifier predicates and retained every failure. The matched timing comparison
is now between two methods meeting the same requested accuracy, not against
uncertified float rejection. The stronger linear-depth finite-work proposal
was outside that certificate and is now covered by the follow-up below.
Core/lab/claims and affected full-law regressions pass; the other six science
suites were not rerun. No manuscript/abstract edits or commits. The open-ended
user goal remains active; this is a bounded result, not a general breakthrough.

**Verified finite-work follow-up progressed: C64/C65 / §FW.** Lower-cost
agents supplied exact-tree/full-law tests and a proof audit; main implemented
and validated both opt-in proposal methods. The skill's stronger-baseline
check exposed a collective late-work invariance, so the binary fixture is
not evidence of persistent observable work mixing. §FW records the original
and observable-initial-rotation tests, the matched timing scope, replay
control and preserved verifier corrections. The initial odd-block probe
provides a more discriminating fixture; current next work lives ONLY in
TODO 25, with the backend-wide bit-cost issue deferred in TODO 24.
Core/lab/claims and affected full-law regressions pass; lab also passes with
the optional backend absent and verified checks explicitly skipped. The other
six science suites were not rerun. Documentation validation is recorded in
§FW. No manuscript/abstract edits or commits; the broader goal stays active.

**Odd-block certification completed: TODO 25 / C61–C64 / §OB.** The exact-input
API now optionally supports three-state work blocks while preserving the
binary default and rejecting its scalar shortcut on the new input. Main
implemented dimension-aware arithmetic and planning; lower-cost agents supplied
full-law tests, edge probes and a read-only proof audit. Main repaired a
same-node replay test and preserved the weaker/failed reports. Complete-law,
edge and same-accuracy comparisons pass. Core/lab/claims, binary full-law
regressions and optional-backend-absent lab checks pass; the other six science
suites were not rerun. See §OB for resource/timing limits and documentation
validation. An initial algebra/literature check motivates TODO 26; it is not
yet a memory-compressed sampler. No manuscript/abstract edits or commits.
The open-ended research goal remains active.

**Forward-memory experiment completed: TODO 26 / C66 / §FM.** Main proved
and tested the prefix-independent state-approximation contract, keeping exact
backward effects and separately charging numerical mass errors. Lower-cost
agents supplied a proof/negative-control audit and initial mixing probes;
main repaired verifier predicates and reran the results. Full tiny laws,
the maximally mixed yet correlated control and binary null pass. No production
sampler stored fewer checkpoints at that stage. TODO 27 subsequently completed
the uniform-phase and baseline question below; §FM records the earlier scoped
validation. Core/lab/claims pass; the other six suites were not rerun. No
manuscript/abstract edits or commits; the broader research goal stays active.

---

## Recent bounded follow-ups

**TODO 27 completed: C66–C67 / §UG.** A certified phase mesh makes the
fixed-mixer gap quantitative and uniform; the sufficient approximation
warm-up remains impractical. Main implemented opt-in exact checkpointing,
with lower-cost initial tests and independent proof audits. Main repaired
verifier predicates, precision replay and matrix-lifetime accounting, then
reran full-law and wide-path comparisons. Core/lab/claims, binary/odd-block
full-law regressions and backend-absent lab checks pass. The other six
science suites were not rerun. §UG owns evidence, failures and documentation
validation. Default storage and both manuscripts/abstract workshops remain
unchanged; no commit was made.

**Stronger lead initially validated: C68 / §RI.** Main derived reverse
instrument sampling with boundary rejection; lower-cost agents checked the
density/pure-trajectory proofs and tiny raw-effect laws. Main fixed further
interval predicates, reran the probe and added a rational nonnormal-matrix
claim regression. Claims pass after that addition. This is an exact
mathematical construction plus tiny diagnostic; a certified reverse production
sampler did not exist at that checkpoint. TODO 28 subsequently completed its
implementation/error test below. No general breakthrough or novelty is established. The
open-ended user goal remains active.

**Reverse implementation completed: TODO 28 / C69 / §RV.** Main reused the
existing mass/acceptance lemmas to certify an unnormalized pure-vector
trajectory. Lower-cost agents audited the proof/code and supplied full-law
and edge/comparison probes; main strengthened predicates and reran them.
Working-state storage shrinks, but the matched timing comparison favors full
storage and checkpointing. Default samplers and manuscripts remain unchanged.
Core/lab/claims and backend-absent lab checks pass; the other six science
suites were not rerun. §RV records evidence, corrected verifier failures and
documentation validation. At this checkpoint TODO 29's representation was
not yet implemented; its subsequent completion is recorded below. No commit
or general-breakthrough claim was made.

**Local-error experiment completed: TODO 29 / C70 / §QD.** Exact integer
projective compression and a joint output/work error argument remove the
full-tree requested-coordinate precision factor. Lower-cost agents audited
the algebra/code and supplied probes; main corrected sector conditioning,
per-boundary law and bit-cost predicates and reran the full-law and serial
comparison experiments. Lower precision did not give faster returned samples.
Core/lab/claims, backend-absent lab and documentation gates pass; the other
six science suites were not rerun. §QD owns reports and preserved failures.

**Coherent merging completed: TODO 30 / C71 / §CM.** Main implemented a
sparse backward vector adapter; lower-cost agents independently tested route
sets, amplitudes, complete tiny laws and instrumented costs. Main strengthened
predicates and corrected history accounting, endpoint work and a sector-modulus
mix-up before rerunning. The skill's stronger-baseline check exposed a smooth
neighboring-route approximation; the separately frozen separated-label test
retains the implementation advantage. C71 owns the proof and limitations,
§CM the evidence and preserved failures. This remains FLOAT-only, not C70's
certificate transferred to a new state. Core/lab/claims and backend-absent
lab pass; other six science suites were not rerun. Documentation validation
is recorded in §CM. TODO 31 subsequently completed the merged-prefix oracle
investigation below. Manuscripts/defaults are unchanged; no commit or general-
breakthrough claim. The user's open-ended research goal remains active.

**Merged prefixes completed: TODO 31 / C72 / §MP.** Main implemented a
batched sparse prefix oracle feeding the unchanged gate-by-gate sampler.
Lower-cost agents supplied initial complex-prefix, full-law and comparison
tests; main strengthened controls/cost attribution and reran them. Removing
rejection trades increased working storage for less retry variability, without
a clear timing win over C71. The failed output-phase control became a correct
positive invariance check; internal phases still matter. Core/lab/claims and
backend-absent lab pass; the other six suites were not rerun. §MP owns raw
reports, the unresolved transient comparison error and documentation validation.
TODO 32 alone owns the broader fixed-alphabet/larger-k question. No manuscript,
abstract, default or commit changes; the open-ended goal remains active.

**Fixed alphabets progressed: TODO 32 / C73 / §FA.** Fixed alphabet size,
not just two labels, permits polynomial reached support in the arithmetic
model. A new input-only class separates larger structural inputs from capped
legacy history methods. Lower-cost agents supplied exact-set and full-law
tests; main strengthened guards, cost predicates and reference comparisons.
The genuine nine-insertion complex-prefix and complete-law checks pass.
Core/lab/claims and backend-absent lab pass; the other six suites were not
rerun. FA records one unresolved exit-139 run followed by two passing repeats,
including crash tracing. A tighter global word envelope is proved and checked
only at set level; TODO 32 retains its production integration and wider
matched comparison. Defaults/manuscripts/abstracts are unchanged; no commit
or general-breakthrough claim. The open-ended goal remains active.

---

**Word envelope implemented, wider run unfinished: C73/WE / TODO 32.**
Main integrated opt-in word support; lower-cost initial testers and main's
stronger reruns verify exact sets, tiny complete laws, complex prefixes and
actual returned-sample accounting. Two main serial timing attempts aborted
with different native faults; a per-row incomplete checkpoint now preserves
completed rows. Nearby standalone and system-Python preflights pass, not
establishing stability. WE owns reports and scoped gates. TODO 34 owns the
read-only host findings, including an old-microcode lead, and the diagnostic
plan. No production arithmetic was altered to mask a crash; no firmware or
system settings changed. Full timing conclusions remain deferred.

**Physical-coordinate algebra progressed: C74/PC / TODO 33.** Main derived
the coordinate/phase-oracle distinction; lower-cost independent checks and
main's exact reruns validate the elementary reductions and controls. This
clarifies the supplied-input promise rather than solving discrete logarithms.
The clean-ancilla physical implementation was still open at this checkpoint;
the subsequent C75/CG follow-up is recorded below. The new claim regression
passes; PC and WE own validation scope. Neither manuscript nor
abstract workshop changed; no commit. The open-ended user goal remains active.

**Clean physical gates completed: TODO 33 / C75 / CG.** Main derived a total
reversible coordinate construction and clean reflection rotation. Lower-cost
agents audited its proof, phase conventions and initial fixtures; main fixed
verifier defects, simplified exact gate expansions to stay under the original
cap, and reran complex-gate and complete-output checks. The stronger static
grouping still wins on structural simplicity: no speedup/breakthrough claim.
CG owns reports, two preserved construction-cap failures and validation.
Core/claims pass; the other seven science suites were not rerun, and shared
production helpers are unchanged. No native crash occurred in these small
checks; this does not resolve TODO 34. No firmware, manuscripts, abstracts,
defaults or commits changed. TODO 35 alone owns the next Gauss-sum phase
experiment and the scope of its lower-cost pilot/prior-art check. The
open-ended user goal remains active.

**Additive phase completed: TODO 35 / C76 / GS.** Main derived an exact
prime-field full-sector-coupling proof and state-specific Gauss bound;
lower-cost agents supplied initial physical/output tests and proof audits.
Main strengthened predicates and reran the physical, exact-polynomial and
complete-output checks. GS owns evidence, failures and resource limits.
This is a boundary of the current representation, not sampling hardness or
a new efficient algorithm. Core/claims and documentation checks pass; other
seven science suites were not rerun. No production helper, manuscript,
abstract, default, host-setting or firmware changes or commits were made.

**State-aware follow-up progressed: TODO 36 / C53 / UT.** A general
early-exponent-unitary criterion is stronger than scalar phase transfer,
but the joint-state optimum loses to omission on the actual measured law.
Lower-cost initial tests caught a synthetic normalization error; main repaired
the verifier and added a physical-suffix comparison, independently audited by
a closed Fourier formula. A parity/no-signalling argument then supplies an
observable obstruction for the entire early-register-only class. UT owns
the numerical evidence, preserved failure and validation; these are not
certified decimal bounds or general hardness. TODO 36 alone owns the next
finite-displacement/uniform-prefix experiment, subsequently completed below. The broader goal stays active;
TODO 34's host issue and TODO 32's timing comparison remain unresolved.

**Uniform-prefix follow-up completed: TODO 36 / C77 / UP.** Finite orbit
displacement certifies a reduced output marginal even with dense sector
coupling. A stateless exact-integer helper checks the sufficient condition
without history enumeration. Lower-cost initial count and two independent
output routes were audited and strengthened by main; UP owns the tests,
resource limits and retained verifier failures. Uniform prefix bits remain
correlated with their suffix. TODO 37 owns the feedback-aware comparison.
Core/lab/claims, backend-absent lab and documentation checks pass; the other
six science suites were not rerun. No existing sampler, manuscript, abstract,
default, host setting, firmware or commit changed. This is not a full sampler
or a general breakthrough, and the broader research goal remains active.

**Feedback and work-first follow-up: TODO 37 / C77–C78 / FB–WF.** The exact
conditional-boundary law is independently validated. Retaining Fourier
feedback improves the approximation but does not meet the frozen accuracy.
A different latent measurement then gives coherent sparse-row sampling and,
for a restricted late-phase schedule, a few-progression representation. C78
owns the proof/input contract; FB/WF own numerical evidence and preserved
verifier failures. Lower-cost agents ran initial tests and independent audits;
main strengthened resource/control predicates and reproduced the final laws.
The progression sampler and end-to-end comparison remain TODO 38, not a
completed production feature. Claims and documentation gates pass at this
checkpoint; earlier C77 core/lab checks remain the scoped production gate.
No manuscripts, abstracts, defaults, host settings, firmware or commits changed.
The open-ended research goal remains active.

**Work-first sampler completed: TODO 38 / C78 / WF.** Main implemented the
opt-in helper using existing progression sampling. Lower-cost initial testers
supplied joint-law, edge and same-output comparisons; main repaired weak
controls, counter/allocation accounting and a leaked-width predicate, then
reproduced them on system Python. Main additionally enumerated the actual
tiny RNG decision table, including gcd lifts, and tested an independent FFT.
WF owns results, failed verifiers and validation. Numerical certification,
arbitrary multi-defect schedules and timing superiority are not established.
The next earlier-phase cycle observation is frozen in TODO 39; its preliminary
proof audit caught a cancellation-revival caveat before implementation.
No manuscripts, abstracts, defaults, host settings, firmware or commits changed.
The user's open-ended research goal remains active.

**Earlier-phase coefficient pilot completed: C79 / ER / TODO 39.** Main
derived the binary-residue cycle; a lower-cost independent audit checked the
proof, and another agent tested literal columns/refined rows. Main caught the
old-cancellation revival caveat, repaired reference/work-accounting predicates
and reproduced the bounded results. ER owns numbers, explicit witnesses and
preserved failures, including two wrong column routes that agreed with each
other but failed the independent row check. Claim regressions pass. The
existing C78 production sampler does NOT support the new schedule; TODO 39
retains the actual-output comparison before any integration or broader claim.
No manuscripts, abstracts, defaults, host settings, firmware or commits changed.
The open-ended user goal remains active.

**Earlier-phase output checkpoint completed: C79 / ER / TODO 39.** Three
lower-cost routes compared geometric Fourier laws, a feedback-aware dephased
candidate and independently assembled Circuit/statevec references. Main
strengthened normalization, controls, cumulative work/allocation accounting
and raw-report retention, then reproduced and cross-compared the full laws.
The strongest tested approximation leaves detectable interference; uniform
first-bit marginals do not eliminate the suffix coherences. ER owns the
accuracy classifications, costs and failed verifiers.

Main also derived a complementary remaining-history cover and a two-history
coherence envelope; independent proof audits confirmed both. The latter is
the existing C59 Cauchy-Schwarz mechanism, not a free efficient sampler.
The dual pilot's first passing flag hid a shared FFT normalization error;
ER records the correction and the stronger normalized/geometric audit.
The exact C79 regressions cover residue partitions, revival and support
separation. The production helper remains C78-only; TODO 39 owns next work.
No manuscripts, abstracts, defaults, host settings, firmware or commits
changed. The user's open-ended research goal remains active.

**Earlier-phase sampler checkpoint completed: C79 / ER / TODO 39.** Main
implemented EarlierPhaseProgressions by specializing the C78 helper, keeping
the same actual work draw/rejection/progression loop and exposing the chosen
row stride. Lower-cost testers supplied joint-law, RNG/edge and long-row
experiments; main corrected verifier normalization, insertion indexing,
nonvacuous controls, cumulative counters and allocation lifetimes, then
independently reproduced them. ER owns reports, failure provenance and limits.
The new permanent lab regression enumerates the weighted actual tiny RNG
transition, including both gcd lifts and the work-redraw bias control.

Core/lab/claims pass; lab passes both with and without the optional backend.
The other six science suites were not rerun. Documentation validation is
recorded in ER. The opt-in prototype is not a finite-bit certificate or a
generic simulation/timing breakthrough. Its C59 mass-weighted proposal
follow-up is now complete at the checkpoint below. No manuscript, abstract, existing default, host
setting, firmware or commit changed. The user's open-ended goal stays active.

**Weighted-proposal checkpoint completed: C80 / CW / TODO39.** Main added
root_mass as an opt-in mode in the shared conditional loop, preserving mass
as default. Lower-cost testers supplied the tiny and physical checks; main
strengthened normalized accepted-law/RNG predicates, failure accounting and
allocation lifetimes, then independently reproduced both. C80 owns the proof
and limitations; CW owns authoritative reports, observed named costs, failure
provenance and scoped gates. Core/lab/claims pass, with lab tested both with
and without the optional backend; the other six suites were not rerun.

**Nested-phase checkpoint completed: TODO40 / C81 / NC.** Main derived the
hybrid coefficient construction and largest-effective-gap bound. Independent
lower-cost proof/integer/amplitude testers supplied initial checks; main fixed
period handling, vacuous/missing controls and work/allocation accounting before
independently reproducing tiny, edge and selected long-row comparisons. The
long-row result saves actual construction terms between exact covers of the
same state, not just an untruncated upper bound. The strongest generic
fixed-r sequential/table baselines are not beaten by this observation.

Main then derived exact finite-truncated T/S counts, independently audited
and enumerated, which select the stated naive construction without building
competing rows or querying phases. C81 owns the proof and its cost-model
limits; NC owns all authoritative reports, numerical witnesses and failures.
A nearly canceled edge row remains a negative conditional-normalization
control, not a silently excluded positive test or a significant joint-law
error. Production helpers still support only C78/C79 schedules.

Core ran first, claims (including the new selector regression) and the legacy
ER reference experiment pass; the other seven science suites were not rerun.
NC records documentation validation. TODO41 alone owns sampler integration,
actual RNG/retry validation and the matched end-to-end comparison. Do not
restart the completed integer/amplitude pilots. Manuscripts, abstracts,
existing defaults, host settings and commits remain untouched. The user's
open-ended research goal remains ACTIVE; this turn made concrete PROGRESS.

**Nested-phase sampler checkpoint completed: TODO41 / C81 / NS.** Main
implemented the opt-in multi-phase helper and audited lower-cost test drafts.
The shared sampling methods/defaults are unchanged. NS owns complete tiny-law
and weighted actual-RNG checks, the matched returned-sample comparison, charged
failure/selector/cache work, numerical caveats and verifier corrections. The
construction saving survives rejection in the bounded fixture; this is not
a native-time win over the fixed-r sequential baseline or a numerical
certificate. TODO42 alone owns the next selector question. No manuscripts,
abstracts, host settings or commits changed. The broader research remains open.

## Paper-refinement context

Decisions taken 2026-08-08, so they need not be re-asked:

1. **Venue.** Shape both papers **as if for a journal** — focused, every
   paragraph earning its place — but the likely first destination is **arXiv**.
   The user has not published before.
2. **Order.** Consistency pass across both (**done**), then Paper A end to end,
   then Paper B.
3. **The third paper** (compilation-dependence) stays **scattered as caveats
   inside A and B**; revisit after the first two are submitted.

**Done so far:** a consistency pass that fixed sixteen defects, and a
structural pass on Paper A that fixed a duplicate Proposition number, cut §9,
and removed a note-to-self from the shipped text. See commits `184f49c`,
`d68eca3`, `0f2bd6e`.

**Still open on Paper A**, and deliberately left for the user:

- **§10.1 / §10.2** (limitations, retractions) are load-bearing and currently
  in-body. Whether they stay there at journal length is an unmade decision.
- The §3.2 / §4.1 / §6.3 trio explains the same invariance relationship three
  times. Defensible signposting, but it is where a "focused" edit would bite.
- **|B| = 2 is measured, not proved** (Paper A §5). The paper says so. Keep it
  that way unless someone proves it.
- Consult the current Appendix A maps before relocating claims; several maps
  and proposition numbers changed during the independent review.

---

## Do NOT redo these

Listed here so a fresh session does not burn time; each has a retracted claim
file or a note carrying the full reasoning.

- **The compilation thesis** ("Toffoli vs Fourier compilation determines
  cost"). Dead twice over — a propagator bug plus an ancilla-count confound.
  Both compilations have exact Z-closure and the identity is
  compilation-invariant.
- **"The QFT is the PPS bottleneck."** Dead — bit-reversal bug plus a toy that
  was not a proxy for Shor.
- **The toy "Shor sandwich"** (H → adder → QFT). Not a proxy. Use `modexp.py`
  or `toffoli_arith.py`.
- **Affineness of the identity block** as a route to the C15 proof. Refuted
  (14336/32768 violations). The proof goes via the parity reduction instead.
- **The (z, z⊕e) pairing** as the C15 mechanism. An n_exp=1 artifact.
- **An intermediate 2-adic law in α.** There is none; the β=1/β>1 split is the
  real structure and the apparent intermediate case was sampling aliasing.
- **The msb↔anc CNOT pairing as the cause of the linear structure.** Refuted;
  see §L2 for the actual three-part mechanism.
- **A "cliff" in the hyperplane deficit between n_exp = 2 and 3.** Refuted the
  same day it was conjectured (§DF); D = 1 − 2·density is monotone decreasing
  in every β > 1 instance.
- **The §I period-ord₂(β) oscillation at CIRCUIT level.** Function-level only;
  it does not transfer. That import is what made 12g look worth a sweep.
- **The identity string Z⁰, or a parity constraint, as the peak deficit.** Both
  were the named candidates in TODO 12d and both are wrong; the measured
  arithmetic deficit is associated with two dominant modes (§PK). Do not
  restore the withdrawn universal peak formula; see C44.
- **Cîrstoiu for TODO 14.** Pulled and read; their group indexes *circuit
  parameters* and their results concern *ensembles*, while Shor's circuit is
  fixed. Useful for framing, not for the technical problem.

---

## Where everything is

| you want | read |
|---|---|
| conventions, traps, discipline | `CLAUDE.md` |
| what is known | `claims/INDEX.md` → `claims/<ID>.md` |
| what is next | `todo/INDEX.md` → `todo/open/*.md` |
| how something was found | `notes/INDEX.md` → `notes/<CODE>-*.md` |
| how the work is conducted | `METHOD.md`, `AGENTS.md`, and the installed `qsim-research` skill |
| why the repo is laid out this way | `RESTRUCTURE.md` |
| GPU capacity and measurements | `notes/GPU-*.md` |
| the resolved 0.716 thread | `notes/HIST-*.md` |

`CLAIMS.md`, `NOTES.md` and every `INDEX.md` are **generated** — edit the
individual files, then run `tools/reindex.py`. `tools/check.py` must pass
before any documentation commit.
