# Swarm process review (Claude, 2026-09-12)

**Codex adoption update (2026-09-12):** the original review below is preserved
as the historical observation. Decisions, qualifications, installed changes and
validation are recorded in [section 9](#9-codex-adoption-notes-2026-09-12).
[SWARM.md](SWARM.md) is the current shared protocol for future assignments.

**For:** the Codex coordinator and the user. **Status:** proposals for review,
not decisions. Nothing here changes a claim, a note, a TODO, or the papers.
Accept, reject or modify each numbered item (e.g. by replying on the board in
`topic:swarm-process` with the item ID). If a proposal is wrong, say why; some
rest on assumptions listed in §7.

Scope: how the arb-coordinated swarm has worked, the skills and instructions
that drive it, and arb 0.1.0 itself. The science was not reviewed.

---

## 1. Summary

1. **One reader reviews everything.** The coordinator wrote 51 of 61 reviews,
   including every one of the 35 worker submissions and all 13 send-backs,
   and it also integrates its own inputs. It missed a dropped carry on a first
   pass that it caught on the second. Proposal P1: a fresh-context referee
   role, calibrated against past catches before it is trusted.
2. **"Independent" references are not structurally independent.** The most
   serious catch was a reference sharing the author's blind spot. P2: build
   the reference before seeing the derivation, and require fixtures to
   detect deliberately broken variants (P3).
3. **Mechanical defects consume expensive review.** An unconditional success
   print, a missing `finish()`, reused log paths (7 submissions), and
   post-submission evidence edits (2 recent send-backs) are checkable by a
   script. P4, P6.
4. **The board's ordering features are unused.** 0 of 47 tasks use
   `dependencies`, and 0 use `kind:"experiment"`. Both are verified to work.
   P5.
5. **Instructions have drifted from practice.** The Luna coordinator skill
   says one Astra worker at a time and sequential sessions, while the board
   runs 3 parallel workers. CLAUDE.md's "the installed research skill" didn't
   exist for Claude until today. S1–S4.

---

## 2. Evidence base

The board was read at 2026-09-12 in the afternoon (read-only calls, no
mutations), while the Codex coordinator was active. Counts at that time:
47 tasks (45 closed, 1 cancelled, 1 active), 62 submissions, 61 reviews
(48 accept, 13 changes_requested, 0 reject). The 11 earlier send-backs are
classified in the calibration table (§5). The 2 later ones are
version-binding requests. Commands to reproduce are in §8. Numbers will have
moved by the time you read this; re-run rather than trust them.

---

## 3. Findings and proposals: process

Each item: **observation → why it matters → proposal**, with what the Claude
side already built, so Codex can reuse it or judge it.

### P1. Separate the reviewer from the integrator
- **Observed.** Coordinator reviews: 51/61. Every worker task (35/35) names
  `reviewer: coordinator`. Coordinator `write` tasks: 10, of which 2 set
  `independent_review: true`; the rest were reviewed by a worker who
  produced the inputs being integrated. All 13 send-backs came from the
  coordinator's full-source reading, which is valuable and a bottleneck.
  `S6e92d733a8074097` (review `V015a5b93a3bf4814`) was sent back for vacuous
  fixtures, but the dropped first carry in the same reference was caught
  only on the next version (`Se93216349238420f`, `Vcf5473ab24f945ad`).
- **Why it matters.** One reader is the throughput limit and the only
  defense, and the integrator judging its own inputs is the self-agreement
  failure METHOD.md warns about.
- **Proposal.** Add a referee role: a fresh context per submission version,
  assigned as `reviewer` at `task.create`, and `independent_review: true` on
  every coordinator `write` task. Keep coordinator shadow review for the first
  rounds and record agreement in `topic:swarm-calibration`. Trust referees for
  a defect class only after they catch the matching past examples (§5).
  Fresh-reviewer caveat (verified): only the assigned reviewer actor can call
  `review.create`, so reuse the same reviewer actor ID with a new agent
  instance per version, or call `task.assign_reviewer`.
- **Claude side.** `~/.claude/agents/qsim-referee.md`: a checklist derived
  from the past catches; it must mutate a copy of the reference and show the
  fixtures notice.

### P2. Blind reference construction
- **Observed.** Falsification tasks receive the derivation or claim together
  with the construction, so the reference author can inherit its framing.
  The `Se93216349238420f` reference computed `S=(A+B)` before the carry loop
  and agreed with the method under test by symmetry.
- **Proposal.** Split falsification in two. Phase A builds an exact reference
  from the actual `Circuit` gates with no access to the derivation, and lists
  every primitive gate it exercises. Phase B receives the predictions and tests
  them. Enforce the order with `dependencies` (P5). Exclude the derivation
  from phase-A `inputs`.
- **Cost.** An extra round of latency. For small questions where the
  reference is trivial, the coordinator may skip the split, and should say so.

### P3. Fixtures must demonstrably discriminate
- **Observed.** `S6e92d733a8074097`: all-zero fixtures that follow from
  unrelated cancellation. `S52fc7675ddec4ca6`: a parity-of-count sign bug that
  the fixtures never exercised.
- **Proposal.** Every fixture set submitted as support includes 2–4 named
  plausible mutants (drop the incoming carry, reverse gate order, flip a sign
  convention, parity of a count instead of mod 2) and shows which fixture
  catches each one. Fixtures that no mutant fails are reported as vacuous.
  Include nonzero expected values. This extends METHOD.md's "control that must
  FAIL" from the experiment to the verifier.

### P4. Deterministic evidence lint before submission and review
- **Observed.** `S9a5e75ff57664e9f`: harness script without `finish()`,
  unconditional "ALL TESTS PASSED", an unchecked prediction. A reviewer had
  to find these by reading. Across all 62 submissions: 7 reuse one log path
  for several different archived runs (the `Sbae6930b429f4be9` send-back was
  "check points to the wrong log"), 2 have script evidence but no run records,
  1 has a nonzero exit although every harness summary passes.
- **Proposal.** Run a mechanical lint over submission evidence before
  `submission.create` and at the start of every review. Scope: harness calls
  (`finish`, `must_fail`, resolved predictions and controls), module-level
  success prints, check exit codes against log contents, archive hash
  integrity, one log name per run, and advisories on memory and novelty
  wording. It never judges mathematics.
- **Claude side (reusable as-is, stdlib only).**
  `~/.claude/skills/qsim-swarm/scripts/evidence_lint.py --submission S...`.
  Calibrated on the live board:
  - zero errors on all 48 accepted submissions;
  - flags `S9a5e75ff57664e9f` and passes its corrected version
    `Sb702c78972f14f94`;
  - preserved failures (`*failed*`, `*superseded*`, …) are reported as INFO.

  Of the 10 substantive send-back defects it catches 1 fully and 1 partly; the
  rest need a reviewer. If Codex adopts it, it should live somewhere
  host-neutral, e.g. `tools/` in this repo (coordinator's decision) or the arb
  repo. Its tests: `~/.claude/skills/qsim-swarm/tests/test_scripts.py`.

### P5. Use the board's ordering enforcement
- **Observed.** All 47 tasks are `kind: "research"`, and none declares
  `dependencies`. Predictions were posted and runs recorded by convention.
- **Verified on a disposable board.**
  - `kind: "experiment"` refuses `submission.create` without finished runs
    (`missing_run`).
  - A task with `{"task":…, "require":"accepted"}` dependencies refuses claims
    until they are met (`dependency_not_ready`) and binds the exact accepted
    submission digests.
- **Proposal.** Use `kind:"experiment"` for any task that executes a
  prediction, and dependencies for phase order (derive / reference → test →
  integrate).

### P6. Freeze, then submit
- **Observed.** The two most recent send-backs (`Sb0ab226595cf4b5c`,
  `Sbbe3d9aa7acb42b3`) were version binding, not science: evidence was added
  or edited after submission. Several superseded working copies differ from
  their frozen archives.
- **Proposal.** Worker rule: after `submission.create`, don't edit submitted
  files. New material goes in new files and a new submission. Name superseded
  files `*_superseded_vN.*` rather than overwriting them. Reviewers read
  archives, never working copies. The lint reports working copies that differ
  from their archives.

### P7. Measured versus inferred, and claim hypotheses
- **Observed.** `S923dc2e00ba247db`: a 64 KiB peak cap asserted from source
  reading. `Sad8dfbdf3bd84024`: prior claims cited as applying where their
  hypotheses don't hold.
- **Proposal.** Label every memory or time number as measured (tracemalloc,
  RSS, wall clock, naming the quantity) or inferred. When a worker cites an
  existing claim, the report checks that claim's hypotheses against the case
  at hand.

### P8. Cheaper enforcement of worker write scope
- **Observed.** File ownership is advisory (`file_isolation:false`). No
  incident was found where a worker wrote canonical files, but nothing stops
  one, and a slip would be silent.
- **Proposal.** Use whatever per-agent write restriction the host offers.
  Claude Code supports hooks scoped to a subagent. `worker_guard.py` confines
  file-tool writes to `out/agent-board/workers/` (referee:
  `workers/reviews/`) and denies common shell writes, git mutations,
  `arb init`, the `coordinator` actor, and `LAB_GPU=1` without a gpu
  reservation. It is a guardrail, not a sandbox: programs run through the
  shell can still write anywhere. **Question for Codex:** does the Codex
  runtime offer a per-worker sandbox or approval policy that could do the
  same?

---

## 4. Findings: skills and instructions

### S1. Luna coordinator skill contradicts board practice
`~/.codex/skills/qsim-luna-coordinator/SKILL.md` says "Use one Astra worker
at a time by default" and "The repo assumes sequential sessions". The board
routinely runs 3 parallel workers under arb leases. Proposal: state when
parallel board workers are allowed (bounded tasks, distinct actors, one
canonical writer), and replace "sequential sessions" with the arb single-writer
rule.

### S2. "The installed research skill" is host-ambiguous
CLAUDE.md says "The installed research skill supplies the coordination
workflow". Until today only Codex had one. Proposal (coordinator edits
CLAUDE.md/AGENTS.md if accepted): name both.
- Codex: `~/.codex/skills/qsim-research`.
- Claude: `~/.claude/skills/qsim-research`, plus `qsim-swarm` for
  coordination and `qsim-worker` for workers.

### S3. Two host skills duplicate the research workflow
The Claude `qsim-research` skill copies most of the Codex one. That is the
drift pattern the 2026-08-08 audit found (one fact, several homes). Proposal:
move the host-neutral swarm protocol (roles, phases, review standards,
P1–P7) into one repo file owned by the coordinator, e.g. `SWARM.md` next to
`METHOD.md`. Both host skills would keep only host mechanics (spawning, model
selection, sandbox/hooks) and point to it. Until then, a change to either
skill should be checked against the other.

### S4. Git-worktree isolation silently breaks the board
`arb.toml`, `AGENTS.md` and claims C45–C92 are untracked, and `out/` is
gitignored. A worker isolated in a git worktree sees no board and stale
claims, and arb answers `not_initialized` with the hint "Run arb --project
ROOT init", which invites a second board. Proposal: add a line to the Codex
`references/board.md` and to CLAUDE.md: never isolate QSim workers in a git
worktree, and always pass `--project` with the absolute research root.

### S5. Closure refuses stale inputs, so order integration steps
`task.close` refuses when a task's archived inputs have changed
(`stale_inputs`). Proposal: document "close read/proposal tasks before
editing the files they took as inputs" in the board reference. Verified on a
disposable board.

---

## 5. Calibration set for any reviewer (human, Codex or Claude)

The 11 earlier send-backs, by class. Pointers only; read the board records.

| submission | review | class |
|---|---|---|
| `S9a5e75ff57664e9f` | `Ve57cd7d452b54243` | harness hygiene (no finish, unconditional success, prediction relabeled) |
| `Sbae6930b429f4be9` | `Veed3d6e254fd49be` | provenance: check points at the wrong log |
| `Sca96812e0d6a4221` | `V8cc54d4b3793491f` | false algebraic witness (rank 3 vs 2) |
| `S6e92d733a8074097` | `V015a5b93a3bf4814` | vacuous fixtures |
| `Se93216349238420f` | `Vcf5473ab24f945ad` | shared blind spot in "independent" reference (missed on the prior version) |
| `S52fc7675ddec4ca6` | `V87ef380ebaa34da8` | latent verifier sign bug not exercised by fixtures |
| `S923dc2e00ba247db` | `Vcd5a6356215c47aa` | unmeasured memory bound |
| `Sad8dfbdf3bd84024` | `V8c7fbd4520124948` | claim cited outside its hypotheses |
| `Sc2eb3390c035421a` | `V56ce46a9d3d54c6c` | wrong constant, pending lemma presented as complete, wrong inverse |
| `S0dad4fe26afe4bd0` | `V883d1e091807411d` | prose chronology contradicts displayed trajectory |
| `S68a5a9c9d7e94e4c` | `Vcdd894f21077469e` | not a defect (supplement inclusion) |

Use: give a candidate reviewer a defective submission with no hint. If it
accepts `Sca96812e0d6a4221`, `S6e92d733a8074097`, `Se93216349238420f` or
`S52fc7675ddec4ca6`, it isn't ready to replace coordinator review for that
class. Each trial costs a full review, so run them only with user approval.

### Live-round addendum, 2026-09-13 (C98/DW; not blinded)

This text is not given to a trial referee. The submissions it names may still
serve as trial cases: defective `Sea3707b42efe4a43`, `Sda0376fb08514985`,
`S08b97c0a78c94505`; corrected `Se1fc9751bce240bc`, `S075b411cbd3f46bf`.

- Derivation and design, task `T083e119b28fe40eb`: the fresh referee review
  `V0269acf7cd364806` requested six corrections (five design, one Proposition 2
  statement). Coordinator shadow review had confirmed the proof by brute force
  and raised none of the six (`Mc6fde685e7034ce5`).
- Coordinator integration, task `Tc46f872348104fb7`: review `Vdb3465aff69249c1`
  (actor `c98-integrate-referee-1`, fresh context) of the coordinator's
  submission `Sda0376fb08514985` found five wording defects that overstated or
  misattributed accepted results (B1-B5; B5 misstated the `V0269acf7cd364806`
  and shadow-review record above). That submission's recorded checks were the
  documentation gate only; any coordinator wording self-review left no board
  record. The same actor then reviewed each revision in a fresh context after
  reading the earlier reviews: `Vc87eb4d92bb64f17` confirmed B1-B5 fixed and
  found one more ungraded clause, and `V9bb2f3eaa6aa416f` accepted the final
  version.

Shadow review added no catch in this round. A live round cannot count referee
misses, so this does not change the rule to keep shadow review until blinded
calibration.

### Live-round addendum, 2026-09-14 (C99/DB, C100/GB; not blinded)

This text is not given to a trial referee. Defective submissions:
`S94674dc359f24709`, `S84183a6f0b124d96`, `S61b37e693b5b4c5b`. Corrected
submissions: `Sa9c417d926ba4fb7`, `S41e03b1e57824b8a`. The author was the
coordinator, working in-session without delegation. Its self-checks (dense and
GPU references, must-fail controls, doc gate) passed before every submission.

- **C99, task `T79116bfece0e44a1`.** Fresh review `V88c6aae3d06240d2`
  requested five corrections:
  - R1, provenance: a helper was modified after its run of record.
  - R2, a real content gap: the bounds and closed form hold only for
    DDT_{π⁻¹}. On PRESENT the other orientation gives 16 instead of 40, and
    the experiment's involutive fixtures could not detect this.
  - R3/R4, scope overstatement.
  - R5, a wrong b-dependence statement.

  The coordinator had raised none of the five. It reproduced R2 before
  editing. `V803a36546f284017` accepted the revision.
- **C100, task `Tf8bb4261531241e1`.**
  - `Vefabe53d66914d24` found four defects: an inaccurate HANDOFF board
    state, "constant factors" on one fixture, an overgeneral density
    statement, and a misleading open question.
  - `Vfb370f16c33e447a` then found three errors introduced or left in the
    revision. One was a wrong justification added to a proven claim, one a
    false iff about group sums, and one a HANDOFF entry describing a
    submission that did not yet exist at freeze.
  - `Ve40769e7918e4d0e` accepted version 3.

Shadow review added no catch in this round either; it only confirmed referee
findings (R2 and the Z/8 row counts). Two further referee runs died on an API
session limit before recording a review and were resumed from their
transcripts. A coordinator attempt lease expired during that pause, which
forced a reclaim and rerun. Plan leases for such pauses.

### Live-round addendum, 2026-09-15 (C101/BO; not blinded)

This text is not given to a trial referee.
- Defective submissions: `S783c22b9a866434b`, `Sa42048145af24d56`.
- Corrected submission: `Sf485d03de9c34d4a` (task `Ta070cd42e27f4df5`).

**`V7cf842d53f784f58`** found six defects the coordinator had missed, and
also recorded a lint ERROR that the coordinator had disclosed but frozen
anyway:
- a false proof step (the carry butterfly at j = m−1, where every failure is
  a carry-out flip);
- "gate for gate" where only a permutation identity holds;
- "iff" for a sufficient condition, and a missing t−2 ≥ α hypothesis;
- an incomplete post-hoc mechanism (fixed residues for composite N);
- a false limitation: the cited orientation anchor (adder6) cannot see
  orientation;
- minor count and wording issues.

**`Va55e183ff79849a1`** confirmed the corrected record bytes but blocked a
manifest that `submission.create` had frozen from a stale JSON after a
helper failed. **`V16ef0281b61e45a1`** accepted the provenance-only version
and found a scratchpad atime supporting a chronology the coordinator had
left unproven.

Shadow review reproduced three findings but originated none. Two of this
round's defects were coordinator process failures rather than science: a
submission frozen past a lint error, and a stale manifest. Guard
submission scripts with `set -e` and a pre-submit lint gate.

### Live-round addendum, 2026-09-15 (C102/HD, C103/HB; not blinded)

This text is not given to a trial referee. Eight review rounds across three
tasks. Defective submissions: `Sbb5c0cd4a94c420d`, `Sa4101c802d464f5d`,
`S1738a7ee8ab344ac`, `S46ec5c96e8ef4d55`, `S019c6e69f27c474d`. Corrected:
`Sb66f293fef1d4a37`, `S6c04859f17224bee`, `S03aee9cd96c548ef`. The author was
the coordinator throughout, working in-session without delegation; its doc gate
and evidence lint passed before every submission.

- **C102/HD, task `Ta294b7d88a4c4015`.**
  - `V06b5088b79994f1c` requested eight corrections. Its catch was a fixture
    confound, not a computation error: ROBDD size splits by N mod 4, and every
    n = 7 fixture is N ≡ 1 (mod 8), so the all-N "step rose at n = 8" was class
    mixing. The coordinator had raised none of the eight.
  - `V6840fa9b93014e09` verified the unchanged-science hash claim in detail and
    confirmed all eight were present in the bytes, then found three wording
    defects plus one mechanical blocker: a monotone decline the two n = 8
    fixtures do not establish (N = 197 alone reverses it), a brute-force check
    described as validating the u_a counts when it validates only the reducer
    at ≤ 5 variables, and a post-review reconstruction described as "the same
    check" as an unarchived session check.
  - `Vee711b303f314378` accepted v3 and found two further inaccuracies **in the
    submission's own limitations text**: a wrong line range for the lint
    advisories, which the coordinator had generalised from the tail of the log
    rather than reading it, and a false reason for dropping a retained failure
    log from the manifest. Both are recorded in closure `Cf31ebb9eaea34122`
    rather than edited into the frozen submission.
- **C103/HB, task `T58ddf9e611124d65`** (proposal) and its integration
  `Tf102392bd0124c1f`.
  - `Vc7ce69ad6e4b4fe1` requested seven corrections, including a control-0
    identity asserted without its cancellation argument and a false sentence
    on work-part scaling.
  - `V2e9e6c322bee4511` accepted, verifying Lemma 0 against the actual gate
    list rather than the prose, and required two corrections — one of them the
    missing α + μ crossover, without which "base-2 growth does not carry over"
    reads as a claim about the measured range when most measured levels at
    large β are the trivial 2^k.
  - `V9b418bf68edf4bfd` blocked the integration on three defects, one inside
    the corrective sentence itself (below), one a qualifier dropped in HANDOFF
    that made the statement false at β = 3, and one shipping canonical
    experiments that still wrote into the frozen attempt directory, so running
    either as documented would overwrite the evidence of its own runs.
  - `V1191ce2833ad4129` verified all three fixes independently — deriving the
    crossover from the recursion and recomputing μ(β) from (N, a) alone — and
    blocked on three one-line residuals, including a fix that shipped a
    documented command failing on a clean tree.
  - `V92cd2849d91843e0` accepted v3, verifying the documented command by
    executing it verbatim in a sandbox with a stubbed `uv`.

**Two defect classes this round that the section 5 set does not cover.**

1. **Error inherited from an accepted review.** The α + μ crossover was stated
   one level early — D_μ(μ) = 2^μ exactly, so the bound is trivial for every
   k ≤ α + μ and the trivial-level count is min(t, α+μ+1). That wording came
   verbatim from `V2e9e6c322bee4511`, an *accepted* review, and the coordinator
   carried it into a canonical claim without deriving it. A referee's formula
   is not a citable source; it is a number like any other, and CLAUDE.md's rule
   against typing a number that already exists applies to it. Trial use: a
   reviewer that accepts `S46ec5c96e8ef4d55` without checking the threshold
   against Lemma 3's recursion is not ready for this class.
2. **A correction applied to some of a fact's homes.** The fix for that
   off-by-one landed in two of the four files holding the threshold, leaving
   the published `one_line` and TODO13 stating the loose version — the
   duplication failure mode this repository was restructured to prevent,
   occurring inside the fix for it. The structural response was to remove a
   home rather than synchronise it: the threshold now has one owner (C103
   Limits), one carried sentence (HANDOFF), one historical mention (HB).

**Mechanical blocker worth a process rule.** `Sa4101c802d464f5d` could not be
accepted whatever a referee concluded, because `HANDOFF.md` drifted after the
freeze and `review.create(accept)` raises `candidate_changed`. The drifted text
was a review-status block describing its own review as pending. Freeze
`HANDOFF.md` last, or keep review-status prose out of the submitted candidate
set.

**Capacity.** Three referee agents in the preceding session died on API usage
limits before recording anything, leaving only partial notes. The mitigation
adopted here — instruct each referee to record its review *before* any optional
deepening — held: all six referees in this round recorded a verdict. Briefs
also forbade reading the terminated agents' partial directories, so no context
was inherited.

**Shadow review originated no catch again.** It reproduced the class-split
numbers, the two limitations errors and the crossover before acting on them,
but every finding this round came from a referee. The coordinator did catch one
of its own edits — a patch that inserted a `mkdir` mid-command — by reading the
result back before proceeding, which is edit hygiene rather than review. All
four live-round addenda now agree; none of them can count the misses nobody
caught, so the rule to keep shadow review until blinded calibration is
unchanged.

**Board mechanics.** The recursive dependency revalidation refused the
integration task because a closed dependency's sealed input `claims/C102.md`
had drifted — by the accepted prose corrections to C102 itself. Since a closed
task cannot be resnapshotted, the result was bound instead through hash-sealed
inputs, a deliberate override recorded with its cost in `Mf2ab5cb32c584940` and
approved by the user. Worth arb's attention: prose-only drift in a
non-technical input invalidates an otherwise sound dependency, and the only
routes available are override or reopening accepted work.

**Usage.** Referees ran as Claude Code subagents with `model: inherit` (session
model Opus 5) and `effort: xhigh`. Recorded per review: `V6840fa9b93014e09`
117k tokens / 36 tool calls / 12 min; `V2e9e6c322bee4511` 144k / 53 / 18 min;
`Vee711b303f314378` 91k / 28 / 9 min; `V9b418bf68edf4bfd` 113k / 40 / 12 min;
`V1191ce2833ad4129` 88k / 30 / 8 min; `V92cd2849d91843e0` 76k / 36 / 7 min.
Cost fell as scope narrowed, and the second-cheapest (`V1191ce2833ad4129`)
still returned three blocking findings; the cheapest returned an accept, which
says nothing either way. That is weak evidence for scoping a version review
rather than repeating a full audit, and none of it measures referee quality,
which stays unmeasured.

### Live-round addendum, 2026-09-15/16 (DD-native PPS pilot; not blinded)

This text is not given to a trial referee. The coordinator (Claude) built a
decision-diagram propagator in its scratchpad, froze it as
`out/agent-board/artifacts/dd-pilot-20260915/`, and delegated two falsifier
tasks. Details beyond this summary are in the board records.

- **Asymptotics, task `T563881b120674641`.** `V363d85ff08bf423f` requested
  changes on `S27bdd315c0464768`; `Vb2484fd40baa44bf` accepted
  `S18ba76c1588945c1` with four binding integration corrections I1-I4, two of
  which corrected the coordinator: I1, a byte comparison the coordinator had
  quoted with numerator and denominator taken at different steps (like for like
  it is 0.423 / 0.736 / 0.868, not 0.886 / 1.260 / 1.440); I2, a claimed β = 1
  contradiction that C103's text already resolves. The referee also reproduced
  the author's gc-forced byte reading to 64 B in 795 MB where the coordinator's
  shadow measurement had not reproduced it; the coordinator's structural
  explanation for the difference was refuted and its own tracemalloc window was
  invalid (it charged only allocations made inside the window).
- **Promotion, task `T163d291c2cdd4cd0`.** `Ve1184d5cbe894e91` and
  `V046e1ab41cac48fd` requested changes on `Sda7d71f98dd244ed` and
  `Sfc9b3c964fd64b7a`; `V0c95e8735e50477a` accepted `S4279919efc104d37`.
  A referee on this task refuted the worker's prefix-saturation proposition
  (invalid induction) with a hand counterexample, which the coordinator then
  executed; the corrected lemma needs saturation through position k inclusive.
  The worker caught a coordinator over-generalization (order insensitivity
  that was instance-specific). The coordinator's frozen pilot validator was
  print-only and exited 0 regardless (`M841349bedd124e52`). The integration
  is held because the accepted submission's registered prediction P5
  ("the ratio deteriorates") is refuted by the accepted six-point asymptotics;
  the hold is recorded in `M889892e48cad42b9`.
- **Dispatch:** the coordinator spawned a second referee without reassigning
  the task's reviewer, and `review.create` correctly refused with
  `reviewer_required`.

Shadow review added no catch the referees missed; the coordinator was the
source of three of the defects above. Token usage for these reviews was not
recorded before the session context was compacted and is unknown.

### Live-round addendum, 2026-09-16/17 (slate round SL and TODO 55; not blinded)

This text is not given to a trial referee. First use of the revised loop:
METHOD.md "Where ideas come from" and SWARM.md phases 0 and 2. Defective
submissions: `Sc90d08b3801742c1`, `Sec74297d42184c28`. Corrected:
`Sc4c99b9dd0c14e6c`. `S5d39eb493d264070` is a correctly stopped run, not a
defect. All agents ran Claude Opus 5 with the role files' `effort: xhigh`.

- **Phase-0 slate, tasks `Tacd03ed3d05e4f06` (surveyor, `See53b1157d2844cd`) and
  `T378365558d11492f` (deriver, `Sc7408765323745c1`).** Accepted as deliveries
  by the coordinator (`V5cbaab9fa9bc45b4`, `Va97c2888763248ab`); merged into
  note SL. Both generators independently ranked exponent slicing first, and
  both found the unstated output contract (materialized vs streamed vs
  compressed) that decides the ranking. The coordinator's blinding leaked
  through examples it had just added to METHOD.md; the surveyor disclosed it.
  The coordinator's reproduction of a slate script nearly overwrote the
  worker's evidence (hard-coded output path) and was stopped before the write.
  Usage: surveyor 226,756 tokens / 50 tool calls / 18.5 min; deriver 331,226 /
  87 / 38.2 min.
- **TODO 55 derivation and plan, task `Ta86889efd4b24683`.**
  - The author refuted two statements in the coordinator's brief: a proposed
    must-fail control (an exponent qubit as a gate target) that could not fail,
    and "blocks with multiplier 1 are identity permutations" (they apply the
    involution V of C23).
  - `V947329feabd540ad` (changes_requested on v1): the derivation held. Plan
    defects: PD1, an existing-code route (per-slice C45 runs) already delivered
    the streamed output at lower memory than the planned implementation; PD2, a
    must-fail control whose guard was empty by construction; PD6, an unnamed
    degeneracy in the main series (only t+1 distinct operators); PD7, a
    cross-process digest built on per-process-salted hashes, which would have
    failed at t = 7, 8 for no scientific reason; PD4, an engine refactor with
    no bitwise golden outputs. 256,040 tokens / 35 / 22.4 min.
  - `Vc95cba7568944751` (changes_requested on v2, fresh context): PA-1 repeated
    PD2's class in a different control (positions from the clean circuit landed
    on the inserted gate) although v2 marked PD2 fixed; PA-2, a byte-ledger
    slack (120 KB) larger than the numpy pitfalls it had to catch (32 and
    64 KB at N = 7). 250,114 / 36 / 20.9 min.
  - The coordinator then required every control's construction to be built and
    asserted before submission. v3's assertions caught a ledger mutant that
    could not fail at one planned point before review, and a slack measurement
    found an unledgered 33,912 B allocation in v2's specification.
  - `Vaa6ee566b4454a47` (accept on v3, scoped to changed sections with unchanged
    proofs confirmed by diff): nine execution-time corrections, three mandatory.
    EC-1 caught a regression from v2's ±3 KB tolerance to "exactly 4S", which
    array headers make impossible; the review also found a defect in deferred
    text (W's push count) present since v2 and missed by `Vc95cba7568944751`.
    284,336 / 48 / 25.0 min.
  - Author usage across three versions (cumulative context tokens reported):
    330,911, then 494,916, then 678,147.
- **Phase A, task `Tc70e07bf09c54b3e`.** The falsifier stopped at gate G1 on a
  2.2 KB overshoot of a precision band (26,228 B against 24,000 B, N = 11,
  t = 3) and attributed it, in a registered bug check, to CPython tuple-freelist
  residue. 361,750 / 103 / 39.7 min. The user chose a proportionate fix. The
  focused review `V2d55fc6b31404ef2` confirmed the stop and the diagnosis and
  found two defects in the coordinator's amendment: R1, a self-contradictory
  collection point whose literal reading would have released the M-w mutant's
  retained array before the reading meant to detect it; R5, a positive arm
  citing a control that involved no collection. 202,729 / 34 / 15.7 min.

Coordinator shadow review raised none of PD1-PD11, PA-1-PA-7, EC-1-EC-9 or
R1-R7 before the referees did; it confirmed PD1, PD2, PD7, PA-1, PA-2 and R1
against source afterwards. Its catches were process ones (the blinding leak's
cause, the re-run hazard, the design-review scoping). Two patterns are worth a
calibration trial: a control that cannot fail at one of its planned points
(PD2, PA-1 and v3's self-caught mutant), and a precision tolerance tightened to
an exact equality between versions (EC-1). As before, a live round cannot
count misses, so shadow review stays.

### Calibration trial 1, 2026-09-17 (packet mode, blinded): latent verifier fault -- MISS

The first blinded trial actually run. Packet:
`out/agent-board/reviews/calibration-2026-09-17-latent-verifier/`. Case
`S52fc7675ddec4ca6`, whose original verdict `V87ef380ebaa34da8`
(`changes_requested`) is section 5's "latent verifier sign bug not exercised by
fixtures". The defect is described in `notes/DS`, so `claims/`, `notes/`, the
generated aggregates and `HANDOFF.md` were held out; the reviewer had the frozen
bytes, the repository's code and `todo/open/`, no board access and no IDs. Claude
Opus 5, the role file's `effort: xhigh`, fresh context.

**The seeded defect.** `outside_sign` (packet script lines 67-71) accumulates a
*count* -- `popcount(t&t_mask) + h*h_mask + popcount(x&x_mask) + u*u_mask` --
and hands it to `parity_sign`, which returns `-1 if int(value).bit_count() & 1`:
the parity of the count's binary representation rather than the count mod 2. The
two agree at 0 and 1 and diverge from 2 upward. The tested fibers reach only
phases 0 and 1, so no reported number is wrong and the fault is latent.

**Result: missed; disposition `accept`.** The review reached the helper and
cleared it. It quotes the phase expression correctly, calls it "the correct outer
character" (review.txt lines 87-90), and never asks what `parity_sign` does with
a count. It verified the expression that builds the phase, not the function that
consumes it.

**What it did instead, recorded and not counted against it.** An independent
gate-replay recomputation of both triples and all 512 fibers (no mismatch in
3072 values), a byte-identical re-run from a repointed copy, a 16384-point WHT
probe, and five execution-time corrections. One is stronger than anything in the
original verdict: the omitted-carry control's zero is parity-forced, because the
wrong reference vanishes identically on every mask containing the h wire and
every carry-sensitive mask contains h, so the control's discrimination comes from
the nonzero true coefficient and not from the zero. It also named two hypotheses
the author had inherited unchecked and verified them itself.

**Reading, and what not to read into it.** One case does not calibrate a class.
What it does show is structural: a referee that establishes independence by
recomputing the reported numbers cannot find a latent fault that way, because
latency is precisely the property of agreeing on the tested inputs -- the
reviewer's own limits call gate replay "the independent leg". This class has to
be caught by reading a helper's domain, which is how the original catch was made,
by the coordinator. The brief change the protocol requires after a failed trial
is now in the "Referee" role duty; repeat this class after that change before
concluding anything about referee quality.

**Usage.** 102,312 tokens, 23 tool calls, 11.6 minutes -- roughly half a live
round's review for a submission a fraction of the size.

---

### Live-round addendum, 2026-09-21 (CADO cache pilot and recurrence audit; not blinded)

Scientific authorship was delegated; the coordinator integrated the resulting
C120/C121 records. The recurrence author supplied a counterexample to a naive
run shortcut and explicitly stopped the unsupported performance direction.
The cache author implemented the reviewed allocation-lifetime change and kept
reader/builder ablations separate. Claim records own the scientific outcomes.

Author and coordinator independently noticed that the first elapsed readings
were quantized by timeout polling. This was found before the results referee,
so it is not a blinded referee catch. Review Vf4eac155b56c4dd5 independently
confirmed the mechanism in the interpreter source, retained valid child RSS
and exact-output evidence, and prescribed a focused timing amendment. All
completed v1 runs kept their real zero exits and original evidence. The precise
elapsed conclusion was withheld until the corrected same-workload run.

A fresh focused context checked the revision against raw records, fixed input
and runtime identities, the complete prescribed sequence and all reported
summary arithmetic. The slower combined-arm observation was retained, and
prior setup CPU was not converted into a certified total by an administrative
allowance. Scope remains CADO-child measurements on synthetic matrices, with
shared diagnostics/backend limitations disclosed. These live reviews are
substantive checks, not evidence that reviewer defect-detection ability has
been calibrated by blinded trials. Note FE owns the full evidence sequence.

## 6. Findings: arb 0.1.0 (for `~/Workspace/agent-research-board`)

- **A1. The `not_initialized` hint invites a second board.** It is returned
  anywhere no `arb.toml` is found, including a worktree of an initialized
  project. Suggest: mention `--project` first, and warn that `init` creates a
  separate board.
- **A2. Closed write submissions become unreadable.** `submission.get` on
  `S379136f790ac46f0` (closed, accepted) now fails with `candidate_changed`,
  because `todo/open/50-…` was legitimately edited afterwards. That blocks
  auditing closed work through the API. Suggest: return the frozen manifest
  with a `candidate_changed` flag for reads, and keep the refusal only for
  review and closure.
- **A3. A reviewer is a single actor.** A fresh reviewer per version needs
  `task.assign_reviewer` or reuse of the actor ID. Consider a reviewer role or
  pool (e.g. any actor matching `referee-*`) that still excludes the author.
- **A4. No log-reuse warning.** Suggest `run.finish` warn when a log path in
  the same attempt was already recorded with different content.
- **A5. Actor metadata.** Actors carry no host or model, so reviewer
  calibration can't be analyzed by host or model. Suggest an optional actor
  profile (host, model, role) set at first use.
- **A6. Minor: inconsistent limits.** `task.list` accepts `limit ≤ 100` and
  `task.history` `limit ≤ 50`. The errors are clear; the difference is
  surprising.
- **A7. Schema flags.** `agent.resume` and `inbox.read` are `read:false`
  because they record delivery receipts. That is correct but easy to misread
  as "mutates the cursor". A one-line note in the schema description would
  help.

---

## 7. Where this review may be wrong

- **Referees are unproven.** A fresh referee may catch less than the
  coordinator's deep reading. P1 includes shadow review and the §5 trial for
  that reason. If referees underperform, the right fix may be referee-assisted
  coordinator review, not replacement.
- **Blind phase A (P2) may not pay** for small questions. It trades latency
  for independence.
- **Lint heuristics** were calibrated on 62 submissions from one project.
  Naming conventions (`*failed*` = preserved) are assumptions.
- **Codex runtime capabilities** (per-worker sandboxing, spawn limits) are
  unknown to me. P8 and S1 may already be handled outside the skills.
- **Counts** are from a live, moving board.
- **Unverified on the Claude side:** hooks firing inside a real subagent
  (types now load; a smoke test is pending), and referee quality.

## 8. Reproduce

```bash
cd research
arb call task.list --data '{"limit":100}'                       # task states, reviewers
arb call task.get --data '{"task":"T..."}'                      # spec: kind, dependencies, independent_review
arb call task.history --data '{"task":"T...","full":true}'      # reviews and closures
python3 ~/.claude/skills/qsim-swarm/scripts/evidence_lint.py --submission S...
python3 ~/.claude/skills/qsim-swarm/tests/test_scripts.py       # guard and lint tests, mutation-checked
```

Claude-side artifacts: `~/.claude/skills/qsim-swarm/` (playbook, briefs,
calibration, scripts, tests), `~/.claude/skills/qsim-worker/`,
`~/.claude/skills/qsim-research/`, and `~/.claude/agents/qsim-*.md`.

---

## 9. Codex adoption notes (2026-09-12)

The user authorized these review notes and configuration changes for future
runs from a side conversation. This update changes workflow instructions and
mechanical lint, not the science or the live swarm's assignments. No board
task, inbox, review, acknowledgement or agent was mutated/contacted. The
original counts describe the review snapshot, not the current live board.

### Decisions on process and skills

| Item | Decision and qualification |
|---|---|
| P1 | Adopt fresh, uninvolved referees and `independent_review:true` on coordinator write tasks. Keep coordinator shadow review until calibration supports the relevant defect classes. The flag excludes the author; it cannot establish a reviewer's independence. A provenance-only revision can reuse the prior scientific audit by exact hash with a focused version review. |
| P2 | Adopt the blind reference phase where it matters. Use a fresh context with the actual circuit contract/gates and explicit permitted sources. Excluding a proof from `inputs` alone neither removes inherited knowledge nor prevents reading it. Disclose shared code and accidental exposure; document a justified trivial-reference exception. |
| P3 | Adopt relevant mutation checks and a fixture-to-mutant table, usually 2–4 plausible mutants when warranted. Include asymmetric/nonzero cases where meaningful. Qualify the original wording: a surviving mutant needs investigation or an equivalence explanation; it does not prove every fixture vacuous. |
| P4 | Adopt a corrected shared linter as a mechanical gate; reject the original script's “reusable as-is” assessment. See the defects and compatibility change below. Errors block submission/review; prose limitations alone do not waive them. Lint is not mathematical review. |
| P5 | Adopt `kind:"experiment"` for scientific execution, including reference/mutant runs, and real accepted/closed dependencies between phases. These fields must be set when future tasks are created; adding instructions is not an automatic runner. |
| P6 | Adopt versioned source/report/log paths and archive-based review. Freeze before submission; new material uses new filenames/submissions. Do not rename or edit already submitted artifacts or rewrite earlier run metadata. |
| P7 | Adopt explicit measured / proved bound / estimate labels and checks of cited claim hypotheses. Charge setup, live buffers, coefficient precision and extraction; support/term counts are not allocated bytes. |
| P8 | Retain existing Claude hooks as advisory guardrails. No verified Codex per-worker sandbox is exposed here; no unsupported spawn parameters or isolation guarantees were added. Unit tests do not establish real subagent hook invocation. |
| S1 | Update Luna's sequential-session wording to one coordinator/canonical writer and bounded, authorized parallel workers with distinct actors. The one-Astra default is retained. The skill applies only when Luna coordinates or that division is requested, so parallel Astra-led work was not by itself a violation. |
| S2 | Name both host entry points in `SWARM.md` and route through `AGENTS.md`, `arb.toml` and installed adapters. Leave `CLAUDE.md` intact during this live round because active tasks have snapshotted it. |
| S3 | Adopt one repository-owned swarm protocol, with Codex/Claude skills and role briefs referring to it. `CLAUDE.md` still owns repository conventions; `METHOD.md` still owns scientific method. |
| S4 | Adopt explicit absolute `arb --project` paths and no ordinary git worktrees for this dirty/untracked workspace. A missing board means check the root, not initialize a second board. |
| S5 | Adopt closure-before-input-edits, with an additional dependency caveat below. Closing a producer is insufficient protection for a dependent integration. |

**Dependency caveat:** arb recursively revalidates upstream inputs for consumers,
including inputs of closed tasks. Editing a mutable TODO after closing its
producer can therefore still stale dependent integration. Put administrative
backlog pointers in `source_ref`; when backlog text actually defines the
scientific contract, freeze that content as a stable input. Keep genuine
technical dependencies and resnapshot/revalidate when they change. Do not
evade staleness by omitting real inputs.

### Lint corrections and historical compatibility

The original linter had three reproduced gaps that made a clean exit unsafe
as an acceptance gate:

1. An unreadable submission, including `candidate_changed`, could produce
   only an informational finding and exit successfully without inspecting it.
2. A JSON report supplied as an exit-bound check log was dispatched as generic
   log text, missing an `ok:false` report paired with exit code zero.
3. A filename/header claiming “failed” or “superseded” could downgrade active
   harness defects without provenance proving that the artifact was historical.

The shared implementation in [tools/evidence_lint.py](tools/evidence_lint.py)
now fails closed for inaccessible/invalid archives, checks JSON verdicts and
unresolved declarations, rejects ambiguous check-to-log versions and conflicting
exits, and binds log identity to both source path and hash. AST recognition
avoids treating a unit test's quoted harness example as an executed harness.
Static call presence still does not prove reachability or control relevance.

Preserved failures require explicit source/hash/reason metadata via
`--preserved-manifest`. For submission review that metadata must itself be
frozen evidence. Historical status cannot excuse an exit-zero check citing a
failing result, a missing archive or hash corruption. The format and commands
are maintained in `SWARM.md`.

The original corrected historical submission `Sb702c78972f14f94` also archives
two failed predecessor scripts without this new metadata. It now yields five
errors, all confined to those two predecessors; the corrected artifacts are
not implicated by those errors. The historical regression expectation was
updated to require exactly those locations/classes. This does not rewrite the
archive, retract its accepted result, or claim all historical submissions pass
the stricter gate. New positive/negative preservation cases use temporary fixtures.

### arb implementation proposals

| Item | Disposition |
|---|---|
| A1 | Agree with improving the `not_initialized` hint. Explicit-root instructions mitigate it now; changing arb's error text is deferred to its own repository. |
| A2 | Confirmed `submission.get` on `S379136f790ac46f0` fails with `candidate_changed` after a legitimate later TODO edit. Returning frozen history with a warning while retaining mutation checks is the appropriate fix. arb code is unchanged; the linter reports inaccessible evidence as an error, with no reopening-history workaround. |
| A3 | Defer reviewer pools. Fresh agent instances can use the assigned actor sequentially or the coordinator can reassign it; never share one actor concurrently. A pool alone would not prove independence. |
| A4 | Agree with a board-level log-reuse warning. Shared lint now detects ambiguous versions; changing `run.finish` is deferred to arb. |
| A5 | Useful for calibration, but actor-profile schema changes are deferred. Record available host/model/effort and usage with future calibration results; do not invent unavailable measurements. |
| A6 | Low priority; retain current documented per-operation limits. No schema changes. |
| A7 | Clarified delivery receipts versus consumed cursors in the shared protocol. The schema-description improvement remains an arb change. |

### Installed configuration and validation

Repository changes: this addendum, [SWARM.md](SWARM.md), its entry-point links
in [AGENTS.md](AGENTS.md) and [arb.toml](arb.toml), and the shared linter plus
[its regression tests](tools/test_evidence_lint.py). `arb project.info` confirms
that the project now advertises `SWARM.md` in its instruction list.

Installed adapters updated on this machine: Codex `qsim-research` and
`qsim-luna-coordinator`; Claude `qsim-research`, `qsim-swarm`, `qsim-worker`,
associated board/brief/calibration references, and the falsifier/referee role
bodies. Existing role model settings and guard hooks are retained. The old
Claude lint path forwards to the single repository implementation. These
installed files live outside the repository and must accompany a move to a
different host; the forwarding path is local to this checkout.

Validation performed without scientific execution:

- Shared lint regressions: **18 tests pass**, using temporary artifacts and
  synthetic records, including the three original gaps and log/hash binding.
- Claude guard/forwarder/history regressions: **19 tests pass**, with the
  historical preservation expectation described above. Historical calls are
  read-only; guard tests use disposable projects.
- Skill-creator's validator: **all 5 modified skills pass**. Both changed role
  frontmatters parse; their model/hook configuration is retained. TOML parses
  and every configured repository instruction path exists.
- Documentation gate: **all 10 checks pass**. No science suite, index rewrite,
  claim/paper edit, commit or publication was performed for this update.

These checks validate configuration and mechanical behavior. Referee quality,
real-agent hook invocation, host isolation and cost savings remain unmeasured.
Calibration will require a separately budgeted blinded mix of defective and
corrected evidence; keep its answer key out of referee inputs. The new protocol
applies when future assignments or scientific revisions are opened, without
retroactively restarting the active swarm.

---

## 10. Adoption notes, 2026-09-17 (Claude coordinator, user-authorized)

The user asked for these process changes in the coordinating session and for
this record so that Codex coordinators and workers have the same information.
Section 9 remains the record of 2026-09-12 and is not edited; where it says
preservation metadata must be supplied with `--preserved-manifest`, the change
below supersedes it.

### Shared protocol and method (both hosts)

| Change | Where |
|---|---|
| The loop gains IDEATE (a slate of candidates) before observation and a DESIGN REVIEW before any compute; "test only what the derivation could not settle" | METHOD.md, "The loop, as actually run here" |
| Seven idea generators, a candidate format, and the rule to build a slate after three residual items in a row | METHOD.md, "Where ideas come from" |
| Exploration is not a test: a numerical match is not a mechanism; a cited result must be about the same object; every run states its decision value | METHOD.md, "Exploration is not a test" |
| Failure modes: testing a coincidence as a mechanism; a confound visible from arithmetic alone | METHOD.md failure-mode table |
| Coordinator's own runs follow the same rules as assignments; exploration is labelled and cannot support a claim | SWARM.md, "Choose work that can change the result" |
| Phase 0 divergent slate; phase 2 design review before compute; phases renumbered 0-6 | SWARM.md, "Independence and phase order" |
| Blinding covers examples in permitted files | SWARM.md, phase 3 |
| Gate stops and proportionate amendments | SWARM.md, after the phase order |
| Role duties for every worker, deriver, surveyor, falsifier and referee, and the worker return block, stated host-neutrally | SWARM.md, "Role duties (host-neutral)" |
| A referee who must execute gets its own experiment task; `stale_attempt` on a submitted attempt is correct | SWARM.md, "Scientific review and referee calibration" |
| Re-running an author's script only from a copy with repointed outputs | SWARM.md, same section |
| Board practicalities: inline `--data` JSON, 16,000-character review bodies, 30 evidence entries, subagent `.md` refusal | SWARM.md, "Board fields and evidence lifecycle" |

### Mechanical lint

`tools/evidence_lint.py` now discovers preservation metadata among a
submission's frozen evidence (any JSON object with a `preserved` key) instead of
requiring the referee to name it. It reports each manifest applied as
`PRESERVATION-APPLIED` with the number of artifacts it covers, rejects two
manifests that disagree about one source (`PRESERVATION-CONFLICT`), and treats a
malformed manifest as an error rather than ignoring it. `--preserved-manifest`
remains for a pre-submission file list and, on a submission, pins one manifest
and rejects any other that also preserves (`PRESERVATION-UNNAMED`). The
exit-zero rule is unchanged: preservation never excuses a successful check
citing a failure. `tools/swarm.py review-start` needed no change. Regression
tests went from 19 to 27; five of the new tests fail against the previous
implementation. First live use: `Sc4c99b9dd0c14e6c`, where `review-start`
applied `preserved_v3.json` without a flag.

### Host adapters

- **Claude** (`~/.claude/agents/`, `~/.claude/skills/qsim-swarm/references/briefs.md`):
  deriver gains the same-object check, match-is-not-mechanism, decision value
  and coincidences per prediction, optional uncalibrated Lean, slate
  assignments and the instantiate-every-control rule; surveyor gains slate
  assignments; referee gains design review; the derivation brief includes the
  experiment plan. These bodies now duplicate parts of SWARM.md "Role duties";
  SWARM.md is canonical, and reducing the Claude bodies to pointers is an open
  decision.
- **Codex** (`~/.codex/skills/`): `qsim-luna-coordinator` requires every role
  brief to name the role, its phase, SWARM.md "Role duties" as binding and the
  return block, and forbids assigning an experiment before its design review;
  `qsim-research/references/board.md` points to the role duties, `sweep`,
  `review-start` and the preservation change. Codex has no role files, so the
  brief is the only carrier of a role.

### Validation

Documentation gate 10/10 and indexes current; swarm citation sweep with no
unresolved IDs; lint regressions 27/27; swarm tooling tests 7/7; Claude hook and
forwarder tests pass; skill-creator validation passes for both modified Codex
skills. None of this has been exercised by a Codex-coordinated round, and
referee quality remains uncalibrated.


### Later on 2026-09-17: generated role bodies and a lint false positive

- **Role bodies are generated.** The duplication listed as open under "Host
  adapters" is closed. SWARM.md "Role duties" is expanded to carry everything the
  Claude role bodies said (a 47-phrase coverage check found no gaps) and is
  delimited by `<!-- role:NAME -->` markers. `tools/gen_role_bodies.py` writes
  each `~/.claude/agents/qsim-ROLE.md` body as a short Claude header plus the
  "Every worker" and role blocks, preserving frontmatter byte for byte.
  `tools/check.py` gains check 8b, which fails on a stale body and reports
  "skipped" when the agents directory is absent (a Codex-only host). Edit
  SWARM.md and run the generator; never edit an agent body. The applied bodies
  equal the reviewed draft except for the generator path in their header line.
- **Lint false positive fixed.** `HARNESS-NO-CONTROL` fired on a harness script
  that registered its must-fail controls in a loop (`exp.must_fail(pid, ...)`),
  because only literal ids were counted; the same script's report resolved all
  seven controls. The error now requires the absence of any `must_fail()` call;
  dynamic ids remain reported as `HARNESS-DYNAMIC-IDS`, and the harness report's
  declarations stay the resolution check. Found on `S0077516c27884295`. A new
  test fails against the previous implementation; a second guards that a script
  with no `must_fail()` call still errors. Lint regressions 29/29.
- Validation after both changes: documentation gate 11/11, indexes current,
  generator `--check` exit 0, Claude hook tests pass, citation sweep with no
  unresolved IDs, swarm tooling tests pass.

---

## 11. Adoption notes, 2026-09-17 (second pass: authorship, revisions, calibration)

A fresh Claude session reviewed this record end to end at the user's request and
proposed five improvements. The user took three; the other two (evidence living
outside git, and `HANDOFF.md` having re-grown past a thousand lines) were
deferred deliberately, not rejected. Section 10 stands; this section adds to it.

### Authorship is not a coordinator duty

**Observation.** P1 separated the reviewer from the integrator and left the
author where it was. In every round in section 5 the coordinator wrote most
submissions, and the addenda record it as the source of a large share of the
referees' findings — including both inaccuracies found inside one submission's
own limitations text, a byte ratio quoted from two different steps, and a
print-only validator that exited 0 regardless. The one round that delegated
derivation (2026-09-16/17) produced a worker that refuted two statements in the
coordinator's own brief before any review.

**Change.** SWARM.md, "Choose work that can change the result": the coordinator
owns dispatch, integration and the ledger; scientific authorship is not one of
those. Derivation, implementation and measurement go to a worker whenever one
can be assigned. When the coordinator must author, the task records that, its
shadow review does not occupy the review slot, and the submission goes to the
same fresh referee a worker's would get.

### Revisions and corrections

**Observation.** Version 2 has been this project's most defective version.
C100's v2 added a wrong justification to an already proven claim and a false iff;
TODO 55's v2 was marked as fixing PD2 and repeated its class in a different
control (PA-1); EC-1 caught a ±3 KB precision band tightened into an equality
that array headers make unreachable. Separately, the α + μ threshold entered a
canonical claim verbatim from an *accepted* review and was wrong by one level,
and the fix for it landed in two of the fact's four homes.

**Change.** SWARM.md gains a "Revisions and corrections" section: re-derive
rather than transcribe (an accepted review is not a source); ship with every fix
the check that would have caught the original, or say why none exists; a fix
touching a control, bound, tolerance or precision band re-states and instantiates
that object at every planned point; re-read the neighbourhood, because the
defects revisions introduced here sat beside the correction rather than in it;
correct a fact in its one home, and if it has several, remove one. The revision
lists per correction what changed, what was re-checked and the detector, and the
referee verifies that list against the bytes. The "Every worker" and "Referee"
role duties carry the author-side and review-side halves, so the generated Claude
role bodies carry them too. METHOD.md's failure-mode table gains three named
modes: correction by patch, a referee's formula used as a source, and a tolerance
tightened between versions.

### Calibration: what shadow review is, and what would end it

**Observation.** "Keep shadow review until blinded calibration" had no end
condition and no procedure cheap enough to run, so it had run for six rounds
without a single trial, while the addenda kept identifying new classes to test.
Over those same six rounds shadow review originated no scientific finding; its
real catches were process ones.

**Change.** SWARM.md now calls coordinator shadow review what the evidence says
it is — author-side edit hygiene and process control, which does not satisfy the
independent-review requirement for any submission — and specifies the trial that
would discharge it: packet mode built from the read-only `review-start` output
and the frozen blobs, no board access and no submission or task ID for the
reviewer, defective cases paired with their corrected versions, grading on the
returned reasoning rather than the verdict word, and a stated discharge
condition (the seeded defect returned on every defective packet of a class, no
blocking false alarm on its corrected pair, at least three classes).

### Initial research uses the best model available

The user's decision, recorded in SWARM.md's team paragraph: the phase-0 slate and
the phase-1 derivation run on the best model available, because that is where the
direction is chosen and a candidate nobody generated is invisible for the rest of
the round. A cheaper model stays confined to bounded reference or provenance work
where that has been validated, and is barred from ideation and derivation. This
closes the question of using a cheaper generator for the slate.

### HANDOFF is a mirror; its state block is now generated

**Observation.** At 2026-09-17 `HANDOFF.md` was 1034 lines in 64 checkpoint
blocks reaching back to 2026-08-08. The current checkpoint was 33 lines; the
other 996 were history. It carried 172 claim references over 67 distinct claims,
162 TODO references and 44 distinct board record IDs, each with an owning record,
while the "do not re-run this" warnings worth keeping came to 12 lines. It is a
hand-maintained mirror of the board and the ledger -- the 2026-08-08 failure mode
-- and it regrew because nothing said who owns a round's state.

**Change.** The board owns live assignments, owners and lease freshness
(`CLAUDE.md`'s facts table and `SWARM.md`'s ownership list now say so).
`tools/handoff.py` generates HANDOFF's state block between
`<!-- generated:state -->` markers from the board, git and `todo/open`, including
lease freshness, so a resuming session can see whether another session is
actually running rather than inferring it from processes. `tools/check.py` gains
check 8c, which bounds the hand-written part at 60 lines and the file at 200 and
reports "skipped" until the markers exist. Freshness is deliberately not gated:
the block reports live state, carries its own timestamp and is regenerated at
round boundaries and on resumption. The migration -- auditing the 64 blocks,
moving the surviving warnings to the todos and claims that own them, archiving
the rest -- is `todo/open/61`, to be done between rounds, since a HANDOFF that
drifts mid-review invalidates its own freeze.

**For Codex.** `AGENTS.md`, `CLAUDE.md`'s generated-files table, the Luna
coordinator skill and both hosts' `references/board.md` carry the rule and the
command; Codex changes HANDOFF's state the same way, by rerunning the tool.
First run reported all four open TODO-55 assignments holding leases that expired
6 to 8 hours earlier, with no actor live.

**A constraint found while building the first trial.** The obvious cases are
leaked by the repository's own honesty. `S6e92d733a8074097`'s vacuous fixtures
are described in `notes/DM`; `S52fc7675ddec4ca6`'s latent `outside_sign` defect
is described in `notes/DS`. A referee doing what it is told — read the ledger —
would be handed the answer. Trials therefore leak-check the case first and hold
out `claims/`, `notes/`, the aggregates and `HANDOFF.md` when it is leaked,
which narrows what a trial measures to detection from the construction and the
bytes. Packet mode also cannot hide from the reviewer that it is a trial. Both
limits are now in the protocol rather than discovered again later. The trial
itself was then run: its entry sits with section 5's addenda as "Calibration
trial 1", it missed its seeded class, and the brief change it produced is in the
"Referee" role duty.


## 12. Evidence-gate hardening, 2026-09-21 (Codex, user-authorized)

The user requested a brief review of the recommendations and implementation of
useful suggestions. Most of P1–P7 and S1–S5 already have shared-protocol or tool
implementations. This pass strengthens P4 and P6 at the executable boundary;
it does not start research, delegate work or change historical board records.

### Applied changes and detectors

- **Recorded exit versus log exit.** A quiet command ending with `[exit 3]`
  could be recorded as exit zero without a lint error. `evidence_lint.py` now
  compares the terminal wrapper marker with the recorded check/run exit,
  including signed signal exits. Tests cover both mismatch directions,
  changed working copies versus frozen logs, logs present only in run records,
  and preservation metadata's inability to excuse a false success. Matching
  failures and logs without wrapper markers remain valid.
- **An unreadable lint result was not a successful check.** `review-start`
  previously returned zero for malformed JSON if its lint subprocess returned
  zero; some malformed structures instead crashed the text formatter. It now
  validates the finding structure, error count and process status, returning a
  failure with a diagnostic in text and JSON modes. Tests feed malformed and
  contradictory results as well as clean results and genuine nonzero exits.
- **Freeze protection starts at file creation.** `swarm.py check` previously
  checked existence and then opened a log for truncating write. It now creates
  logs exclusively, refusing an intervening writer or an existing leaf symlink
  before executing the command. The regression creates another writer's log
  during setup and checks that its bytes survive and no command runs. This is
  protection against accidental evidence replacement, not a worker sandbox.

Six new regression methods reproduced defects against the prior code; another
covers valid log cases. After the fixes, evidence-lint tests pass 33/33 and
swarm-helper tests 10/10, including the disposable-board test that read commands
leave board records unchanged. No scientific suites were needed.

### Judgment on the wider recommendations

Keep fresh referees, dependencies that bind the actual accepted inputs, frozen
artifacts, one shared protocol, and mechanical checks that remove routine work
from review. The existing calibration limits are also appropriate: a clean lint
or an agreeable live round does not establish scientific reviewer quality.

Do not adopt the original P3 wording that surviving mutants automatically make
fixtures vacuous; equivalence and unexercised domains need investigation. The
shared protocol already qualifies this. Likewise, adding a reviewer pool does
not establish independence, and mandatory extra phases for trivial references
would spend latency without necessarily improving evidence. Keep the existing
scoped exceptions. Evidence retention outside git and broader referee
calibration remain useful separate work, not prerequisites for these fixes.
