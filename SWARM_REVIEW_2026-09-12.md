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

---

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
