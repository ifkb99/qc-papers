# QSim swarm protocol

This is the shared protocol for future QSim swarm assignments on Codex and
Claude Code. CLAUDE.md owns repository conventions and validation; METHOD.md
owns scientific method. Claims, notes and TODOs retain their existing homes.
SWARM_REVIEW_2026-09-12.md records the review and adoption decisions.

These rules do not restart, reassign or retroactively relint live/frozen work.
Apply them when opening the next task or scientific revision. User intent and
existing authorization govern scope; loading this file grants no permission
to delegate, execute science, publish or commit. A side conversation remains
an observer unless the user explicitly authorizes its particular mutations.
Do not take over the main coordinator or contact its agents from a side chat.

## Entrypoints and ownership

- Codex loads `~/.codex/skills/qsim-research`; the Luna coordinator adapter
  applies only when Luna is coordinating or the user requests that division.
- Claude loads `~/.claude/skills/qsim-research`, with `qsim-swarm` and
  `qsim-worker` supplying host mechanics. Both hosts read this file instead
  of maintaining separate copies of the protocol.
- One coordinating session owns the `coordinator` actor and canonical writes.
  Parallel workers use distinct actor IDs and write only in their assigned
  attempt output directories. Reviewer notes/mutants go in an assigned review
  directory `out/agent-board/reviews/NAME/`, never into submitted evidence. Actor identity is cooperative,
  not an authentication or filesystem boundary.
- Pass `arb --project /absolute/research/root --actor ACTOR ...` explicitly.
  Use only an actor assigned to this session. Observers use read operations,
  not another session's resume/inbox/ack. Resume and inbox create delivery
  receipts; only ack advances that actor's consumed cursor.
- Do not put current QSim workers in ordinary git worktrees: this dirty
  workspace contains untracked claims/config and ignored board artifacts.
  `not_initialized` means check the supplied root; do not create a second board.
- Scope instructions and Claude hooks are guardrails. Do not claim a per-worker
  sandbox without verified host enforcement. Never add fictitious sandbox
  fields to spawning APIs. Respect actual process handles when reclaiming;
  a lease or stale state file alone does not establish that work stopped.

## Choose work that can change the result

A conceptual answer or process-only edit needs no research swarm or science
suite. For scientific work, state the question, exact output/accuracy, actual
circuit and dirty/clean inputs, baseline, resource cap and what each result
would change. Reuse reviewed evidence by exact identity when it still applies.
Do not schedule a pilot solely because a prior check passed.

For cost or invariance claims, distinguish the full-space permutation from the
clean logical code, pre-QFT expectations from post-QFT probability/sampling,
final support from intermediate peak, and exact statements from floating-point
diagnostics or approximation. State the scaled parameter and what else varies
with the circuit. Small scalar tensor rank is not joint-state rank or a memory
lower bound. Equal logical behavior need not give equal full-space spectra.
Before reordering arithmetic, check commutation and last use on the actual
reachable subspace. For sampling baselines, charge order discovery, orbit
enumeration, transition construction, preprocessing and precision. A known
reduction may settle one approach without closing the broader research question.

Use the smallest useful team. Keep a coherent author for sustained derivation;
reserve capacity for independent reference and referee work instead of filling
every slot with builders. Honour user-selected models. When model choice is
delegated, use stronger reasoning for proofs and substantive scientific review,
and a less expensive model for bounded references/provenance where validated.
Judge savings by completed-task quality, aggregate usage, latency and rework.

## Independence and phase order

1. **Contract and derivation.** Freeze the mathematical object and candidate
   predictions before their test. Give each assignment a concrete acceptance
   condition, sufficient exact technical inputs and a budget.
2. **Blind reference, where independence matters.** Start a fresh context with
   the circuit contract and actual gates, before revealing the candidate proof,
   expected answers or implementation. Excluding a derivation from `inputs`
   alone does not remove inherited knowledge or prevent reading it. Give an
   explicit permitted-source brief; disclose accidental exposure/shared code.
   List primitive gates exercised. A trivial existing reference can skip this
   split with a stated reason and its remaining independence limits.
3. **Reference validation.** Freeze the reference before testing the candidate.
   Demonstrate that relevant fixtures detect plausible faults, including faults
   in the verifier. Use roughly 2–4 named mutants when warranted: incoming
   carry omission, chronology error, wrong sign, count-parity error. Record a
   fixture-to-mutant table and nonzero/asymmetric cases where appropriate.
   A surviving mutant needs investigation or an equivalence explanation;
   survival alone does not prove that every fixture is vacuous.
4. **Candidate experiment.** Dependencies bind the accepted derivation and
   reference versions before the task can be claimed. No result-driven
   adjustment of the reference without preserving the earlier version and
   registering the new prediction. Stop or redirect if the premise is refuted
   or an already known reduction settles the question.
5. **Referee, then integration.** A fresh referee who did not produce the
   scientific inputs judges the exact submission. The coordinator integrates
   only supported results and checks that the wording retains their scope.

Blindness is a workflow property, not enforced filesystem isolation. A shared
bug in the gate interpreter remains possible; disclose which components are
shared and use an independent route at the layer being claimed.

## Board fields and evidence lifecycle

Discover one operation at a time with `arb --project ROOT schema OPERATION`.
Global flags may come before or after the subcommand. Mutations use an explicit
`--request-id`; retry a lost response with that ID and identical arguments, and
use a new ID for changed arguments. Simple calls use fields
(`arb call task.get task=T_ID`; `key:=json` for lists/objects); multiline
content goes through quoted stdin or `--data`.

Find earlier work before creating a task: `arb search TERMS` (tasks, messages,
submissions, reviews, closures, runs), `arb show ID` for any cited board ID,
and `uv run python tools/swarm.py todo KEY --records` for every assignment
against one TODO. `tools/swarm.py` wraps the repeated steps below; it only
reads the board and formats checks, and never creates, reviews or closes.
The installed CLI's agent guide is at
`/home/djneko/Workspace/agent-research-board/docs/AGENT_GUIDE.md`; its design
documents may describe features that the running CLI does not implement.

At round start the coordinator runs `uv run python tools/swarm.py sweep`. It
lists board attention (stalled leases, blocked tasks, submitted or accepted
work whose inputs changed, and work waiting longer than the `arb.toml`
threshold) and board IDs cited in hand-written records that this board cannot
resolve. Clear attention before opening new work: review, reclaim after
confirming the worker stopped, or cancel with a reason. Board state is local
and not versioned, so a note citing a board ID must still state the finding
itself; an unresolved citation is a prompt to restate it, not to restore it.

An assigned worker calls `task.claim`, retains its attempt ID, lease,
`output_dir` and input hashes, and checks for changed inputs before using them.
Renew with `task.heartbeat` before the lease can expire, especially before long
commands; a running process does not renew a lease. Use `task.block` for a
concrete blocker and preserve partial output. Read the assigned actor's inbox
and task thread at handoff boundaries; acknowledge consumed pages and continue
while `has_more`. Put durable coordination in `task:T_ID` or `topic:NAME` and
resolve answered board questions explicitly. Native host messages can point to
board IDs. Expired attempts cannot renew or submit; only the owning coordinator
may reclaim after verifying that the worker stopped. Reopening accepted work
makes its consumers stale. These operations require the session's authorized
role; an observer does not perform them.

Use `arb schema OPERATION` for current field shapes. Research execution tasks,
including reference/mutant runs, use `kind:"experiment"`; pure derivation uses
`kind:"research"`, engineering-only changes `kind:"implementation"`. Mode
still controls writes: workers `read`/`proposal`, coordinator `write` with exact
`write_paths`. Every coordinator write task sets `independent_review:true` and
names a referee. That flag only excludes the submission author; it does not
prove that the reviewer was uninvolved or used a fresh context.

Declare real dependencies, for example
`{"task":"T_REFERENCE","require":"accepted"}` and
`{"task":"T_DERIVATION","require":"accepted"}` for a candidate test.
Use `require:"closed"` when integrated code is needed. A run record proves
recorded ordering and terminal status, not mathematical validity.

Keep changing backlog text in `source_ref` unless its content is a real
scientific input; if needed, freeze that content in a task-owned contract.
Do not omit genuine inputs to evade staleness. Close read/proposal tasks before
editing their inputs. Closing them does NOT protect dependent integration from
later staleness: downstream tasks still validate upstream inputs recursively.
Plan separate stable scientific dependencies and administrative write scope;
if an actual technical input changes, resnapshot/revalidate rather than bypassing
the check. No dependent task may count an unreviewed replacement as the old result.

Core first for scientific execution as required by CLAUDE.md. Post prediction,
call `run.start`, execute within the cap, observe the actual terminal exit,
then `run.finish`. Every run has its own source revision, log and report path,
such as `reference_v2.py`, `run_v2.log`, `report_v2.json`. Finish all started
runs and retain failed runs. Reserve/release GPU resources only when applicable.

Before submission, run the shared mechanical check from the research root:

```bash
uv run python tools/evidence_lint.py --project /absolute/research/root SOURCE.py REPORT.json
uv run python tools/evidence_lint.py --project /absolute/research/root RUN.log --exit-code 0
```

Use the actual exit code. `uv run python tools/swarm.py check LOG -- COMMAND`
runs a command into a project log and prints the `{command,exit_code,log}`
object that `checks` expects; `tools/swarm.py doc-gate OUTDIR` does this for the
documentation gate (`--reindex` regenerates indexes, coordinator only). With
`--append FILE` both collect a list to pass as `checks:=@FILE`. Include
all source, report and log paths used by checks/runs. Submit once the complete version is ready; after submission do
not edit or rename its files. Revisions and new material get new filenames and
a new reviewed submission. Never rewrite an old run's metadata to match a
replacement. Record hashes, physical normalization and stated limitations.

At review start, inspect the frozen manifest and run the evidence lint:

```bash
uv run python tools/swarm.py --project /absolute/research/root review-start S_ID
```

This prints task acceptance, reviewer and candidate status, the frozen
manifest with archive paths, and grouped `tools/evidence_lint.py --submission`
findings (`--verbose` lists every line). It replaces the separate lint call.

Errors block mechanical acceptance until corrected or explicitly established
as historical with the provenance mechanism below. Warnings need a disposition.
Unreadable archives are an error, even for `candidate_changed`; do not call
that a clean lint or reopen historical work just to recover a read. A closed
write task whose canonical files later changed stays readable: `submission.get`
returns its frozen manifest and lists those files in `candidate_changed`, which
the evidence lint reports as a `CANDIDATE-CHANGED` warning.
Lint examines syntax, declarations, reports and provenance. It does not prove
that `finish()` is reachable or that a mutation/control or theorem is meaningful.

For deliberately preserved failures, supply `--preserved-manifest roles.json`
alongside the file list or submission. The metadata format is

```json
{"schema_version":1,"preserved":[{"source":"out/PATH/reference_v1.py","sha256":"FULL_SHA256","reason":"Superseded by v2; defect and replacement are documented in the report"}]}
```

The linter verifies the exact source/hash identity. For a frozen submission,
`roles.json` must itself be frozen evidence. Filename/header conventions never
establish historical status. Historical status cannot excuse an exit-zero check
citing a failing result or an invalid archive. The referee still checks whether
the artifact is historical rather than current scientific support. Existing
archives remain intact; stricter lint does not retroactively retract claims.

## Scientific review and referee calibration

Review archived source/log/report bytes; compare current candidates with those
bytes before integration. Check the actual object, normalization, proof steps,
boundary cases, independent reference, mutant results, provenance and scope.
For each cited claim, check its hypotheses against this case. Label numbers as
measured (name wall time/RSS/tracemalloc/other quantity), proved bounds or
estimates. Count setup, input, live buffers, coefficient bits and extraction.
State-vector/Pauli term counts are not allocated bytes. Novelty requires a
scoped primary-source audit and a fair same-output baseline.

Use fresh referees for initial submissions and substantive scientific revisions.
An unchanged scientific body with a provenance-only addition may reuse the
previous scientific audit by hash and receive a focused version review; explain
that scope. Reusing a reviewer actor for sequential versions does not reuse the
agent's context. Never share one actor concurrently.

Keep coordinator shadow review for initial rounds and uncalibrated defect
classes. Before replacing it, use blinded historical trials containing both
defective and corrected submissions from the review's calibration pointers.
Do not give the trial reviewer the answer key, old verdict or defect-specific
hint. Record its evidence and compare misses, false alarms, time and usage;
do not mutate/reopen historical tasks. Calibration is separately budgeted work,
not automatically executed by a configuration edit. No referee capability or
host isolation has been validated merely by installing these instructions.

Live coordinated rounds do not calibrate referees: they cannot count the
misses nobody caught, and they never substitute for blinded trials. The
coordinator may log a round's referee and shadow-review catches as a dated
addendum to section 5 of `SWARM_REVIEW_2026-09-12.md`, which is calibration
material and not a referee input.

After a supported referee verdict, the coordinator integrates, runs the checks
appropriate to the actual change, submits the write result for independent
review and closes the task with exact files/checks. Integration evidence names
the hand-written files and the gate logs; generated `CLAIMS.md`, `NOTES.md` and
`INDEX.md` files are covered by the reindex check and need not be archived.

`arb.toml` requires every closure to name the records it created or changed
(`records`, matched against its `[records]` patterns and hashed into the
closure) or to say in `records_note` why none changed, for example "negative
audit; finding restated in DW" or "integrated by T_ID". A derivation closed
before its integration cites that integration task; a finding that belongs in
the ledger is not closed with only a board thread as its home. Reuse unchanged validated
work; repeat tests only for changes, failures or unresolved concerns. No commit
or publication follows without the user's instruction.
