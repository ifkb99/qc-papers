# Codex skill sync, 2026-09-14

Written by the Claude coordinator for the Codex `qsim-research` and
`qsim-luna-coordinator` skills. It lists what changed in the shared
protocol and tooling on 2026-09-13/14, and the edits the Codex skills need to
match. `SWARM.md` still owns the protocol; the skills should point to it, not
copy it. Board records: tasks `T28cea8e4305f44e8`, `T5721a746cbf04d9d` and
`T8c4636167d8c4e03`, each closed after independent review.

## What changed

**arb (`~/Workspace/agent-research-board`, `docs/AGENT_GUIDE.md`)**

- Global flags (`--actor`, `--project`, `--request-id`, `--pretty`) may come
  before or after the subcommand. Flag abbreviations such as `--act` are
  rejected.
- Simple calls take fields instead of JSON:
  - `key=value`, typed per the schema.
  - `key:=json`.
  - `key=@file` for text and `key:=@file` for JSON, e.g.
    `summary=@report.md`, `checks:=@checks.json`.
  - A call that needs fields but gets none fails with `missing_input` instead
    of waiting on stdin. Use `--file -` for piped JSON.
- Finding earlier work:
  - `arb search TERMS` searches tasks, messages, submissions, reviews, closures
    and runs.
  - `arb show ID` resolves any cited board ID.
  - `record.check` checks many cited IDs at once.
  - `task.list source_ref=...` lists every task for one TODO.
- `board.read` and the coordinator's `resume` include `attention`: stalled
  leases, blocked or rejected tasks, ready tasks with a cancelled dependency,
  submitted or accepted tasks whose inputs changed, and work waiting longer than
  24 hours.
- `task.close` accepts `records` (the project records the work created or
  changed, hashed into the closure) or `records_note` (why none changed). This
  project's `arb.toml` sets `[closure] require_records = true`, so one of them
  is required.
- `submission.get` now returns a write-task submission's frozen manifest even
  when its canonical files have changed since, pending or closed, and lists
  those files in `candidate_changed`. Before, it refused. Accepting a review
  and closing still refuse changed candidates.

**Project tooling**

- `tools/swarm.py` (tests: `tools/test_swarm.py`) only reads the board or
  formats checks:
  - `sweep`: round-start attention, open questions, held resources, and board
    IDs cited in hand-written records that the board cannot resolve.
  - `todo KEY --records`: every task for one TODO, with its closure records.
  - `review-start S_ID`: the referee packet — task acceptance, frozen manifest,
    and grouped evidence lint.
  - `check [--append FILE] LOG -- CMD ...` and `doc-gate [--append FILE] OUTDIR`:
    run commands into logs and print the `{command,exit_code,log}` objects that
    `checks` expects.
- `tools/evidence_lint.py` warns `CANDIDATE-CHANGED` when a submitted canonical
  file differs from its frozen version.

**`SWARM.md`**

The relevant sections now cover:

- the round-start sweep, and clearing attention before new work;
- searching earlier work before creating a task;
- a note that cites a board ID must still state the finding itself;
- the check and doc-gate helpers, and `review-start` for referees;
- the referee review directory `out/agent-board/reviews/NAME/`;
- integration evidence that names hand-written files and gate logs, not the
  generated `CLAIMS.md`, `NOTES.md` or `INDEX.md`;
- the closure records rule.

Live referee observations go in a dated addendum to `SWARM_REVIEW_2026-09-12.md`
section 5 (calibration material, not a referee input), not in `SWARM.md`.

## Edits for the Codex skills

1. **`~/.codex/skills/qsim-research/references/board.md`**
   - Replace "global flags precede the subcommand" with "global flags may
     precede or follow the subcommand".
   - Add one pointer: before new work, look for earlier work with
     `arb search TERMS`, `arb show ID` and
     `uv run python tools/swarm.py todo KEY --records`; the coordinator starts
     a round with `tools/swarm.py sweep` (see `SWARM.md`).
   - This matches the Claude copy, `~/.claude/skills/qsim-research/references/board.md`.
2. **`~/.codex/skills/qsim-luna-coordinator/SKILL.md`**
   - In "Integrate without repeating the work": coordinators run
     `tools/swarm.py sweep` at round start and clear its attention items, and
     close tasks with `records` or `records_note` as `SWARM.md` requires.
   - In the brief list: a referee brief names `tools/swarm.py review-start S_ID`
     and a review directory under `out/agent-board/reviews/`.
   - Point to `SWARM.md` for details rather than repeating them.
3. **`~/.codex/skills/qsim-research/SKILL.md`**
   - Optional one clause beside "Check relevant retractions before reviving a
     direction": also check board history for that TODO
     (`tools/swarm.py todo KEY --records`).

No Codex edit is needed for the following:

- The Claude-only `worker_guard.py` hook.
  - It now examines commands behind shell keywords, wrappers, `sh -c`,
    substitutions and continuations, and refuses when a command's structure is
    ambiguous.
  - If Codex ever adds similar command checking, keep that design: attempts
    to parse shell more precisely kept opening new bypasses over 11 review
    rounds.
- The arb web UI.

## Validation after editing

- Run the skill validator Codex normally uses on the two changed skills.
- Run `uv run python tools/check.py` from the research root; it should report
  all 10 checks passing.
- Confirm that `arb --project ROOT schema task.close` shows `records` and
  `records_note`.
- These are process edits; no science suite or experiment is needed.
