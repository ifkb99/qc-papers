---
id: 61
state: done
title: "Migrate HANDOFF.md to a generated state block plus a bounded hand-written intent section"
outcome: "HANDOFF.md migrated 2026-09-18: 56 hand-written lines (direction, C103 headline caveat, standing instructions, paper decisions, do-not-redo list) plus the generated state block; the 1034-line history archived verbatim in archive/HANDOFF-2026-08-to-09.md; tools/handoff.py now reads open tasks from board.read (task.list's first page missed every task after the hundredth)"
claims: []
---
# HANDOFF is a mirror; make the mirrored part generated

## The measurement

At 2026-09-17 HANDOFF.md was 1034 lines in 64 checkpoint blocks, 2026-08-08 to
2026-09-17. The current checkpoint was 33 lines and the other 996 were history.
It carried 172 claim references over 67 distinct claims, 162 TODO references and
44 distinct board record IDs -- every one of which has an owning record. The
"do not re-run this" warnings, the content worth keeping, were 12 lines.

So HANDOFF is a hand-maintained mirror of the board and the ledger: the failure
mode the 2026-08-08 audit found and the restructure was built to remove. It grew
back because nothing said who owns a round's state, so each round appended.

## Already built (2026-09-17)

- `tools/handoff.py` generates the state block between
  `<!-- generated:state -->` markers from the board (open assignments, owners,
  lease freshness, attention, questions, resources), git (modified/untracked,
  with untracked claims, notes and experiments named) and `todo/open`
  frontmatter. `--print` works today; write mode refuses until the markers exist.
- `tools/check.py` check 8c enforces the budgets: 60 lines outside the block,
  200 in the file. It reports "skipped" until the markers exist, so this
  migration can land without breaking the gate meanwhile.
- Freshness is deliberately not gated. The block reports live state and carries
  its own timestamp; regenerate when a round opens or closes and on resumption.

## What remains: the audit

Between rounds, not during one -- a HANDOFF that drifts mid-review invalidates
its own freeze, which is what blocked `Sa4101c802d464f5d`.

1. Read the 64 blocks. For each, decide whether it holds a live instruction with
   no other home. Expect roughly 12 such lines.
2. Move each survivor to the record that owns it: a "do not re-run X" goes to
   the TODO that owns the question (`todo/open/50` already does this well) or to
   a claim's Limits. Moving, not copying -- the point is one home.
3. Move the rest verbatim to `archive/HANDOFF-2026-08-to-09.md`. Nothing is
   lost; it stops being loaded into every session.
4. Leave a hand-written intent section within budget: the current direction and
   why, and the one or two caveats that change how a headline reads (HANDOFF's
   C103 crossover note is the model). Then place the markers and run the tool.

## Known risk

A generated HANDOFF is only as current as the board. A round run without board
records leaves the block quiet while work happens. That is a change of failure
mode, not its removal -- and it makes board neglect visible instead of hiding it
behind prose, which is the trade this item accepts.
