---
id: 65
state: open
title: "Small fixes to tools/handoff.py, check 8c and HANDOFF pointers carried from the TODO 61 reviews"
outcome: ""
claims: []
---
# Carried from the TODO 61 migration (reviews Vdb88f55fd7bc4180, Vde45154df0a04791)

None changes what HANDOFF says or what check 8c measures today. Do them together at
the next edit of `tools/handoff.py`, with one reviewed submission.

1. **Negative count (required at the next edit).** The truncation line prints
   `+{len(open_tasks) - MAX_TASKS}` whenever board.read sets `tasks_truncated`, so
   3 open tasks with the flag set print "+-9". Unreachable with `limit=100` and
   `MAX_TASKS = 12`. Fix with `max(..., 0)`. Detector: the v1 referee's probe
   (`out/agent-board/reviews/referee-handoff-61/probe_board_lines.py`), whose
   3-open-plus-truncated shape must show no minus sign.
2. The non-truncated overflow line lost the word "more".
3. `tools/check.py` 8c still has a `"skipped"` branch that can no longer fire, and
   `todo/done/61` describes the superseded skip behaviour.
4. Reversed markers pass `check()` on short text (the budget catches it on the
   real file), and the write path then reports "rewrote" without changing anything.
5. HANDOFF: its standing-instruction line on Paper A is one 143-character line;
   the compilation-thesis pointer should also name F4 (the propagator bug) and F6.
6. `README.md` lines 72-74 still say HANDOFF states "how to get running".
