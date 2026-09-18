# Codex entrypoint for the research workspace

Read [CLAUDE.md](CLAUDE.md) for repository conventions and validation rules;
it remains the shared instruction source. Read [METHOD.md](METHOD.md) before
designing or running research, and [HANDOFF.md](HANDOFF.md) when resuming it.
Current claim status and evidence live in the individual `claims/` records.

`HANDOFF.md`'s state block is generated between `<!-- generated:state -->`
markers by `uv run python tools/handoff.py`; never hand-edit inside them, and
refresh it when a round opens or closes and when resuming. Live assignments,
their owners and their lease freshness belong to the board, not to prose:
`arb --project ROOT board`, or the generated block, will tell you whether
another session is actually running. Everything outside the markers is
hand-written and bounded by `tools/check.py`; `CLAUDE.md` lists every generated
file with its tool.

For authorized swarm work, read [SWARM.md](SWARM.md), the shared Codex/Claude
protocol for future assignments, independent references/referees, dependencies
and evidence lint. Installed host skills are adapters to that protocol.

Use the installed `qsim-research` skill for this project's research discussion,
claim review, experiments, and paper refinement. It is a workflow guide, not
a second ledger or an instruction to continue research during a review-only
request. If skill discovery is unavailable, its source is
`$CODEX_HOME/skills/qsim-research/SKILL.md` (default
`~/.codex/skills/qsim-research/SKILL.md`); read it directly when relevant.

Match validation to the requested work as specified in `CLAUDE.md`. Do not
start science suites for a short conceptual answer. Preserve the current dirty
worktree; do not commit or publish unless requested.
